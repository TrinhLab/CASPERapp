from Bio import SeqIO, Entrez
from bs4 import BeautifulSoup
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from ftplib import FTP
import gzip
import pandas as pd
import shutil
import threading
from multiprocessing import Pool
import os
import sys
import ssl
import GlobalSettings

ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "casper2informatics@gmail.com"


#decompress file function - has to be global to run in sub-processes
def decompress_file(filename):
    with gzip.open(str(filename), 'rb') as f_in:
        with open(str(filename).replace('.gz', ''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(str(filename))
    #QtWidgets.QApplication.processEvents()


#model for filtering columns in ncbi table
class CustomProxyModel(QtCore.QSortFilterProxyModel):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._filters = dict()

    @property
    def filters(self):
        return self._filters

    def setFilter(self, expresion, column):
        if expresion:
            self.filters[column] = expresion
        elif column in self.filters:
            del self.filters[column]
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        for column, expresion in self.filters.items():
            text = self.sourceModel().index(source_row, column, source_parent).data()
            regex = QtCore.QRegExp(
                expresion, QtCore.Qt.CaseInsensitive, QtCore.QRegExp.RegExp
            )
            if regex.indexIn(text) == -1:
                return False
        return True


#model for the data in the ncbi search table
class PandasModel(QtCore.QAbstractTableModel):

    def __init__(self, df=pd.DataFrame(), parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent=parent)
        self._df = df.copy()

    def toDataFrame(self):
        return self._df.copy()

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role != QtCore.Qt.DisplayRole:
            return QtCore.QVariant()

        if orientation == QtCore.Qt.Horizontal:
            try:
                return self._df.columns.tolist()[section]
            except (IndexError, ):
                return QtCore.QVariant()
        elif orientation == QtCore.Qt.Vertical:
            try:
                return self._df.index.tolist()[section]
            except (IndexError, ):
                return QtCore.QVariant()

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.TextAlignmentRole:
            return QtCore.Qt.AlignCenter
        if role != QtCore.Qt.DisplayRole:
            return QtCore.QVariant()
        if not index.isValid():
            return QtCore.QVariant()

        return QtCore.QVariant(str(self._df.iloc[index.row(), index.column()]))

    def setData(self, index, value, role):
        row = self._df.index[index.row()]
        col = self._df.columns[index.column()]
        if hasattr(value, 'toPyObject'):
            # PyQt4 gets a QVariant
            value = value.toPyObject()
        else:
            # PySide gets an unicode
            dtype = self._df[col].dtype
            if dtype != object:
                value = None if value == '' else dtype.type(value)
        self._df.set_value(row, col, value)
        return True

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self._df.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self._df.columns)

    def sort(self, column, order):
        colname = self._df.columns.tolist()[column]
        self.layoutAboutToBeChanged.emit()
        self._df.sort_values(colname, ascending=order, inplace=True)
        self._df.reset_index(inplace=True, drop=True)
        self.layoutChanged.emit()


#ncbi
class NCBI_search_tool(QtWidgets.QWidget):

    def __init__(self):
        super(NCBI_search_tool, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'ncbi.ui', self)
        self.logicalIndex = 0
        self.filters = dict()
        self.organism_line_edit.setText(GlobalSettings.mainWindow.orgChoice.currentText())
        self.download_button.clicked.connect(self.download_files_wrapper)
        self.search_button.clicked.connect(self.query_db)
        self.ncbi_table.verticalHeader().hide()
        self.all_rows.clicked.connect(self.select_all)
        self.ncbi_table.setFocusPolicy(QtCore.Qt.NoFocus)
        self.loading_window = loading_window()

    @QtCore.pyqtSlot()
    def query_db(self):
        #setup table
        self.comboBox = QtWidgets.QComboBox(self)
        self.horizontalHeader = self.ncbi_table.horizontalHeader()
        self.horizontalHeader.sectionClicked.connect(self.on_view_horizontalHeader_sectionClicked)

        #Build Query commands
        retmax = int(self.ret_max_line_edit.text())
        if retmax == "":
            retmax = 1000
        org = self.organism_line_edit.text()
        term = '"' + org + '"[Organism]'
        if self.yes_box.isChecked():
            term += ' AND "Complete Genome"[Assembly Level]'
        if self.infra_name_line_edit.text() != "":
            term += ' AND "' + self.infra_name_line_edit.text() + '"[Infraspecific name]'

        #Search DB for IDs
        handle = Entrez.esearch(db="assembly", retmax=retmax, term=term)
        content = handle.readlines()
        content = "".join(str(content))
        bs_content = BeautifulSoup(content, "lxml")

        #Extract IDs
        idlist = bs_content.find('idlist')
        ids = idlist.find_all('id')
        ids = [i.text for i in ids]

        # Get Details on IDs
        handle = Entrez.esummary(db="assembly", id=','.join(ids))
        content = handle.readlines()
        handle.close()
        content = "".join(str(content))
        bs_content = BeautifulSoup(content, 'lxml')

        #print(bs_content.prettify())

        #Prep Data for Table
        assembly_name = bs_content.find_all('assemblyname')
        genbank_ids = bs_content.find_all('genbank')
        refseq_ids = bs_content.find_all('refseq')
        assembly_status = bs_content.find_all('assemblystatus')
        species_name = bs_content.find_all('speciesname')
        temp_strains = bs_content.find_all('infraspecieslist')
        assembly_name = [i.text for i in assembly_name]
        genbank_ids = [i.text for i in genbank_ids]
        refseq_ids = [i.text for i in refseq_ids]
        assembly_status = [i.text for i in assembly_status]
        species_name = [i.text for i in species_name]
        ids = [int(i) for i in ids]
        strains = []
        for i in range(len(genbank_ids)):
            temp_str = str(temp_strains[i])
            temp_str = BeautifulSoup(temp_str, 'lxml')
            temp_str = temp_str.find('sub_value')
            if temp_str != None:
                strains.append(temp_str.text)
            else:
                strains.append('N/A')

        #Get ftp links
        genbank_links = bs_content.find_all('ftppath_genbank')
        refseq_links = bs_content.find_all('ftppath_refseq')
        genbank_links = [i.text for i in genbank_links]
        refseq_links = [i.text for i in refseq_links]
        self.genbank_ftp_dict = {}
        self.refseq_ftp_dict = {}
        for i in range(len(ids)):
            if genbank_ids[i] == '':
                self.genbank_ftp_dict[ids[i]] = ''
            else:
                self.genbank_ftp_dict[ids[i]] = genbank_links[i] + '/'
            if refseq_ids[i] == '':
                self.refseq_ftp_dict[ids[i]] = ''
            else:
                self.refseq_ftp_dict[ids[i]] = refseq_links[i] + '/'


        #Build dataframe
        self.df = pd.DataFrame({'ID': ids,
                            'Species Name' : species_name,
                            'Strain' : strains,
                           'Assembly Name' : assembly_name,
                           'GenBank assembly accession': genbank_ids,
                           'RefSeq assembly accession': refseq_ids,
                           'Assembly Status': assembly_status})

        #Build table view
        self.df.replace('', 'N/A', inplace=True)
        self.model = PandasModel(self.df)
        self.proxy = CustomProxyModel(self)
        self.proxy.setSourceModel(self.model)
        self.ncbi_table.setModel(self.proxy)
        self.ncbi_table.resizeColumnsToContents()
        self.comboBox.addItems(["{0}".format(col) for col in self.model._df.columns])
        self.ncbi_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.ncbi_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)

    @QtCore.pyqtSlot(int)
    def on_view_horizontalHeader_sectionClicked(self, logicalIndex):
        self.logicalIndex = logicalIndex
        self.menuValues = QtWidgets.QMenu(self)
        self.signalMapper = QtCore.QSignalMapper(self)
        self.comboBox.blockSignals(True)
        self.comboBox.setCurrentIndex(logicalIndex)
        self.comboBox.blockSignals(True)
        valuesUnique = self.model._df.iloc[:, logicalIndex].unique()
        #print(valuesUnique)
        if logicalIndex == 0:
            valuesUnique = ['Sort: 0-9', 'Sort: 9-0']
        elif logicalIndex == 2:
            valuesUnique = ['Exclude N/A', 'Only N/A', 'Sort: A-Z', 'Sort: Z-A']
        elif logicalIndex == 3:
            valuesUnique = ['Sort: A-Z', 'Sort: Z-A']
        elif logicalIndex == 4 or logicalIndex == 5:
            valuesUnique = ['Exclude N/A', 'Only N/A', 'Sort: A-Z', 'Sort: Z-A']

        actionAll = QtWidgets.QAction("All", self)
        actionAll.triggered.connect(self.on_actionAll_triggered)
        self.menuValues.addAction(actionAll)
        self.menuValues.addSeparator()
        for actionNumber, actionName in enumerate(sorted(list(set(valuesUnique)))):
            action = QtWidgets.QAction(actionName, self)
            self.signalMapper.setMapping(action, actionNumber)
            action.triggered.connect(self.signalMapper.map)
            self.menuValues.addAction(action)
        self.signalMapper.mapped.connect(self.on_signalMapper_mapped)
        headerPos = self.ncbi_table.mapToGlobal(self.horizontalHeader.pos())
        posY = headerPos.y() + self.horizontalHeader.height()
        posX = headerPos.x() + self.horizontalHeader.sectionViewportPosition(logicalIndex)

        self.menuValues.exec_(QtCore.QPoint(posX, posY))

    @QtCore.pyqtSlot()
    def on_actionAll_triggered(self):
        filterColumn = self.logicalIndex
        self.proxy.setFilter("", filterColumn)

    @QtCore.pyqtSlot(int)
    def on_signalMapper_mapped(self, i):
        indexes = self.ncbi_table.selectionModel().selectedRows()
        #stringAction = self.signalMapper.mapping(i).text()
        if self.logicalIndex == 0:
            if i == 0:
                self.model.sort(self.logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.model.sort(self.logicalIndex, QtCore.Qt.AscendingOrder)
        elif self.logicalIndex == 3:
            if i == 0:
                self.model.sort(self.logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.model.sort(self.logicalIndex, QtCore.Qt.AscendingOrder)
        elif self.logicalIndex == 2 or self.logicalIndex == 4 or self.logicalIndex == 5:
            if i == 0:
                stringAction = "(?!^N/A$)(^.*$)"
                filterColumn = self.logicalIndex
                self.proxy.setFilter(stringAction, filterColumn)
            elif i == 1:
                stringAction = "N/A"
                filterColumn = self.logicalIndex
                self.proxy.setFilter(stringAction, filterColumn)
            elif i == 2:
                self.model.sort(self.logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.model.sort(self.logicalIndex, QtCore.Qt.AscendingOrder)
        elif self.logicalIndex == 6:
            stringAction = self.signalMapper.mapping(i).text()
            filterColumn = self.logicalIndex
            self.proxy.setFilter(stringAction, filterColumn)

    @QtCore.pyqtSlot()
    def download_files_wrapper(self):
        #self.loading_window.loading_bar.setValue(0)
        #self.loading_window.info_label.setText("Downloading Files")
        self.loading_window.show()
        thread = threading.Thread(target=self.download_files)
        thread.start()

    def download_files(self):
        print("Downloading Files")
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        indexes = self.ncbi_table.selectionModel().selectedRows()
        increment = 50 / len(indexes)
        progress_val = 0
        files = []
        for index in indexes:
            NewIndex = self.ncbi_table.model().index(index.row(), 0)
            id = self.ncbi_table.model().data(NewIndex)
            dirs = []
            if self.genbank_checkbox.isChecked():
                genbank_ftp = self.genbank_ftp_dict[int(id)]
                dirs.append(genbank_ftp)
            if self.refseq_checkbox.isChecked():
                refseq_ftp = self.refseq_ftp_dict[int(id)]
                dirs.append(refseq_ftp)

            for dir in dirs:
                link = str(dir).replace('ftp://ftp.ncbi.nlm.nih.gov', '')
                ftp.cwd(link)
                dir_files = ftp.nlst()
                for file in dir_files:
                    if self.feature_table_checkbox.isChecked():
                        if file.find('feature_table.txt') != -1:
                            with open(file, 'wb') as f:
                                ftp.retrbinary(f"RETR {file}", f.write)
                            files.append(file)
                    if self.gbff_checkbox.isChecked():
                        if file.find('genomic.gbff') != -1:
                            with open(file, 'wb') as f:
                                ftp.retrbinary(f"RETR {file}", f.write)
                            files.append(file)
                    if self.gff_checkbox.isChecked():
                        if file.find('genomic.gff') != -1:
                            with open(file, 'wb') as f:
                                ftp.retrbinary(f"RETR {file}", f.write)
                            files.append(file)
            progress_val += increment
            #self.loading_window.loading_bar.setValue(progress_val)


        print('starting decompression')
        #self.loading_window.info_label.setText("Decompressing Files")
        #self.loading_window.loading_bar.setValue(50)
        p = Pool(3)
        p.map(decompress_file, files)
        self.loading_window.hide()
        print('Finished')

    @QtCore.pyqtSlot()
    def select_all(self):
        for i in range(self.df.shape[0]):
            self.ncbi_table.selectRow(i)


#progress bar gui
class loading_window(QtWidgets.QWidget):
    def __init__(self):
        super(loading_window, self).__init__()
        uic.loadUi(GlobalSettings.appdir + "loading_data_form.ui", self)
        #self.loading_bar.setValue(0)
        self.loading_bar.hide()
        self.setWindowTitle("Downloading/Decompressing Files")
        self.info_label.setText("Downloading/Decompressing Files")
        self.hide()