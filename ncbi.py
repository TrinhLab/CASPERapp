from Bio import SeqIO, Entrez
from bs4 import BeautifulSoup
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from ftplib import FTP
import gzip
import pandas as pd
import shutil
import os
import ssl
import GlobalSettings
import platform
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "casper2informatics@gmail.com"

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
        self.download_button.clicked.connect(self.download_files_wrapper)
        self.search_button.clicked.connect(self.query_db)
        self.ncbi_table.verticalHeader().hide()
        self.all_rows.clicked.connect(self.select_all)
        self.back_button.clicked.connect(self.go_back)
        self.ncbi_table.setFocusPolicy(QtCore.Qt.NoFocus)
        self.progressBar.setValue(0)
        self.rename_window = rename_window()
        self.rename_window.submit_button.clicked.connect(self.submit_rename)
        self.rename_window.go_back.clicked.connect(self.go_back)
        self.df = pd.DataFrame()
        groupbox_style = """
                QGroupBox:title{subcontrol-origin: margin;
                                left: 10px;
                                padding: 0 5px 0 5px;}
                QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                                border-radius: 9px;
                                font: 15pt "Arial";
                                font: bold;
                                margin-top: 10px;}"""
        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1","Step2"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1","Step3"))

    def go_back(self):
        """ Clear table """
        self.df = pd.DataFrame() ###Make empty DF
        self.model = PandasModel(self.df)
        self.proxy = CustomProxyModel(self)
        self.proxy.setSourceModel(self.model)
        self.ncbi_table.setModel(self.proxy)
        self.ncbi_table.verticalHeader().hide()
        """ Clear all line edits """
        self.organism_line_edit.clear()
        self.infra_name_line_edit.clear()
        self.ret_max_line_edit.setText("1000")
        self.infra_name_line_edit.clear()
        """ Reset all checkboxes """
        self.yes_box.setChecked(False)
        self.genbank_checkbox.setChecked(False)
        self.refseq_checkbox.setChecked(False)
        self.gbff_checkbox.setChecked(False)
        self.fna_checkbox.setChecked(False)
        """ Hide window """
        self.close()

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
        self.progressBar.setValue(0)
        if self.df.shape[0] != 0:
            self.download_files()

    def download_files(self):
        indexes = self.ncbi_table.selectionModel().selectedRows()
        if len(indexes) == 0:
            return
        if self.genbank_checkbox.isChecked() == False and self.refseq_checkbox.isChecked() == False:
            return
        self.progressBar.setValue(5)
        QtWidgets.QApplication.processEvents()
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
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
            else:
                refseq_ftp = self.refseq_ftp_dict[int(id)]
                dirs.append(refseq_ftp)

            for dir in dirs:
                link = str(dir).replace('ftp://ftp.ncbi.nlm.nih.gov', '')
                ftp.cwd(link)
                dir_files = ftp.nlst()

                for file in dir_files:
                    if self.gbff_checkbox.isChecked():
                        if file.find('genomic.gbff') != -1:
                            # check OS for output path
                            if platform.system() == "Windows":
                                output_file = GlobalSettings.CSPR_DB + "\\GBFF\\" + file
                            else:
                                output_file = GlobalSettings.CSPR_DB + "/GBFF/" + file

                            with open(output_file, 'wb') as f:
                                ftp.retrbinary(f"RETR {file}", f.write)
                            files.append(output_file)

                    if self.fna_checkbox.isChecked():
                        if file.find('genomic.fna') != -1 and file.find('_cds_') == -1 and file.find('_rna_') == -1:
                            # check OS for output path
                            if platform.system() == "Windows":
                                output_file = GlobalSettings.CSPR_DB + "\\FNA\\" + file
                            else:
                                output_file = GlobalSettings.CSPR_DB + "/FNA/" + file

                            with open(output_file, 'wb') as f:
                                ftp.retrbinary(f"RETR {file}", f.write)
                            files.append(output_file)

            progress_val += increment
            self.progressBar.setValue(progress_val)

        QtWidgets.QApplication.processEvents()
        for file in files:
            self.decompress_file(file)
            progress_val += increment
            self.progressBar.setValue(progress_val)

        self.progressBar.setValue(100)
        QtWidgets.QApplication.processEvents()

        for i in range(len(files)):
            files[i] = files[i].replace('.gz', '')
            if platform.system() == 'Windows':
                files[i] = files[i][files[i].rfind("\\")+1:]
            else:
                files[i] = files[i][files[i].rfind("/") + 1:]

        self.rename_files(files)

    @QtCore.pyqtSlot()
    def select_all(self):
        if self.all_rows.isChecked():
            self.ncbi_table.selectAll()
        else:
            self.ncbi_table.clearSelection()

    # decompress file function
    def decompress_file(self, filename):
        with gzip.open(str(filename), 'rb') as f_in:
            with open(str(filename).replace('.gz', ''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(str(filename))

    def rename_files(self, files):
        ans = QtWidgets.QMessageBox.question(self, "Rename Files", "Would you like to rename the downloaded files?", QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)
        if ans == QtWidgets.QMessageBox.Yes:
            self.rename_window.rename_table.setRowCount(len(files))
            cnt = 0
            for file in files:
                self.rename_window.rename_table.setItem(cnt, 0, QtWidgets.QTableWidgetItem(file))
                self.rename_window.rename_table.setCellWidget(cnt, 1, QtWidgets.QLineEdit())
                cnt += 1

            header = self.rename_window.rename_table.horizontalHeader()
            #header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
            #self.rename_window.rename_table.setColumnWidth(1, self.rename_window.rename_table.columnWidth(0))
            self.rename_window.rename_table.resizeColumnsToContents()
            header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
            header.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
            self.rename_window.resize(self.rename_window.sizeHint())
            self.rename_window.show()
        
    def rename_go_back(self):
        self.rename_window.close()
        
    def submit_rename(self):
        #loop through columns and rename the files
        for row in range(self.rename_window.rename_table.rowCount()):
            orig = str(self.rename_window.rename_table.item(row, 0).text())
            new = str(self.rename_window.rename_table.cellWidget(row, 1).text())
            if new != "":
                if orig.find(".gbff")!= -1:
                    if new.find(".") != -1:
                        new = new[:new.find(".")]
                    new = new + ".gbff"
                elif orig.find(".fna") != -1:
                    if new.find(".") != -1:
                        new = new[:new.find(".")]
                    new = new + ".fna"

                if platform.system() == "Windows":
                    if new.find(".gbff") != -1:
                        os.rename(GlobalSettings.CSPR_DB + "\\GBFF\\" + orig, GlobalSettings.CSPR_DB + "\\GBFF\\" + new)
                    else:
                        os.rename(GlobalSettings.CSPR_DB + "\\FNA\\" + orig, GlobalSettings.CSPR_DB + "\\FNA\\" + new)

                else:
                    #unix cannot have spaces in paths
                    new = new.replace(" ","")
                    if new.find(".gbff") != -1:
                        os.rename(GlobalSettings.CSPR_DB + "/GBFF/" + orig, GlobalSettings.CSPR_DB + "/GBFF/" + new)
                    else:
                        os.rename(GlobalSettings.CSPR_DB + "/FNA/" + orig, GlobalSettings.CSPR_DB + "/FNA/" + new)

        self.rename_window.rename_table.setRowCount(0)
        self.rename_window.close()
        GlobalSettings.mainWindow.fill_annotation_dropdown()
        QtWidgets.QMessageBox.question(self, "Download Completed",
                                       "Successfully downloaded file(s) to " + GlobalSettings.CSPR_DB + "\n\nYou may close this window or download more files.",
                                       QtWidgets.QMessageBox.Ok)


                                       
class rename_window(QtWidgets.QWidget):
    def __init__(self):
        super(rename_window, self).__init__()
        uic.loadUi(GlobalSettings.appdir + "rename_window.ui", self)
        self.setWindowTitle("Rename Files")
        self.rename_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.rename_table.setColumnCount(2)
        self.rename_table.setHorizontalHeaderLabels(['Original Filename', 'New Filename'])
        self.hide()

