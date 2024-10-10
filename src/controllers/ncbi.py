from Bio import Entrez
from bs4 import BeautifulSoup
from PyQt5 import QtWidgets, Qt, QtCore, uic
from ftplib import FTP
import gzip
import pandas as pd
import shutil
import os, time
import ssl
import models.GlobalSettings as GlobalSettings
import platform
import traceback
import math
from utils.ui import show_message, show_error, scale_ui, center_ui

#global logger
logger = GlobalSettings.logger

ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "casper2informatics@gmail.com"

#model for filtering columns in ncbi table
class CustomProxyModel(QtCore.QSortFilterProxyModel):
    def __init__(self, parent=None):
        print("Initializing CustomProxyModel")
        try:
            super().__init__(parent)
            self._filters = dict()
        except Exception as e:
            show_error("Error initializing CustomProxyModel class in ncbi tool.", e)

    @property
    def filters(self):
        try:
            return self._filters
        except Exception as e:
            show_error("Error in filter() in custom proxy model in ncbi tool.", e)

    def setFilter(self, expresion, column):
        try:
            if expresion:
                self.filters[column] = expresion
            elif column in self.filters:
                del self.filters[column]
            self.invalidateFilter()
        except Exception as e:
            show_error("Error in setFilters() in custom proxy model in ncbi tool.", e)

    def filterAcceptsRow(self, source_row, source_parent):
        try:
            for column, expresion in self.filters.items():
                text = self.sourceModel().index(source_row, column, source_parent).data()
                regex = QtCore.QRegExp(
                    expresion, QtCore.Qt.CaseInsensitive, QtCore.QRegExp.RegExp
                )
                if regex.indexIn(text) == -1:
                    return False
            return True

        except Exception as e:
            show_error("Error in filterAcceptsRow() in custom proxy model in ncbi tool.", e)

#model for the data in the ncbi search table
class PandasModel(QtCore.QAbstractTableModel):
    def __init__(self, df=pd.DataFrame(), parent=None):
        try:
            QtCore.QAbstractTableModel.__init__(self, parent=parent)
            self._df = df.copy()
        except Exception as e:
            show_error("Error initializing PandasModel class in ncbi tool.", e)

    def toDataFrame(self):
        try:
            return self._df.copy()

        except Exception as e:
            show_error("Error in toDataFrame() in Pandas Model in ncbi tool.", e)

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        try:
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
        except Exception as e:
            show_error("Error in headerData() in Pandas Model in ncbi tool.", e)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        try:
            if role == QtCore.Qt.TextAlignmentRole:
                return QtCore.Qt.AlignCenter
            if role != QtCore.Qt.DisplayRole:
                return QtCore.QVariant()
            if not index.isValid():
                return QtCore.QVariant()

            return QtCore.QVariant(str(self._df.iloc[index.row(), index.column()]))

        except Exception as e:
            show_error("Error in data() in Pandas Model in ncbi tool.", e)

    def setData(self, index, value, role):
        try:
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
        except Exception as e:
            show_error("Error in setData() in Pandas Model in ncbi tool.", e)

    def rowCount(self, parent=QtCore.QModelIndex()):
        try:
            return len(self._df.index)

        except Exception as e:
            show_error("Error in rowCount() in Pandas Model in ncbi tool.", e)

    def columnCount(self, parent=QtCore.QModelIndex()):
        try:
            return len(self._df.columns)
        except Exception as e:
            show_error("Error in columnCount() in Pandas Model in ncbi tool.", e)

    def sort(self, column, order):
        try:
            colname = self._df.columns.tolist()[column]
            self.layoutAboutToBeChanged.emit()
            self._df.sort_values(colname, ascending=order, inplace=True)
            self._df.reset_index(inplace=True, drop=True)
            self.layoutChanged.emit()
        except Exception as e:
            show_error("Error in sort() in Pandas Model in ncbi tool.", e)

## Taken from StackOverflow: https://stackoverflow.com/questions/57607072/update-pyqt-progress-from-another-thread-running-ftp-download
class DownloadThread(QtCore.QThread):
    """ Overall signals """
    finished = QtCore.pyqtSignal(object)
    started = QtCore.pyqtSignal(object)

    """ Download specific signals """
    data_progress = QtCore.pyqtSignal(object) # This signal emits progress data for the progress bar
    data_size = QtCore.pyqtSignal(object) # This signal emits the size of the file being downloaded
    file_started = QtCore.pyqtSignal(object) # This signal emits when the file starts downloading 
    file_finished = QtCore.pyqtSignal(object) # This signal emits the name of the file downloaded

    def __init__(self,parent,url,id):
        try:
            QtCore.QThread.__init__(self,parent)
            self.id = id
            self.url = url # Initialize URL
            self.ftp = FTP('ftp.ncbi.nlm.nih.gov') # Initialize FTP object
            self.ftp.login()
        except Exception as e:
            show_error("Error initializing Thread class in ncbi tool.", e)

    def run(self):
        try:
            ### Start by making sure the url is valid
            if self.url == "":
                self.finished.emit((self.id, False))
                return
            else:
                self.started.emit(self.id)
                self.ftp.cwd(self.url) # Change to appropriate directory
                dir_files = self.ftp.nlst() # Get list of files in directory
                for file in dir_files: # Loop through every file in the directory
                    if GlobalSettings.mainWindow.ncbi.fna_checkbox.isChecked():
                        if file.find('genomic.fna') != -1 and file.find('_cds_') == -1 and file.find('_rna_') == -1:
                            # check OS for output path
                            if platform.system() == "Windows":
                                output_file = GlobalSettings.CSPR_DB + "\\FNA\\" + file
                            else:
                                output_file = GlobalSettings.CSPR_DB + "/FNA/" + file

                            self.ftp.voidcmd('TYPE I')
                            totalsize = self.ftp.size(file)
                            self.file_started.emit((self.id,'Downloading FNA: ' + str(round(totalsize/1e6,2)) + 'MB...',str(totalsize)))
                            with open(output_file, 'wb') as self.f:
                                self.ftp.retrbinary(f"RETR {file}", self.file_write)

                            self.decompress_file(output_file)
                            self.file_finished.emit((self.id,'FNA Downloaded!',output_file))

                    if GlobalSettings.mainWindow.ncbi.gbff_checkbox.isChecked():
                        if file.find('genomic.gbff') != -1:
                            # check OS for output path
                            if platform.system() == "Windows":
                                output_file = GlobalSettings.CSPR_DB + "\\GBFF\\" + file
                            else:
                                output_file = GlobalSettings.CSPR_DB + "/GBFF/" + file
                            
                            self.ftp.voidcmd('TYPE I')
                            totalsize = self.ftp.size(file)
                            self.file_started.emit((self.id,'Downloading GBFF: ' + str(round(totalsize/1e6,2)) + 'MB...',str(totalsize)))
                            with open(output_file, 'wb') as self.f:
                                self.ftp.retrbinary(f"RETR {file}", self.file_write)
                            self.decompress_file(output_file)
                            self.file_finished.emit((self.id,'GBFF Downloaded!',output_file))
                self.finished.emit((self.id,True))
                self.ftp.quit() # Stop the FTP connection once everything has been downloaded
        except Exception as e:
            show_error("Error downloading file within DownloadThread class.", e)

    def file_write(self, data):
        self.f.write(data) # Write the downloaded data to a file 
        # The other signals increase a progress
        self.data_progress.emit((self.id,str(len(data)))) # Emit a signal updating progress

    def decompress_file(self, filename):
        try:
            block_size = 65536
            with gzip.open(filename, 'rb') as f_in:
                with open(str(filename).replace('.gz', ''), 'wb') as f_out:
                    while True:
                        block = f_in.read(block_size)
                        if not block:
                            break
                        else:
                            f_out.write(block)
            os.remove(str(filename))
        except Exception as e:
            show_error("Error in decompress_file() in ncbi tool.", e)

class NCBI_search_tool(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            print("Initializing NCBI search tool...")
            super(NCBI_search_tool, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/ncbi.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("NCBI Download Tool")
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
            self.rename_window.go_back.clicked.connect(self.rename_go_back)
            self.df = pd.DataFrame()
            groupbox_style = """
                    QGroupBox:title{subcontrol-origin: margin;
                                    left: 10px;
                                    padding: 0 5px 0 5px;}
                    QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                                    border-radius: 9px;
                                    font: bold 14pt 'Arial';
                                    margin-top: 10px;}"""
            self.Step1.setStyleSheet(groupbox_style)
            self.Step2.setStyleSheet(groupbox_style.replace("Step1","Step2"))
            self.Step3.setStyleSheet(groupbox_style.replace("Step1","Step3"))
            self.ncbi_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.ncbi_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            #navigation page
            self.goToPrompt = goToPrompt()
            self.goToPrompt.stay.clicked.connect(self.stay)
            self.goToPrompt.close.clicked.connect(self.close)

            #loading label
            self.loading_window = loading_window()
            
            self.genbank_checkbox.toggled.connect(self.check_genbank)

            self.first_show = True
            scale_ui(self, custom_scale_width=1000, custom_scale_height=750)
            self.fontSize = 12
        except Exception as e:
            show_error("Error initializing NCBI_search_tool class.", e)
    
    def check_genbank(self):
        if self.genbank_checkbox.isChecked():
            show_message(
                fontSize=self.fontSize,
                title="Warning!",
                message="The GenBank collection may contain poorly or partially annotated annotation files. We highly recommend using the RefSeq collection if it is available.",
                icon=QtWidgets.QMessageBox.Icon.Warning
            )
            return
        else:
            pass 

    def go_back(self):
        try:
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
            self.ret_max_line_edit.setText("100")
            self.infra_name_line_edit.clear()
            """ Reset all checkboxes """
            self.yes_box.setChecked(False)
            self.genbank_checkbox.setChecked(False)
            self.refseq_checkbox.setChecked(False)
            self.gbff_checkbox.setChecked(False)
            self.fna_checkbox.setChecked(False)
            """ Hide window """
            self.close()
        except Exception as e:
            show_error("Error in go_back() in ncbi tool.", e)

    @QtCore.pyqtSlot()
    def query_db(self):
        try:
            self.loading_window.loading_bar.setValue(5)
            center_ui(self.loading_window)
            self.loading_window.show()
            QtCore.QCoreApplication.processEvents()

            #setup table
            self.comboBox = QtWidgets.QComboBox(self)
            self.horizontalHeader = self.ncbi_table.horizontalHeader()
            self.horizontalHeader.sectionClicked.connect(self.on_view_horizontalHeader_sectionClicked)

            #Build Query commands
            retmax = int(self.ret_max_line_edit.text())
            if retmax == "":
                retmax = 100
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
            bs_content = BeautifulSoup(content, "html.parser")

            self.loading_window.loading_bar.setValue(20)
            QtCore.QCoreApplication.processEvents()

            #Extract IDs
            idlist = bs_content.find('idlist')
            ids = idlist.find_all('id')
            ids = [i.text for i in ids]

            self.loading_window.loading_bar.setValue(35)
            QtCore.QCoreApplication.processEvents()

            # Get Details on IDs
            handle = Entrez.esummary(db="assembly", id=','.join(ids))
            content = handle.readlines()
            handle.close()
            content = "".join(str(content))
            bs_content = BeautifulSoup(content, "html.parser")

            self.loading_window.loading_bar.setValue(55)
            QtCore.QCoreApplication.processEvents()

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
                temp_str = BeautifulSoup(temp_str, 'html.parser')
                temp_str = temp_str.find('sub_value')
                if temp_str != None:
                    strains.append(temp_str.text)
                else:
                    strains.append('N/A')

            self.loading_window.loading_bar.setValue(65)
            QtCore.QCoreApplication.processEvents()

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

            self.loading_window.loading_bar.setValue(80)
            QtCore.QCoreApplication.processEvents()

            #Build dataframe
            data = [{'ID': ids[i],
                     'Species Name': species_name[i],
                     'Strain': strains[i],
                     'Assembly Name': assembly_name[i],
                     'RefSeq assembly accession': refseq_ids[i],
                     'GenBank assembly accession': genbank_ids[i],
                     'Assembly Status': assembly_status[i]} for i in range(len(ids))]

            self.loading_window.loading_bar.setValue(90)
            QtCore.QCoreApplication.processEvents()

            #Build table view
            self.populate_table(data)

            #close loading gif
            self.loading_window.hide()
            self.loading_window.loading_bar.setValue(0)
            QtCore.QCoreApplication.processEvents()
        except Exception as e:
            show_error("Error in query_db() in ncbi tool.", e)

    def populate_table(self, data):
        try:
            # Create a DataFrame from all the data
            self.df = pd.DataFrame(data)
            
            # Reorder columns to the specified order
            columns_order = [
                'ID', 
                'Species Name', 
                'Strain', 
                'Assembly Name', 
                'RefSeq assembly accession', 
                'GenBank assembly accession', 
                'Assembly Status'
            ]
            self.df = self.df[columns_order]
            
            self.df.replace('', 'N/A', inplace=True)
            self.model = PandasModel(self.df)
            self.proxy = CustomProxyModel(self)
            self.proxy.setSourceModel(self.model)
            self.ncbi_table.setModel(self.proxy)
            self.ncbi_table.resizeColumnsToContents()
            self.comboBox.addItems(["{0}".format(col) for col in self.model._df.columns])
            self.ncbi_table.verticalHeader().hide()
            self.activateWindow()
        except Exception as e:
            show_error("Error in populate_table() in ncbi tool.", e)

    @QtCore.pyqtSlot(int)
    def on_view_horizontalHeader_sectionClicked(self, logicalIndex):
        try:
            self.logicalIndex = logicalIndex
            self.menuValues = QtWidgets.QMenu(self)
            self.signalMapper = QtCore.QSignalMapper(self)
            self.comboBox.blockSignals(True)
            self.comboBox.setCurrentIndex(logicalIndex)
            self.comboBox.blockSignals(True)
            valuesUnique = self.model._df.iloc[:, logicalIndex].unique()
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
        except Exception as e:
            show_error("Error in on_view_horizontalHeader_sectionClicked() in ncbi tool.", e)
            
    @QtCore.pyqtSlot()
    def on_actionAll_triggered(self):
        try:
            filterColumn = self.logicalIndex
            self.proxy.setFilter("", filterColumn)
        except Exception as e:
            show_error("Error in on_actionAll_triggered() in ncbi tool.", e)

    @QtCore.pyqtSlot(int)
    def on_signalMapper_mapped(self, i):
        try:
            indices = self.ncbi_table.selectionModel().selectedRows()
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
                    self.model.sort(self.logicalIndex, QtCore.Qt.AscendingOrder)
                else:
                    self.model.sort(self.logicalIndex, QtCore.Qt.AscendingOrder)
            elif self.logicalIndex == 6:
                stringAction = self.signalMapper.mapping(i).text()
                filterColumn = self.logicalIndex
                self.proxy.setFilter(stringAction, filterColumn)
        except Exception as e:
            show_error("Error in on_signalMapper_mapped() in ncbi tool.", e)

    @QtCore.pyqtSlot()
    def download_files_wrapper(self):
        try:
            self.progressBar.setValue(0)

            #make sure rows are present in table
            if self.df.shape[0] == 0:
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="No Query Results",
                    message="Please run an NCBI query to fill the table with results to choose from!"
                )
                return

            #make sure user has selected at least one row
            indices = self.ncbi_table.selectionModel().selectedRows()
            if len(indices) == 0:
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="No Rows Selected",
                    message="Please select rows from the table!"
                )
                return

            threadCount = QtCore.QThreadPool.globalInstance().maxThreadCount() # Get thread count
            if len(indices) > threadCount:
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Too Many Selections!",
                    message="You only have " + str(threadCount) + " threads avaiable to download with.\n\nPlease select " + str(threadCount) + " or fewer rows."
                )
                return

            #make sure file type is selected
            if self.fna_checkbox.isChecked() == False and self.gbff_checkbox.isChecked() == False:
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="No File Type Selected",
                    message="No file type selected. Please select the file types you want to download!"
                )
                return
            self.download_files()
        except Exception as e:
            show_error("Error in download_files_wrapper() in ncbi tool.", e)

    ### When a new thread is started, spawn a label and progress bar for that thread and add it to the form layout.
    def on_thread_start(self, data):
        id = data # This is the id for the thread
        tmp_lbl = QtWidgets.QLabel()
        tmp_lbl.setText("Download(s) Started...")
        tmp_bar = QtWidgets.QProgressBar()
        self.labels[id] = tmp_lbl # Append thread label to list
        self.progressbars[id] = tmp_bar # Append thread label to list
        self.formLayout.addRow(tmp_lbl,tmp_bar) # Add label and bar to form layout

    ### When a thread is finished, update its label and progressbar for the last time
    def on_thread_finish(self,data):
        id = data[0]
        my_bool = data[1]
        if my_bool: # If thread finished succesfully
            self.progressbars[id].setValue(int(self.progressbars[id].maximum())) #Make sure progress bar is full
            self.labels[id].setText("Download(s) Complete!") #Make sure progress bar is full
            self.progressBar.setValue(int(self.progressBar.value()+1)) # Increment overall progress bar when a thread finishes
            QtWidgets.QApplication.processEvents() # Allow the progress bar to update
        else:
            self.error_occurred = True
            if not self.error_shown:
                self.error_shown = True
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Link Failed!",
                    message="Failed to find a valid link for one or more IDs. Please make sure the selected IDs are available in the selected database."
                )
            return

    ### When a file download starts, update the label and progress bar
    def on_file_start(self, data):
        id = data[0] # This is the id for the thread
        prompt = data[1] # This is the prompt to set the label to
        maxval = int(data[2]) # This is the file size
        self.progressbars[id].setValue(0) # Set max value of progressbar
        self.progressbars[id].setMaximum(int(maxval/1e3)) # Set max value of progressbar

        self.labels[id].setText(str(prompt)) # Set label text 

    ## Taken from StackOverflow: https://stackoverflow.com/questions/57607072/update-pyqt-progress-from-another-thread-running-ftp-download
    def on_progress_ready(self, data):
        id = data[0] # This is the id for the thread
        val = int(data[1]) # This is the increment
        self.progressbars[id].setValue(int(self.progressbars[id].value() + val/1e3)) # Increment progress bar

    ### When a file download finishes, update the label and progress bar and add file name to list 
    def on_file_finish(self,data):
        id = data[0] # This is the id for the thread
        prompt = data[1] # This is the prompt to set the label to
        file_name = data[2] # This is the filename
        self.progressbars[id].setValue(int(self.progressbars[id].maximum())) # Make sure progress bar is set to max after finishing download
        self.labels[id].setText(str(prompt))
        self.files.append(str(file_name)) # Add filename to list of downloaded files

    def clean_bars(self):
        for key in self.progressbars:
            label = self.formLayout.labelForField(self.progressbars[key])
            if label is not None:
                label.deleteLater()
            self.progressbars[key].deleteLater()

        self.labels.clear()
        self.progressbars.clear()
        self.threads.clear()
    def clear_layout(self):
        for i in reversed(range(self.formLayout.count())): 
            self.formLayout.itemAt(i).widget().deleteLater()

    ## Taken from StackOverflow: https://stackoverflow.com/questions/57607072/update-pyqt-progress-from-another-thread-running-ftp-download
    def download_files(self):
        try:
            self.labels = {}
            self.progressbars = {}
            self.threads = {}
            self.error_occurred = False
            self.error_shown = False
            indices = self.ncbi_table.selectionModel().selectedRows()
            len_ind = len(indices)
            self.progressBar.setMaximum(int(len_ind))
            if len_ind == 0:
                return
            if not self.genbank_checkbox.isChecked() and not self.refseq_checkbox.isChecked():
                return
            self.progressLabel.setText("Download(s) Started...")
            self.files = []

            for index in indices:
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
                    url = str(dir).replace('ftp://ftp.ncbi.nlm.nih.gov', '')
                    downloader = DownloadThread(parent=self, url=url, id=id)
                    downloader.started.connect(self.on_thread_start)
                    downloader.finished.connect(self.on_thread_finish)
                    downloader.file_started.connect(self.on_file_start)
                    downloader.file_finished.connect(self.on_file_finish)
                    downloader.data_progress.connect(self.on_progress_ready)
                    self.threads[url] = downloader
                    downloader.start()

            while len([self.threads[t] for t in self.threads if self.threads[t].isRunning()]) > 0:
                time.sleep(0.1)
                QtWidgets.QApplication.processEvents()

            if not self.error_occurred and len(self.files) > 0:
                self.progressBar.setValue(int(self.progressBar.maximum()))
                self.progressLabel.setText("Download(s) Complete!")
            else:
                self.progressBar.setValue(0)
                self.progressLabel.setText("")

            self.clean_bars()

            for i in range(len(self.files)):
                self.files[i] = self.files[i].replace('.gz', '')
                if platform.system() == 'Windows':
                    self.files[i] = self.files[i][self.files[i].rfind("\\")+1:]
                else:
                    self.files[i] = self.files[i][self.files[i].rfind("/") + 1:]

            if len(self.files) > 0:
                self.rename_files(self.files)
            else:
                self.clear_layout()
                self.clean_bars()
                if not self.error_shown:
                    show_message(
                        fontSize=self.fontSize,
                        icon=QtWidgets.QMessageBox.Icon.Critical,
                        title="No Files Downloaded",
                        message="No files were downloaded from the selected NCBI files. Please make sure the selected files are available in the database selected."
                    )
                return

        except Exception as e:
            show_error("Error in download_files() in ncbi tool.", e)

    @QtCore.pyqtSlot()
    def select_all(self):
        try:
            if self.all_rows.isChecked():
                self.ncbi_table.selectAll()
            else:
                self.ncbi_table.clearSelection()
        except Exception as e:
            show_error("Error in select_all() in ncbi tool.", e)

    def rename_files(self, files):
        try:
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
            msgBox.setWindowTitle("Rename Files")
            msgBox.setText("Would you like to rename the downloaded files?")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
            msgBox.exec()

            if msgBox.result() == QtWidgets.QMessageBox.Yes:
                self.rename_window.rename_table.setRowCount(len(files))
                cnt = 0
                for file in files:
                    item = QtWidgets.QTableWidgetItem(file)
                    item.setFlags(QtCore.Qt.ItemIsEnabled)
                    self.rename_window.rename_table.setItem(cnt, 0, item)

                    self.rename_window.rename_table.setCellWidget(cnt, 1, QtWidgets.QLineEdit())
                    cnt += 1

                header = self.rename_window.rename_table.horizontalHeader()
                self.rename_window.rename_table.resizeColumnsToContents()
                header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
                header.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
                center_ui(self.rename_window)
                self.rename_window.show()
                self.rename_window.activateWindow()
            else:
                GlobalSettings.mainWindow.fill_annotation_dropdown()
                self.goToPrompt.show()
                self.goToPrompt.activateWindow()
        except Exception as e:
            show_error("Error in rename_files() in ncbi tool.", e)

    def rename_go_back(self):
        try:
            self.rename_window.close()
        except Exception as e:
            show_error("Error in rename_go_back() in ncbi tool.", e)

    def submit_rename(self):
        try:
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
                            if os.path.isfile(GlobalSettings.CSPR_DB + "\\GBFF\\" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "\\GBFF\\" + orig, GlobalSettings.CSPR_DB + "\\GBFF\\" + new)
                            else:
                                show_message(
                                    fontSize = self.fontSize,
                                    icon = QtWidgets.QMessageBox.Icon.Critical,
                                    title = "Renaming Error",
                                    message = "The filename: " + str(new) + " already exists. Please use a different name."
                                )
                                return
                        else:
                            if os.path.isfile(GlobalSettings.CSPR_DB + "\\FNA\\" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "\\FNA\\" + orig, GlobalSettings.CSPR_DB + "\\FNA\\" + new)
                            else:
                                show_message(
                                    fontSize = self.fontSize,
                                    icon = QtWidgets.QMessageBox.Icon.Critical,
                                    title = "Renaming Error",
                                    message = "The filename: " + str(new) + " already exists. Please use a different name."
                                )
                                return
                    else:
                        #unix cannot have spaces in paths
                        new = new.replace(" ","")
                        if new.find(".gbff") != -1:
                            if os.path.isfile(GlobalSettings.CSPR_DB + "/GBFF/" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "/GBFF/" + orig, GlobalSettings.CSPR_DB + "/GBFF/" + new)
                            else:
                                show_message(
                                    fontSize = self.fontSize,
                                    icon = QtWidgets.QMessageBox.Icon.Critical,
                                    title = "Renaming Error",
                                    message = "The filename: " + str(new) + " already exists. Please use a different name."
                                )
                                return
                        else:
                            if os.path.isfile(GlobalSettings.CSPR_DB + "/FNA/" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "/FNA/" + orig, GlobalSettings.CSPR_DB + "/FNA/" + new)
                            else:
                                show_message(
                                    fontSize = self.fontSize,
                                    icon = QtWidgets.QMessageBox.Icon.Critical,
                                    title = "Renaming Error",
                                    message = "The filename: " + str(new) + " already exists. Please use a different name."
                                )
                                return

            self.rename_window.rename_table.setRowCount(0)
            self.rename_window.close()
            GlobalSettings.mainWindow.fill_annotation_dropdown()
            center_ui(self.goToPrompt)
            self.goToPrompt.show()
            self.goToPrompt.activateWindow()
        except Exception as e:
            show_error("Error in submit_rename() in ncbi tool.", e)

    def stay(self):
        try:
            self.goToPrompt.hide()
        except Exception as e:
            show_error("Error in stay() in ncbi tool.", e)

    def close(self):
        try:
            self.hide()
            self.goToPrompt.hide()
        except Exception as e:
            show_error("Error in close() in ncbi tool.", e)

class goToPrompt(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(goToPrompt, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/ncbi_nav_page.ui', self)
            self.setWindowTitle("Navigate")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.label.setText("Successfully downloaded file(s) to the CASPERdb directory:\n" + GlobalSettings.CSPR_DB)
            self.label.setAlignment(QtCore.Qt.AlignCenter)
            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            font: bold 14pt 'Arial';
                            margin-top: 10px;}"""
            self.groupBox.setStyleSheet(groupbox_style)

            # self.scaleUI()
            scale_ui(self, custom_scale_width=575, custom_scale_height=225)
        except Exception as e:
            show_error("Error initializing goToPrompt class in ncbi tool.", e)

class rename_window(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(rename_window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "ui/ncbi_rename_window.ui", self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Rename Files")
            self.rename_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
            self.rename_table.setColumnCount(2)
            self.rename_table.setHorizontalHeaderLabels(['Original Filename', 'New Filename'])

            scale_ui(self, custom_scale_width=575, custom_scale_height=300)
        except Exception as e:
            show_error("Error initializing rename_window class in ncbi tool.", e)
   
#loading window UI class for when data is loading
class loading_window(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(loading_window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "ui/loading_data_form.ui", self)
            self.loading_bar.setValue(0)
            self.setWindowTitle("Loading NCBI Data")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            scale_ui(self, custom_scale_width=450, custom_scale_height=125)
        except Exception as e:
            show_error("Error initializing loading_window class in ncbi tool.", e)
