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
import traceback
import math

#global logger
logger = GlobalSettings.logger

ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "casper2informatics@gmail.com"

#model for filtering columns in ncbi table
class CustomProxyModel(QtCore.QSortFilterProxyModel):

    def __init__(self, parent=None):
        try:
            super().__init__(parent)
            self._filters = dict()
        except Exception as e:
            logger.critical("Error initializing CustomProxyModel class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    @property
    def filters(self):
        try:
            return self._filters
        except Exception as e:
            logger.critical("Error in filter() in custom proxy model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def setFilter(self, expresion, column):
        try:
            if expresion:
                self.filters[column] = expresion
            elif column in self.filters:
                del self.filters[column]
            self.invalidateFilter()
        except Exception as e:
            logger.critical("Error in setFilters() in custom proxy model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
            logger.critical("Error in filterAcceptsRow() in custom proxy model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)


#model for the data in the ncbi search table
class PandasModel(QtCore.QAbstractTableModel):

    def __init__(self, df=pd.DataFrame(), parent=None):
        try:
            QtCore.QAbstractTableModel.__init__(self, parent=parent)
            self._df = df.copy()
        except Exception as e:
            logger.critical("Error initializing PandasModel class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def toDataFrame(self):
        try:
            return self._df.copy()

        except Exception as e:
            logger.critical("Error in toDataFrame() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
            logger.critical("Error in headerData() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
            logger.critical("Error in data() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
            logger.critical("Error in setData() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def rowCount(self, parent=QtCore.QModelIndex()):
        try:
            return len(self._df.index)

        except Exception as e:
            logger.critical("Error in rowCount() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def columnCount(self, parent=QtCore.QModelIndex()):
        try:
            return len(self._df.columns)
        except Exception as e:
            logger.critical("Error in columnCount() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def sort(self, column, order):
        try:
            colname = self._df.columns.tolist()[column]
            self.layoutAboutToBeChanged.emit()
            self._df.sort_values(colname, ascending=order, inplace=True)
            self._df.reset_index(inplace=True, drop=True)
            self.layoutChanged.emit()
        except Exception as e:
            logger.critical("Error in sort() in Pandas Model in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)


#ncbi
class NCBI_search_tool(QtWidgets.QMainWindow):

    def __init__(self):
        try:
            super(NCBI_search_tool, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ncbi.ui', self)
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
            self.rename_window.go_back.clicked.connect(self.go_back)
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

            #scale UI
            self.first_show = True
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing NCBI_search_tool class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.fontSize = fontSize
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            # CASPER header scaling
            fontSize = 30
            self.title.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 850x750
            scaledWidth = int((width * 1000) / 1920)
            scaledHeight = int((height * 750) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(scaledWidth / 2))
            y = y - (math.ceil(scaledHeight / 2))
            self.setGeometry(x, y, scaledWidth, scaledHeight)

            self.repaint()
            QtWidgets.QApplication.processEvents()

        except Exception as e:
            logger.critical("Error in scaleUI() in NCBI tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            # center window on current screen
            width = self.width()
            height = self.height()
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in centerUI() in NCBI tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
            logger.critical("Error in go_back() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    @QtCore.pyqtSlot()
    def query_db(self):
        try:
            #show loading
            self.loading_window.loading_bar.setValue(5)
            self.loading_window.centerUI()
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
            #bs_content = BeautifulSoup(content, "lxml")
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
            #bs_content = BeautifulSoup(content, 'lxml')
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
                #temp_str = BeautifulSoup(temp_str, 'lxml')
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
            self.df = pd.DataFrame({'ID': ids,
                                'Species Name' : species_name,
                                'Strain' : strains,
                               'Assembly Name' : assembly_name,
                               'GenBank assembly accession': genbank_ids,
                               'RefSeq assembly accession': refseq_ids,
                               'Assembly Status': assembly_status})

            self.loading_window.loading_bar.setValue(90)
            QtCore.QCoreApplication.processEvents()

            #Build table view
            self.df.replace('', 'N/A', inplace=True)
            self.model = PandasModel(self.df)
            self.proxy = CustomProxyModel(self)
            self.proxy.setSourceModel(self.model)
            self.ncbi_table.setModel(self.proxy)
            self.ncbi_table.resizeColumnsToContents()
            self.comboBox.addItems(["{0}".format(col) for col in self.model._df.columns])
            self.activateWindow()

            #close loading gif
            self.loading_window.hide()
            self.loading_window.loading_bar.setValue(0)
            QtCore.QCoreApplication.processEvents()
        except Exception as e:
            logger.critical("Error in query_db() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
        except Exception as e:
            logger.critical("Error in on_view_horizontalHeader_sectionClicked() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    @QtCore.pyqtSlot()
    def on_actionAll_triggered(self):
        try:
            filterColumn = self.logicalIndex
            self.proxy.setFilter("", filterColumn)
        except Exception as e:
            logger.critical("Error in on_actionAll_triggered() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    @QtCore.pyqtSlot(int)
    def on_signalMapper_mapped(self, i):
        try:
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
        except Exception as e:
            logger.critical("Error in on_signalMapper_mapped() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    @QtCore.pyqtSlot()
    def download_files_wrapper(self):
        try:
            self.progressBar.setValue(0)

            #make sure rows are present in table
            if self.df.shape[0] == 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No Query Results")
                msgBox.setText("Please run an NCBI query to fill the table with results to choose from!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            #make sure user has selected at least one row
            indexes = self.ncbi_table.selectionModel().selectedRows()
            if len(indexes) == 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No Rows Selected")
                msgBox.setText("Please select rows from the table!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            #make sure file type is selected
            if self.gbff_checkbox.isChecked() == False and self.fna_checkbox.isChecked() == False:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No File Type Selected")
                msgBox.setText("No file type selected. Please select the file types you want to download!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            self.download_files()
        except Exception as e:
            logger.critical("Error in download_files_wrapper() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def download_files(self):
        try:
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


            if len(files) > 0:
                self.rename_files(files)
            else:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No Files Downloaded")
                msgBox.setText("No files were downloaded from the selected NCBI files. Please make sure the selected files are available in the database selected.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()
                return

        except Exception as e:
            logger.critical("Error in download_files() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    @QtCore.pyqtSlot()
    def select_all(self):
        try:
            if self.all_rows.isChecked():
                self.ncbi_table.selectAll()
            else:
                self.ncbi_table.clearSelection()
        except Exception as e:
            logger.critical("Error in select_all() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    # decompress file function
    def decompress_file(self, filename):
        try:
            with gzip.open(str(filename), 'rb') as f_in:
                with open(str(filename).replace('.gz', ''), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(str(filename))
        except Exception as e:
            logger.critical("Error in decompress_file() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
                    self.rename_window.rename_table.setItem(cnt, 0, QtWidgets.QTableWidgetItem(file))
                    self.rename_window.rename_table.setCellWidget(cnt, 1, QtWidgets.QLineEdit())
                    cnt += 1

                header = self.rename_window.rename_table.horizontalHeader()
                self.rename_window.rename_table.resizeColumnsToContents()
                header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
                header.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
                self.rename_window.centerUI()
                self.rename_window.show()
                self.rename_window.activateWindow()
            else:
                GlobalSettings.mainWindow.fill_annotation_dropdown()
                self.goToPrompt.show()
                self.goToPrompt.activateWindow()
        except Exception as e:
            logger.critical("Error in rename_files() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def rename_go_back(self):
        try:
            self.rename_window.close()
        except Exception as e:
            logger.critical("Error in rename_go_back() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

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
                                msgBox = QtWidgets.QMessageBox()
                                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                                msgBox.setWindowTitle("Renaming Error")
                                msgBox.setText("The filename: " + str(new) + " already exists. Please use a different name." )
                                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                                msgBox.exec()
                                return
                        else:
                            if os.path.isfile(GlobalSettings.CSPR_DB + "\\FNA\\" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "\\FNA\\" + orig, GlobalSettings.CSPR_DB + "\\FNA\\" + new)
                            else:
                                msgBox = QtWidgets.QMessageBox()
                                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                                msgBox.setWindowTitle("Renaming Error")
                                msgBox.setText("The filename: " + str(new) + " already exists. Please use a different name." )
                                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                                msgBox.exec()
                                return

                    else:
                        #unix cannot have spaces in paths
                        new = new.replace(" ","")
                        if new.find(".gbff") != -1:
                            if os.path.isfile(GlobalSettings.CSPR_DB + "/GBFF/" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "/GBFF/" + orig, GlobalSettings.CSPR_DB + "/GBFF/" + new)
                            else:
                                msgBox = QtWidgets.QMessageBox()
                                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                                msgBox.setWindowTitle("Renaming Error")
                                msgBox.setText(
                                    "The filename: " + str(new) + " already exists. Please use a different name.")
                                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                                msgBox.exec()
                                return
                        else:
                            if os.path.isfile(GlobalSettings.CSPR_DB + "/FNA/" + new) == False:
                                os.rename(GlobalSettings.CSPR_DB + "/FNA/" + orig, GlobalSettings.CSPR_DB + "/FNA/" + new)
                            else:
                                msgBox = QtWidgets.QMessageBox()
                                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                                msgBox.setWindowTitle("Renaming Error")
                                msgBox.setText(
                                    "The filename: " + str(new) + " already exists. Please use a different name.")
                                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                                msgBox.exec()
                                return

            self.rename_window.rename_table.setRowCount(0)
            self.rename_window.close()
            GlobalSettings.mainWindow.fill_annotation_dropdown()
            self.goToPrompt.centerUI()
            self.goToPrompt.show()
            self.goToPrompt.activateWindow()
        except Exception as e:
            logger.critical("Error in submit_rename() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def stay(self):
        try:
            self.goToPrompt.hide()
        except Exception as e:
            logger.critical("Error in stay() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    def close(self):
        try:
            self.hide()
            self.goToPrompt.hide()
        except Exception as e:
            logger.critical("Error in close() in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)


class goToPrompt(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(goToPrompt, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ncbi_nav_page.ui', self)
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

            self.scaleUI()
        except Exception as e:
            logger.critical("Error initializing goToPrompt class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 550x200
            scaledWidth = int((width * 575) / 1920)
            scaledHeight = int((height * 225) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(scaledWidth / 2))
            y = y - (math.ceil(scaledHeight / 2))
            self.setGeometry(x, y, scaledWidth, scaledHeight)

            self.repaint()
            QtWidgets.QApplication.processEvents()

        except Exception as e:
            logger.critical("Error in scaleUI() in goToPrompt in NCBI tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            # center window on current screen
            width = self.width()
            height = self.height()
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in centerUI() in goToPrompt in NCBI tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)


class rename_window(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(rename_window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "ncbi_rename_window.ui", self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Rename Files")
            self.rename_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
            self.rename_table.setColumnCount(2)
            self.rename_table.setHorizontalHeaderLabels(['Original Filename', 'New Filename'])

            self.scaleUI()
        except Exception as e:
            logger.critical("Error initializing rename_window class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 550x200
            scaledWidth = int((width * 575) / 1920)
            scaledHeight = int((height * 300) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(scaledWidth / 2))
            y = y - (math.ceil(scaledHeight / 2))
            self.setGeometry(x, y, scaledWidth, scaledHeight)

            self.repaint()
            QtWidgets.QApplication.processEvents()

        except Exception as e:
            logger.critical("Error in scaleUI() in rename window in NCBI tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            # center window on current screen
            width = self.width()
            height = self.height()
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in centerUI() in rename window in NCBI tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)


#loading window UI class for when data is loading
class loading_window(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(loading_window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "loading_data_form.ui", self)
            self.loading_bar.setValue(0)
            self.setWindowTitle("Loading NCBI Data")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.scaleUI()
        except Exception as e:
            logger.critical("Error initializing loading_window class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # scale/center window
            scaledWidth = int((width * 450) / 1920)
            scaledHeight = int((height * 125) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(scaledWidth / 2))
            y = y - (math.ceil(scaledHeight / 2))
            self.setGeometry(x, y, scaledWidth, scaledHeight)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in scaleUI() in loading_window() class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            width = self.width()
            height = self.height()
            #scale/center window
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in centerUI() in loading_window() class in ncbi tool.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)
