from Algorithms import get_table_headers
import sys
import os
import io
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from CoTargeting import CoTargeting
from closingWin import closingWindow
from Results import Results
from NewGenome import NewGenome
from NewEndonuclease import NewEndonuclease
import genomeBrowser
import webbrowser
import requests
import GlobalSettings
import multitargeting
from AnnotationParser import Annotation_Parser
from export_tool import export_tool
from generateLib import genLibrary
from CSPRparser import CSPRparser
import populationAnalysis
import platform
import ncbi
import glob
import traceback
import math
import logging
from annotation_functions import *

#logger alias for global logger
logger = GlobalSettings.logger

#Annotation file and search query from MainWindow
class AnnotationsWindow(QtWidgets.QMainWindow):
    #init annotation window class
    def __init__(self, info_path):
        super(AnnotationsWindow, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'annotation_details.ui', self)
        self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
        self.Submit_button.clicked.connect(self.submit)
        self.Go_Back_Button.clicked.connect(self.go_Back)
        self.select_all_checkbox.stateChanged.connect(self.select_all_genes)
        self.mainWindow = ""
        self.type = ""
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.switcher_table = [1, 1, 1, 1, 1, 1, 1, 1]
        self.tableWidget.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.tableWidget.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.tableWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.tableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tableWidget.setAutoScroll(False)
        self.tableWidget.horizontalHeader().sectionClicked.connect(self.table_sorting)

        #scale UI
        self.scaleUI()

    def table_sorting(self,logicalIndex):
        try:
            self.switcher_table[logicalIndex] *= -1
            if self.switcher_table[logicalIndex] == -1:
                self.tableWidget.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.tableWidget.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)
        except Exception as e:
            logger.critical("Error in table_sorting() in Annotation Window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
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
            # 16px is used for 92 dpi / 1920x1080
            fontSize = 12
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )
            self.menuBar().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )

            #CASPER header scaling
            fontSize = 30
            self.label_8.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            scaledWidth = int((width * 900) / 1920)
            scaledHeight = int((height * 600) / 1080)

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
            logger.critical("Error in scaleUI() in AnnotationWindow.")
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

            #center UI on current screen
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
        except Exception as e:
            logger.critical("Error in centerUI() in AnnotationWindow.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

            QtWidgets.QApplication.processEvents()
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #submit selected rows for results to process
    def submit(self):
        try:
            self.mainWindow.collect_table_data_nonkegg()
            self.mainWindow.show() # Open main window back up
            self.hide() # Close annotation window upon Submit
        except Exception as e:
            logger.critical("Error in submit() in AnnotationsWindow.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

            self.hide()
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #go back to main
    def go_Back(self):
        try:
            self.tableWidget.clear()
            self.mainWindow.searches.clear()
            self.tableWidget.setColumnCount(0)

            self.mainWindow.show()
            self.mainWindow.progressBar.setValue(0)
            self.hide()
        except Exception as e:
            logger.critical("Error in go_Back() in AnnotationsWindow.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

            self.mainWindow.checkBoxes.clear()
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # this function is very similar to the other fill_table, it just works with the other types of annotation files
    def fill_table_nonKegg(self, mainWindow,results_list):
        try:
            self.tableWidget.clearContents()
            self.mainWindow = mainWindow
            self.tableWidget.setColumnCount(5)
            self.mainWindow.progressBar.setValue(85)
            self.tableWidget.setHorizontalHeaderLabels(["Feature Type","Chromosome/Scaffold #","Feature ID/Locus Tag","Feature Name","Feature Description"])
            header = self.tableWidget.horizontalHeader()
            mainWindow.checkBoxes = []
            self.type = "nonkegg"
            index = 0
            for result in results_list:
                self.tableWidget.setRowCount(index + 1) # Increment table row count

                ### Set temp values
                chrom = result[0]
                feature = result[1]

                ### Create table items
                feature_type = QtWidgets.QTableWidgetItem(feature.type)
                chrom_number = QtWidgets.QTableWidgetItem(str(chrom))
                feature_id = QtWidgets.QTableWidgetItem(get_id(feature))
                feature_name = QtWidgets.QTableWidgetItem(get_name(feature))
                feature_desc = QtWidgets.QTableWidgetItem(get_description(feature))

                ### Set table items 
                self.tableWidget.setItem(index, 0, feature_type)
                self.tableWidget.setItem(index, 1, chrom_number)
                self.tableWidget.setItem(index, 2, feature_id)
                self.tableWidget.setItem(index, 3, feature_name)
                self.tableWidget.setItem(index, 4, feature_desc)

                mainWindow.checkBoxes.append((result[0],feature,index)) # Append chromosome number, SeqFeature object, and index in table to checkBoxes list. This list will be indexed later to return selected rows in the table by matching indices

                index +=1
                if index >= 1000:
                    break

            index = 0
            self.tableWidget.resizeColumnsToContents()
            mainWindow.hide()

            #center on current screen
            self.centerUI()
            self.show()
            self.activateWindow()

            return 0
        except Exception as e:
            logger.critical("Error in fill_table_nonKegg() in AnnotationsWindow.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

            self.mainWindow = mainWindow
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # this is the connection for the select all checkbox - selects/deselects all the genes in the table
    def select_all_genes(self):
        try:
            # check to see if we're selecting all of them or not
            if self.select_all_checkbox.isChecked():
                select_all = True
            else:
                select_all = False

            # # go through and do the selection
            # for i in range(self.tableWidget.rowCount()):
            #     #self.tableWidget.cellWidget(i, 4).setChecked(select_all)
            if select_all == True:
                self.tableWidget.selectAll()
            else:
                self.tableWidget.clearSelection()
        except Exception as e:
            logger.critical("Error in select_all_genes() in AnnotationsWindow.")
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

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        try:
            GlobalSettings.mainWindow.closeFunction()
        except Exception as e:
            logger.critical("Error in closeEvent() in AnnotationsWindow.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

            event.accept()
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)


# =========================================================================================
# CLASS NAME: CMainWindow
# Inputs: Takes in the path information from the startup window and also all input parameters
# Outputs: The results of the target search process by generating a new Results window
# =========================================================================================
class CMainWindow(QtWidgets.QMainWindow):

    def __init__(self, info_path):
        super(CMainWindow, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'CASPER_main.ui', self)
        self.dbpath = ""
        self.inputstring = "" # This is the search string
        self.info_path = info_path
        self.anno_name = ""
        self.endo_name = ""
        self.org = ""
        self.TNumbers = {}  # the T numbers from a kegg search
        self.orgcodes = {}  # Stores the Kegg organism code by the format {full name : organism code}
        self.gene_list = {}  # list of genes (no ides what they pertain to
        self.searches = {}
        self.checkBoxes = []
        self.genlib_list = [] # This list stores selected SeqFeatures from annotation window
        self.checked_info = {}
        self.check_ntseq_info = {}  # the ntsequences that go along with the checked_info
        self.annotation_parser = Annotation_Parser()
        self.link_list = list()  # the list of the downloadable links from the NCBI search
        self.organismDict = dict()  # the dictionary for the links to download. Key is the description of the organism, value is the ID that can be found in link_list
        self.results_list = list()
        self.organismData = list()
        self.ncbi = ncbi.NCBI_search_tool()

        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        margin-top: 10px;
                        font: bold 14pt 'Arial';}
                        """

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1", "Step3"))
        self.CASPER_Navigation.setStyleSheet(groupbox_style.replace("Step1", "CASPER_Navigation").replace("solid","dashed").replace("rgb(111,181,110)","rgb(88,89,91)"))

        self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
        self.pushButton_FindTargets.clicked.connect(self.gather_settings)
        self.pushButton_ViewTargets.clicked.connect(self.view_results)
        self.pushButton_ViewTargets.setEnabled(False)
        self.GenerateLibrary.setEnabled(False)
        self.radioButton_Gene.clicked.connect(self.toggle_annotation)
        self.radioButton_Position.clicked.connect(self.toggle_annotation)
    
        """ Connect functions to buttons """ 
        self.newGenome_button.clicked.connect(self.launch_newGenome) # Connect launch function to New Genome
        self.newEndo_button.clicked.connect(self.launch_newEndonuclease) # Connect launch function to New Endonuclease
        self.multitargeting_button.clicked.connect(self.changeto_multitargeting) # Connect launch function to Multitargeting
        self.populationAnalysis_button.clicked.connect(self.changeto_population_Analysis) # Connect launch function to PA
        self.GenerateLibrary.clicked.connect(self.prep_genlib)

        """ Connect functions to actions (toolbar) """ 
        self.actionOpen_Genome_Browser.triggered.connect(self.launch_newGenomeBrowser)
        self.actionExit.triggered.connect(self.close_app)
        self.visit_repo.triggered.connect(self.visit_repo_func)
        self.actionChange_Directory.triggered.connect(self.change_directory)
        self.actionNCBI.triggered.connect(self.open_ncbi_web_page)
        # self.actionCasper2.triggered.connect(self.open_casper2_web_page)
        self.actionNCBI_BLAST.triggered.connect(self.open_ncbi_blast_web_page)



        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()
        self.Annotation_Window = AnnotationsWindow(info_path)
        self.geneEntryField.setPlaceholderText("Example Inputs: \n\n"
                                         "Option 1: Feature (ID, Locus Tag, or Name)\n"
                                         "Example: 854068/YOL086C/ADH1 for S. cerevisiae alcohol dehydrogenase 1\n\n"
                                         "Option 2: Position (chromosome,start,stop)\n"
                                         "Example: 1,1,1000 for targeting chromosome 1, base pairs 1 to 1000\n\n"
                                         "Option 3: Sequence (must be within the selected organism)\n"
                                         "Example: Any nucleotide sequence between 100 and 10,000 base pairs.\n\n"
                                         "*Note: to multiplex, separate multiple queries by new lines*\n"
                                         "Example:\n"
                                         "1,1,1000\n"
                                         "5,1,500\n"
                                         "etc.")

        # show functionalities on window
        self.newGenome = NewGenome(info_path)
        self.newEndonuclease = NewEndonuclease()
        self.CoTargeting = CoTargeting(info_path)
        self.Results = Results()
        self.export_tool_window = export_tool()
        self.genLib = genLibrary()
        self.myClosingWindow = closingWindow()
        self.genomebrowser = genomeBrowser.genomebrowser()
        self.launch_ncbi_button.clicked.connect(self.launch_ncbi)

        #scale UI
        self.first_show = True
        self.scaleUI()

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
            # 16px is used for 92 dpi / 1920x1080
            fontSize = 12
            self.fontSize = fontSize
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )
            self.menuBar().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )

            #CASPER header scaling
            fontSize = 30
            self.label_8.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            #window resize and center
            scaledWidth = int((width * 1150) / 1920)
            scaledHeight = int((height * 650) / 1080)

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
            logger.critical("Error in scaleUI() in main.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

# that define the search for targets e.g. endonuclease, organism genome, gene target etc.
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

            width = self.width()
            height = self.height()
            # scale/center window
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
        except Exception as e:
            logger.critical("Error in centerUI() in main.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

            QtWidgets.QApplication.processEvents()
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # this function prepares everything for the generate library function
    # it is very similar to the gather settings, how ever it stores the data instead of calling the Annotation Window class
    # it moves the data onto the generateLib function, and then opens that window
    def prep_genlib(self):
        # make sure the user actually inputs something
        try:
            inputstring = str(self.geneEntryField.toPlainText())
            if (inputstring.startswith("Example Inputs:") or inputstring == ""):
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Error")
                msgBox.setText("No gene has been entered.  Please enter a gene.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

            else:
                # standardize the input
                inputstring = inputstring.lower()
                found_matches_bool = True
                # call the respective function
                self.progressBar.setValue(10)
                if self.radioButton_Gene.isChecked():
                    if len(self.genlib_list) > 0:
                        found_matches_bool = True
                    else:
                        found_matches_bool = False
                elif self.radioButton_Position.isChecked() or self.radioButton_Sequence.isChecked():
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("Error")
                    msgBox.setText("Generate Library can only work with feature searches.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    return
                """
                elif self.radioButton_Position.isChecked():
                    pinput = inputstring.split(';')
                    found_matches_bool = self.run_results("position", pinput,openAnnoWindow=False)
                elif self.radioButton_Sequence.isChecked():
                    sinput = inputstring
                    found_matches_bool = self.run_results("sequence", sinput, openAnnoWindow=False)
                """
                # if matches are found
                if found_matches_bool == True:
                    # get the cspr file name
                    cspr_file = self.organisms_to_files[self.orgChoice.currentText()][self.endoChoice.currentText()][0]
                    if platform.system() == 'Windows':
                        cspr_file = GlobalSettings.CSPR_DB + '\\' + cspr_file
                    else:
                        cspr_file = GlobalSettings.CSPR_DB + '/' + cspr_file
                    kegg_non = 'non_kegg'

                    # launch generateLib
                    self.progressBar.setValue(100)

                    # calculate the total number of matches found
                    tempSum = len(self.genlib_list)

                    # warn the user if the number is greater than 50
                    if tempSum > 50:
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                        msgBox.setWindowTitle("Many Matches Found")
                        msgBox.setText("More than 50 matches have been found. Continuing could cause a slow down...\n\n Do you wish to continue?")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
                        msgBox.exec()

                        if (msgBox.result() == QtWidgets.QMessageBox.No):
                            self.searches.clear()
                            self.progressBar.setValue(0)
                            return -2

                    self.genLib.launch(self.genlib_list,cspr_file, kegg_non)
                else:
                    self.progressBar.setValue(0)
        except Exception as e:
            logger.critical("Error in prep_genlib() in main.")
            logger.critical(e)
            logger.critical(traceback.format_exc())

        #print("prep genlib")
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # Function for collecting the settings from the input field and transferring them to run_results
    def gather_settings(self):
        try:
            ### If user searches multiple times for the same thing, this avoids re-searching the entire annotation file
            check_org = self.orgChoice.currentText().lower()
            check_endo = self.endoChoice.currentText().lower()
            check_anno_name = self.annotation_files.currentText().lower()
            check_input = str(self.geneEntryField.toPlainText()).lower()
            if (check_input == self.inputstring and check_org == self.org and check_anno_name == self.anno_name and check_endo == self.endo_name):
                same_search = True
            else:
                self.org = check_org
                self.anno_name = check_anno_name 
                self.inputstring = check_input
                self.endo_name = check_endo
                same_search = False

            # Error check: make sure the user actually inputs something
            if (self.inputstring.startswith("Example Inputs:") or self.inputstring == ""):
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Error")
                msgBox.setText(
                    "No feature has been searched for. Please enter a search.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

            else:

                ### Remove additional scoring columns if necessary
                header = get_table_headers(self.Results.targetTable) # Returns headers of the target table in View Targets window 
                col_indices = [header.index(x) for x in GlobalSettings.algorithms if x in header] # Returns the index(es) of the alternative scoring column(s) in the target table of View Targets window
                if len(col_indices) > 0: # If alternative scoring has been done
                    for i in col_indices:
                        self.Results.targetTable.removeColumn(i)
                self.Results.targetTable.resizeColumnsToContents()

                self.progressBar.setValue(10)
                if self.radioButton_Gene.isChecked():
                    ginput = [x.strip() for x in self.inputstring.split('\n')] # Split search based on newline character and remove deadspace
                    self.run_results("feature", ginput, same_search)
                elif self.radioButton_Position.isChecked():
                    pinput = [x.strip() for x in self.inputstring.split('\n')] # Split search based on newline character and remove deadspace
                    self.run_results("position", pinput, same_search)
                elif self.radioButton_Sequence.isChecked():
                    sinput = self.inputstring
                    self.run_results("sequence", sinput, same_search)
        except Exception as e:
            logger.critical("Error in gather_settings() in main.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # ---- Following functions are for running the auxillary algorithms and windows ---- #
    # this function is parses the annotation file given, and then goes through and goes onto results
    # it will call other versions of collect_table_data and fill_table that work with these file types
    # this function should work with the any type of annotation file, besides kegg.
    # this assumes that the parsers all store the data the same way, which gff and feature table do
    # please make sure the genbank parser stores the data in the same way
    # so far the gff files seems to all be different. Need to think about how we want to parse it
    def run_results_own_ncbi_file(self, inputstring, fileName, same_search, openAnnoWindow=True):
        try:
            self.progressBar.setValue(35)
            ### Now actually search for inputs in annotation file
            self.results_list = self.annotation_parser.genbank_search(inputstring,same_search)

            ### Quick error check to make sure the chromosome numbers match
            cspr_file = self.organisms_to_files[self.orgChoice.currentText()][self.endoChoice.currentText()][0]
            if platform.system() == 'Windows':
                cspr_file = GlobalSettings.CSPR_DB + '\\' + cspr_file
            else:
                cspr_file = GlobalSettings.CSPR_DB + '/' + cspr_file

            own_cspr_parser = CSPRparser(cspr_file)
            own_cspr_parser.read_first_lines()

            if len(own_cspr_parser.karystatsList) != self.annotation_parser.max_chrom:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Warning:")
                msgBox.setText(
                    "The number of chromosomes do not match. This could cause errors.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

            self.progressBar.setValue(60)

            # now go through and search for the actual locus tag, in the case the user input that
            searchValues = self.separate_line(inputstring[0])
            self.searches.clear()

            self.progressBar.setValue(75)
            if len(self.results_list) <= 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No Matches Found")
                msgBox.setText(
                    "No matches found with that search, please try again.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                self.progressBar.setValue(0)
                if openAnnoWindow:
                    return
                else:
                    return False

            # if we get to this point, that means that the search yieleded results, so fill the table
            self.progressBar.setValue(80)
            # check whether this function call is for Annotation Window, or for generate Lib
            if openAnnoWindow:
                self.Annotation_Window.fill_table_nonKegg(self,self.results_list)
            else:
                return True
        except Exception as e:
            logger.critical("Error in run_results_own_ncbi_file() in main.")
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

    def run_results(self, inputtype, inputstring, same_search, openAnnoWindow=True):
        try:
            fileName = self.annotation_files.currentText()
            # self.annotation_parser = Annotation_Parser()
            #get complete path of file
            for file in glob.glob(GlobalSettings.CSPR_DB + "/**/*.gb*", recursive=True):
                if file.find(fileName) != -1:
                    self.annotation_parser.annotationFileName = file
                    break
            self.Results.annotation_path = self.annotation_parser.annotationFileName  ### Set annotation path
            #print("run results")

            progvalue = 15
            self.searches = {}
            self.gene_list = {}
            self.progressBar.setValue(progvalue)

            try:
                self.Results.endonucleaseBox.currentIndexChanged.disconnect()
            except Exception as e:
                pass
            # set Results endo combo box
            self.Results.endonucleaseBox.clear()

            # set the results window endoChoice box menu
            # set the mainWindow's endoChoice first, and then loop through and set the rest of them
            self.Results.endonucleaseBox.addItem(self.endoChoice.currentText())
            for item in self.organisms_to_endos[str(self.orgChoice.currentText())]:
                if item != self.Results.endonucleaseBox.currentText():
                    self.Results.endonucleaseBox.addItem(item)

            self.Results.endonucleaseBox.currentIndexChanged.connect(self.Results.changeEndonuclease)
            self.Results.get_endo_data()

    #        self.Results.change_start_end_button.setEnabled(False)
            self.Results.displayGeneViewer.setChecked(0)

            if inputtype == "feature":

                fileType = self.annotation_parser.find_which_file_version()

                # if the parser retuns the 'wrong file type' error
                if fileType == -1:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("Error:")
                    msgBox.setText("Feature search requires a GenBank formatted annotation file. Please select a file from the dropdown menu or search by position")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    self.progressBar.setValue(0)
                    return

                # make sure an annotation file has been selected
                if self.annotation_files.currentText() == "None":
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("No Annotation")
                    msgBox.setText("Search by feature requires a GenBank annotation file. Please select one from the dropdown menu or search by position.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    self.progressBar.setValue(0)
                    return
                # this now just goes onto the other version of run_results
                myBool = self.run_results_own_ncbi_file(inputstring, self.annotation_files.currentText(), same_search, openAnnoWindow=openAnnoWindow)
                if not openAnnoWindow:
                    return myBool
                else:
                    self.progressBar.setValue(0)
                    return

            # position code below
            if inputtype == "position":
                full_org = str(self.orgChoice.currentText())
                self.checked_info.clear()
                self.check_ntseq_info.clear()

                for item in inputstring: # Loop through each search
                    searchIndices = [x.strip() for x in item.split(',')] # Parse input query
                    ### Make sure the right amount of arguments were passed
                    if len(searchIndices) != 3:
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Position Error: Invalid Input")
                        msgBox.setText(
                            "There are 3 arguments required for this function: chromosome, start position, and end position.")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        self.progressBar.setValue(0)
                        return

                    ### Make sure user inputs digits
                    if not searchIndices[0].isdigit() or not searchIndices[1].isdigit() or not searchIndices[2].isdigit():
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Position Error: Invalid Input")
                        msgBox.setText(
                            "The positions given must be integers. Please try again.")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        self.progressBar.setValue(0)
                        return

                    ### Make sure start is less than end
                    elif int(searchIndices[1]) >= int(searchIndices[2]):
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Position Error: Start Must Be Less Than End")
                        msgBox.setText(
                            "The start index must be less than the end index.")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        self.progressBar.setValue(0)
                        return

                    ### Make sure range isn't too large
                    elif abs(int(searchIndices[2])-int(searchIndices[1])) > 50000:
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Position Error: Range Too Large")
                        msgBox.setText(
                            "The search range must be less than 50,000 nt.")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        self.progressBar.setValue(0)
                        return

                    ### Make sure chromosome exists
                    elif int(searchIndices[0]) > self.annotation_parser.get_max_chrom():
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Position Error: Chromsome Doesn't Exist")
                        msgBox.setText(
                            "Chromosome %s does not exist in the selected annotation file." % searchIndices[0])
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        self.progressBar.setValue(0)
                        return

                    # append the data into the checked_info
                    tempString = 'chrom: ' + str(searchIndices[0]) + ',start: ' + str(searchIndices[1]) + ',end: ' + str(searchIndices[2])
                    self.checked_info[tempString] = (int(searchIndices[0]), int(searchIndices[1])-1, int(searchIndices[2]))

                self.progressBar.setValue(50)
                self.Results.transfer_data(full_org, self.organisms_to_files[full_org], [str(self.endoChoice.currentText())], os.getcwd(), self.checked_info, self.check_ntseq_info,inputtype)
                self.Results.load_gene_viewer()
                self.progressBar.setValue(100)
                self.pushButton_ViewTargets.setEnabled(True)
                self.GenerateLibrary.setEnabled(True)

            # sequence code below
            if inputtype == "sequence":
                fileType = self.annotation_parser.find_which_file_version()
                # if the parser retuns the 'wrong file type' error
                if fileType == -1:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("Error:")
                    msgBox.setText("Search by sequence requires a GenBank annotation file. Please select one from the dropdown menu or search by position.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    self.progressBar.setValue(0)
                    return
                if self.annotation_files.currentText() == "None":
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("Error:")
                    msgBox.setText("Search by sequence requires a GenBank annotation file. Please select one from the dropdown menu or search by position.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    self.progressBar.setValue(0)
                    return


                checkString = 'AGTCN'
                full_org = str(self.orgChoice.currentText())
                self.checked_info.clear()
                self.progressBar.setValue(10)
                inputstring = inputstring.replace('\n','').upper().strip()

                # make sure all the chars are one of A, G, T, C, or N
                for letter in inputstring:
                    if letter not in checkString:
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Sequence Error")
                        msgBox.setText(
                            "The sequence must consist of A, G, T, C, or N. No other characters are allowed.")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        self.progressBar.setValue(0)
                        return

                # check to make sure that the use gave a long enough sequence
                if len(inputstring) < 100:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("Error")
                    msgBox.setText(
                        "The sequence given is too small. At least 100 characters are required.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    self.progressBar.setValue(0)
                    return

                # give a warning if the length of the sequence is long
                if len(inputstring) > 10000:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                    msgBox.setWindowTitle("Large Sequence Detected")
                    msgBox.setText(
                        "The sequence given is too large one.\n\nPlease input a sequence less than 10kb in length.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                    msgBox.exec()

                    self.progressBar.setValue(0)
                    return

                self.progressBar.setValue(30)

                # Check the GBFF file for the sequence
                my_check = self.annotation_parser.get_sequence_info(inputstring)

                self.progressBar.setValue(55) # Update progress bar

                if type(my_check) == bool: # This means the sequence was not found
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                    msgBox.setWindowTitle("Sequence Not Found")
                    msgBox.setText(
                        "The sequence entered was not found.\n\nPlease input a sequence that is in the selected organism.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                    msgBox.exec()
                    self.progressBar.setValue(0)
                    return

                else: # This means the sequence was found
                    # append the data into the checked_info
                    tempString = 'chrom: ' + str(my_check[0]) + ',start: ' + str(my_check[1]) + ',end: ' + str(my_check[2])
                    self.checked_info[tempString] = (int(my_check[0]), int(my_check[1])-1, int(my_check[2]))

                self.progressBar.setValue(75) # Update progress bar

                self.Results.transfer_data(full_org, self.organisms_to_files[full_org], [str(self.endoChoice.currentText())], os.getcwd(), self.checked_info, self.check_ntseq_info, inputtype)
                self.Results.load_gene_viewer()
                self.progressBar.setValue(100)
                self.pushButton_ViewTargets.setEnabled(True)
                self.GenerateLibrary.setEnabled(True)


        except Exception as e:
            logger.critical("Error in run_results() in main.")
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

    #launch new genome tool
    def launch_newGenome(self):
        try:
            # update endo list
            self.newGenome.fillEndo()
            if self.newGenome.first_show == True:
                self.newGenome.centerUI()
                self.newGenome.first_show = False
            self.hide()
            self.newGenome.show()
        except Exception as e:
            logger.critical("Error in launch_newGenome() in main.")
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

    #launch new endo tool
    def launch_newEndonuclease(self):
        try:
            self.newEndonuclease.centerUI()
            self.newEndonuclease.show()
            self.newEndonuclease.activateWindow()
        except Exception as e:
            logger.critical("Error in launch_newEndonuclease() in main.")
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

    #launch genome browser tool
    def launch_newGenomeBrowser(self):
        try:
            self.genomebrowser.createGraph(self)
        except Exception as e:
            logger.critical("Error in launch_newGenomeBrowser() in main.")
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

    #launch ncbi tool
    def launch_ncbi(self):
        try:
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Information)
            msgBox.setWindowTitle("Note:")
            msgBox.setText(
                "NCBI Annotation Guidelines:\n\nDownload annotation files of the exact species and strain used in Analyze New Genome.\n\nMismatched annotation files will inhibit downstream analyses.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
            msgBox.exec()

            if self.ncbi.first_show == True:
                self.ncbi.first_show = False
                self.ncbi.centerUI()

            self.ncbi.show()
            self.ncbi.activateWindow()
        except Exception as e:
            logger.critical("Error in launch_ncbi() in main.")
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

    # this function does the same stuff that the other collect_table_data does, but works with the other types of files
    def collect_table_data_nonkegg(self):
        try:
            # start out the same as the other collect_table_data
            self.checked_info.clear()
            self.genlib_list.clear()
            self.check_ntseq_info.clear()
            full_org = str(self.orgChoice.currentText())
            holder = ()
            selected_indices = []
            selected_rows = self.Annotation_Window.tableWidget.selectionModel().selectedRows()
            for ind in sorted(selected_rows):
                selected_indices.append(ind.row())

            for item in self.checkBoxes:
                feature = item[1]
                # If inidices of checkBoxes list and selected rows in table match...
                if item[2] in selected_indices:
                    holder = (item[0],int(feature.location.start),int(feature.location.end)) # Tuple order: Feature chromosome/scaffold number, feature start, feature end
                    ### If locus tag available, combine with gene name to create dict key
                    if 'locus_tag' in feature.qualifiers:
                        tag = feature.qualifiers['locus_tag'][0]
                        key = tag + ": " + get_name(feature)
                    else:
                        key = get_name(feature)
                    self.checked_info[key] = holder
                    self.genlib_list.append((item[0],feature)) # Tuple order: Feature chromosome/scaffold number, SeqFeature object
                else:
                    # If item was not selected in the table, go to the next item
                    continue

            # now call transfer data
            self.progressBar.setValue(95)
            self.Results.transfer_data(full_org, self.organisms_to_files[full_org], [str(self.endoChoice.currentText())], os.getcwd(),
                                       self.checked_info, self.check_ntseq_info,inputtype="feature")
            self.Results.load_gene_viewer()

            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)
            self.GenerateLibrary.setEnabled(True)
        except Exception as e:
            logger.critical("Error in collect_table_data_nonkegg() in main.")
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

    def separate_line(self, input_string):
        try:
            export_array = []
            while True:
                index = input_string.find('\n')
                if index == -1:
                    if len(input_string) == 0:
                        return export_array
                    else:
                        export_array.append(input_string)
                        return export_array
                export_array.append(input_string[:index])
                input_string = input_string[index + 1:]
        except Exception as e:
            logger.critical("Error in seperate_line() in main.")
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

    def removeWhiteSpace(self, strng):
        try:
            while True:
                if len(strng) == 0 or (strng[0] != " " and strng[0] != "\n"):
                    break
                strng = strng[1:]
            while True:
                if len(strng) == 0 or (strng[len(strng) - 1] != " " and strng[0] != "\n"):
                    return strng
                strng = strng[:len(strng) - 1]
        except Exception as e:
            logger.critical("Error in removeWhiteSpace() in main.")
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

    # Function to enable and disable the Annotation function if searching by position or sequence
    def toggle_annotation(self):
        try:
            if self.radioButton_Gene.isChecked():
                self.Step2.setEnabled(True)
            else:
                self.Step2.setEnabled(True)
        except Exception as e:
            logger.critical("Error in toggle_annotation() in main.")
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

    def fill_annotation_dropdown(self):
        try:
            #recursive search for all GenBank files in casper db folder
            self.annotation_files.clear()
            annotation_files = glob.glob(GlobalSettings.CSPR_DB + "/**/*.gb*", recursive=True)
            if platform.system() == "Windows":
                for i in range(len(annotation_files)):
                    annotation_files[i] = annotation_files[i].replace("/","\\")
                    annotation_files[i] = annotation_files[i][annotation_files[i].rfind("\\") + 1:]
            else:
                for i in range(len(annotation_files)):
                    annotation_files[i] = annotation_files[i].replace("\\","/")
                    annotation_files[i] = annotation_files[i][annotation_files[i].rfind("/") + 1:]

            annotation_files.sort(key=str.lower)
            self.annotation_files.addItems(annotation_files)
            self.annotation_files.addItems(["None"])
        except Exception as e:
            logger.critical("Error in fill_annotation_dropdown() in main.")
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

    def make_dictonary(self):
        try:
            url = "https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:" + self.TNumbers[
                self.Annotations_Organism.currentText()]
            source_code = requests.get(url, verify=False)
            plain_text = source_code.text
            buf = io.StringIO(plain_text)

            while True:
                line = buf.readline()
                if line[0] == "-":
                    break
            while True:
                line = buf.readline()
                if line[1] != "a":
                    return
                line = line[line.find(">") + 1:]
                seq = line[line.find(":") + 1:line.find("<")]
                line = line[line.find(">") + 1:]

                i = 0
                while True:
                    if line[i] == " ":
                        i = i + 1
                    else:
                        break
                key = line[i:line.find("\n") - 1]
                if key in self.gene_list:
                    if seq not in self.gene_list[key]:
                        self.gene_list[key].append(seq)
                else:
                    self.gene_list[key] = [seq]
                z = 5
        except Exception as e:
            logger.critical("Error in make_dictionary() in main.")
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

    def organism_finder(self, long_str):
        try:
            semi = long_str.find(";")
            index = 1
            while True:
                if long_str[semi - index] == " ":
                    break
                index = index + 1
            return long_str[:semi - index]
        except Exception as e:
            logger.critical("Error trying in organism_finder() in main.")
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

    # This method is for testing the execution of a button call to make sure the button is linked properly
    def testexe(self):
        try:
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
            msgBox.setWindowTitle("Extract!")
            msgBox.setText(
                "Are you sure you want to quit?")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
            msgBox.exec()

            if msgBox.result() == QtWidgets.QMessageBox.Yes:
                # print(self.orgChoice.currentText())
                sys.exit()
            else:
                pass
        except Exception as e:
            logger.critical("Error in testexe() in main.")
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

    def getData(self):
        try:
            try:
                self.orgChoice.currentIndexChanged.disconnect()
            except Exception as e:
                 pass

            self.orgChoice.clear()
            self.endoChoice.clear()
            mypath = os.getcwd()
            found = False
            self.dbpath = mypath
            onlyfiles = [str(f) for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
            onlyfiles.sort(key=str.lower)
            self.organisms_to_files = {}
            self.organisms_to_endos = {}
            first = True
            for file in onlyfiles:
                if file.find('.cspr') != -1:
                    if first == True:
                        first = False
                    found = True
                    newname = file[0:-4]
                    endo = newname[newname.rfind("_")+1:-1]
                    hold = open(file, 'r')
                    buf = (hold.readline())
                    buf = str(buf)
                    buf = buf.strip()
                    species = buf.replace("GENOME: ",'')

                    if species in self.organisms_to_files:
                        self.organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]
                    else:
                        self.organisms_to_files[species] = {}
                        self.organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]

                    if species in self.organisms_to_endos:
                        self.organisms_to_endos[species].append(endo)
                    else:
                        self.organisms_to_endos[species] = [endo]
                        if self.orgChoice.findText(species) == -1:
                            self.orgChoice.addItem(species)

            #self.orgChoice.addItem("Custom Input Sequences")
            # auto fill the kegg search bar with the first choice in orgChoice
            if found == False:
                return False

            self.endoChoice.clear()
            self.endoChoice.addItems(self.organisms_to_endos[str(self.orgChoice.currentText())])
            self.orgChoice.currentIndexChanged.connect(self.changeEndos)
        except Exception as e:
            logger.critical("Error in getData() in main.")
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

    def changeEndos(self):
        try:
            if self.orgChoice.currentText() != "Custom Input Sequences":
                self.Step2.setEnabled(True)
                self.endoChoice.setEnabled(True)
                self.radioButton_Gene.show()
                self.radioButton_Position.show()
                self.endoChoice.clear()
                self.endoChoice.addItems(self.organisms_to_endos[str(self.orgChoice.currentText())])
            else:
                self.Step2.setEnabled(False)
                self.endoChoice.clear()
                self.endoChoice.setEnabled(False)
                self.radioButton_Gene.hide()
                self.radioButton_Position.hide()
        except Exception as e:
            logger.critical("Error in changeEndos() in main.")
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

    def change_directory(self):
        try:
            filed = QtWidgets.QFileDialog()
            mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a folder...",
                                                               self.dbpath, QtWidgets.QFileDialog.ShowDirsOnly)

            if os.path.isdir(mydir) == False:
                #check if directory is a valid directory
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Not a directory")
                msgBox.setText(
                    "The directory you selected does not exist.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return


            #check if directory contains CSPR files
            found = False
            for file in os.listdir(mydir):
                if (file.find(".cspr") != -1):
                    found = True
                    break
            if (found == False):
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Directory is invalid!")
                msgBox.setText(
                    "You must select a directory with CSPR Files!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return


            os.chdir(mydir)
            if platform.system() == "Windows":
                mydir = mydir.replace("/","\\")
            GlobalSettings.CSPR_DB = mydir
            #update dropdowns in main, MT, pop
            self.getData()

            GlobalSettings.MTWin.directory = mydir
            GlobalSettings.MTWin.get_data()
            GlobalSettings.pop_Analysis.get_data()
            self.fill_annotation_dropdown()
        except Exception as e:
            logger.critical("Error in change_directory() in main.")
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

    #change to multi-targeting window
    def changeto_multitargeting(self):
        try:
            os.chdir(os.getcwd())
            if GlobalSettings.MTWin.first_show == True:
                GlobalSettings.MTWin.show()
                GlobalSettings.MTWin.first_show = False
            else:
                GlobalSettings.MTWin.show()
                GlobalSettings.mainWindow.hide()

        except Exception as e:
            logger.critical("Error in changeto_multitargeting() in main.")
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

    #change to population analysis window
    def changeto_population_Analysis(self):
        try:
            GlobalSettings.pop_Analysis.launch()
            if GlobalSettings.pop_Analysis.first_show == True:
                GlobalSettings.pop_Analysis.centerUI()
                GlobalSettings.pop_Analysis.first_show = False
            GlobalSettings.pop_Analysis.show()
            GlobalSettings.mainWindow.hide()
        except Exception as e:
            logger.critical("Error in changeto_population_Analysis() in main.")
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

    def annotation_information(self):
        try:
            info = "Annotation files are used for searching for spacers on a gene/locus basis and can be selected here using either " \
                   "NCBI databases or a local file."
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Annotation Information")
            msgBox.setText(
                info)
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
            msgBox.exec()

        except Exception as e:
            logger.critical("Error in annotation_information() in main.")
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

    def open_ncbi_blast_web_page(self):
        try:
            webbrowser.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi', new=2)
        except Exception as e:
            logger.critical("Error in open_ncbi_blast_web_page() in main.")
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

    def open_ncbi_web_page(self):
        try:
            webbrowser.open('https://www.ncbi.nlm.nih.gov/', new=2)
        except Exception as e:
            logger.critical("Error in open_ncbi_web_page() in main.")
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

    # def open_casper2_web_page(self):
    #     try:
    #         webbrowser.open('http://casper2.org/', new=2)
    #     except Exception as e:
    #         logger.critical("Error in open_casper2_web_page() in main.")
    #         logger.critical(e)
    #         logger.critical(traceback.format_exc())
    #         exit(-1)

    def visit_repo_func(self):
        try:
            webbrowser.open('https://github.com/TrinhLab/CASPERapp')
        except Exception as e:
            logger.critical("Error in visit_repo_func() in main.")
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
    def view_results(self):
        try:
            #center results window on current screen
            if self.Results.first_show == True:
                self.Results.first_show = False
                self.Results.centerUI()

            self.Results.show()
            self.hide()
        except Exception as e:
            logger.critical("Error in view_results() in main.")
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

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        try:
            self.closeFunction()
            event.accept()
        except Exception as e:
            logger.critical("Error in closeEvent() in main.")
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

    def closeFunction(self):
        try:
            try:
                self.ncbi.close()
            except Exception as e:
                print("no ncbi window to close")
            self.myClosingWindow.get_files()
            self.myClosingWindow.centerUI()
            self.myClosingWindow.show()
        except Exception as e:
            logger.critical("Error in closeFunction() in main.")
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

    def close_app(self):
        try:
            try:
                self.ncbi.close()
            except Exception as e:
                print("no ncbi window to close")

            self.closeFunction()
            self.close()
        except Exception as e:
            logger.critical("Error in close_app() in main.")
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


#startup window class
class StartupWindow(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(StartupWindow, self).__init__()

            #load UX files
            try:
                uic.loadUi(GlobalSettings.appdir + 'startupCASPER.ui', self)
                self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
            except Exception as e:
                logger.critical("Unable to load UX files for Startup Window.")
                logger.critical(e)
                logger.critical(traceback.format_exc())
                exit(-1)

            #set "Main" button to be the default highlighted button on startup
            self.goToMain.setDefault(True)

            #get current directory, and update based on current operating system
            self.currentDirectory = os.getcwd()
            self.databaseDirectory = self.loadDatabaseDirectory()
            GlobalSettings.CSPR_DB = self.databaseDirectory
            if platform.system() == "Windows":
                GlobalSettings.CSPR_DB = GlobalSettings.CSPR_DB.replace("/","\\")
            else:
                GlobalSettings.CSPR_DB = GlobalSettings.CSPR_DB.replace("\\","/")

            #setup event handlers for startup buttons
            self.currentDirText.setText(self.databaseDirectory)
            self.changeDir.clicked.connect(self.changeDirectory)
            self.goToMain.clicked.connect(self.launchMainWindow)
            self.goToNewGenome.clicked.connect(self.launchNewGenome)

            self.setWindowTitle("CASPER")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing StartupWindow class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #function for scaling the font size and logo size based on resolution and DPI of screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = math.ceil(screen.physicalDotsPerInch())
            width = screen.geometry().width()
            height = screen.geometry().height()

            #log width x height
            logger.info("Resolution: %s x %s", width, height)

            # log DPI information
            logger.info("DPI = %d" % (dpi))

            # font scaling
            # 16px is used for 92 dpi 1920x1080
            fontSize = 12
            self.fontSize = fontSize
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )

            #set logo image
            pixmapOriginal = QtGui.QPixmap(GlobalSettings.appdir + "CASPER-logo.jpg")
            self.logo.setPixmap(pixmapOriginal)

            #scale and center UI
            scaledWidth = int((width * 850) / 1920)
            scaledHeight = int((height * 550) / 1080)
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
            logger.critical("Error in scaleUI() in startup window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            #center UI
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
            logger.critical("Error in centerUI() in startup window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #event handler for user clicking the "Change..." button - used for changing CASPER database directory
    def changeDirectory(self):
        try:
            #launch OS file browser
            fileBrowser = QtWidgets.QFileDialog()
            newDirectory = QtWidgets.QFileDialog.getExistingDirectory(fileBrowser, "Open a folder...",
                                                               self.databaseDirectory, QtWidgets.QFileDialog.ShowDirsOnly)
            #check if selected path is a directory in the system
            if (os.path.isdir(newDirectory) == False):
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Not a directory")
                msgBox.setText(
                    "The directory you selected does not exist.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            #make sure directory contains correct filepath format based on OS
            if platform.system() == "Windows":
                newDirectory = newDirectory.replace("/","\\")

            #update text edit showing the current selected database directory
            self.currentDirText.setText(newDirectory)

            #update casper database directories
            self.databaseDirectory = newDirectory
            GlobalSettings.CSPR_DB = newDirectory
        except Exception as e:
            logger.critical("Error in changeDirectory() in startup window.")
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

    #function for loading the default database directory specified in CASPERinfo
    #returns: default database parsed from CASPERinfo
    def loadDatabaseDirectory(self):
        try:
            #variable to hold default directory from CASPERinfo
            defaultDirectory = ""

            try:
                #open CASPERinfo file in application directory
                CASPERInfo = open(GlobalSettings.appdir + "CASPERinfo", 'r+')
                #read file, parse for the default database directory
                CASPERInfo = CASPERInfo.read()
                lines = CASPERInfo.split('\n')
                for item in lines:
                    if 'DIRECTORY:' in item:
                        #default directory found
                        defaultDirectory = item
                        break

                #remove lines meta-data
                defaultDirectory = defaultDirectory.replace("DIRECTORY:", "")

                #make sure directory is formatted properly based on OS
                if platform.system() == "Windows":
                    defaultDirectory = defaultDirectory.replace("/","\\")
                else:
                    defaultDirectory = defaultDirectory.replace("\\", "/")

                logger.debug("Successfully parsed CASPERinfo for default database directory.")
            except Exception as e:
                logger.error("Unable to read CASPERinfo file to get default database directory.")
                logger.error(e)
                logger.error(traceback.format_exc())
                return "Where would you like to store CASPER database files?"

            return defaultDirectory

        except Exception as e:
            logger.critical("Error in loadDatabaseDirectory() in startup window.")
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

    #function for saving the currently selected database directory to CASPERinfo to be the new default value on startup
    def saveDatabaseDirectory(self):
        try:
            #variable to hold the CASPERinfo data with new default directory change
            CASPERInfoNewData = ""

            #new default directory string for CASPERinfo
            newDefaultDirectory = "DIRECTORY:" + str(self.databaseDirectory)

            #open CASPERinfo file to read in the files data and add in new change
            try:
                CASPERInfo = open(GlobalSettings.appdir + "CASPERinfo", 'r+')
                CASPERinfoData = CASPERInfo.read()
                CASPERinfoData = CASPERinfoData.split('\n')
                for line in CASPERinfoData:
                    #if directory line found, use new default directory string instead
                    if 'DIRECTORY:' in line:
                        CASPERInfoNewData = CASPERInfoNewData + "\n" + newDefaultDirectory
                    else:
                        CASPERInfoNewData = CASPERInfoNewData + "\n" + line
                CASPERInfoNewData = CASPERInfoNewData[1:]

                #close CASPERinfo
                CASPERInfo.close()

                #re-open the file and re-write it with current changes
                CASPERInfo = open(GlobalSettings.appdir + "CASPERinfo", 'w+')
                CASPERInfo.write(CASPERInfoNewData)
                CASPERInfo.close()
                logger.debug("Successfully updated CASPERinfo with new default database directory.")
            except Exception as e:
                logger.critical("Unable to write to CASPERinfo file to update database directory.")
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
        except Exception as e:
            logger.critical("Error in saveDatabaseDirectory() in startup window.")
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

    # even handler for user clicking the "New Genome" button - used for launching New Genome
    def launchNewGenome(self):
        try:
            # make sure database directory variable is up-to-date based on what the user has in the text edit
            self.databaseDirectory = str(self.currentDirText.text())

            # make sure the path is a valid path before launching New Genome
            if (os.path.isdir(self.databaseDirectory)):

                # change directories to the specified database directory provided
                os.chdir(self.databaseDirectory)

                # write out the database directory to CASPERinfo to be the new default loaded value
                self.saveDatabaseDirectory()

                # update global database variable
                GlobalSettings.CSPR_DB = self.databaseDirectory

                # make sure FNA and GBFF subdirectories are present, if not create them
                subdirs = os.listdir(self.databaseDirectory)
                if "FNA" not in subdirs and os.path.isdir("FNA") == False:
                    try:
                        os.mkdir("FNA")
                        logger.debug("Successfully created FNA subdirectory in database directory.")
                    except Exception as e:
                        logger.critical("Unable to make 'FNA' subdirectory in database directory")
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
                if "GBFF" not in subdirs and os.path.isdir("GBFF") == False:
                    try:
                        os.mkdir("GBFF")
                        logger.debug("Successfully created GBFF subdirectory in database directory.")
                    except Exception as e:
                        logger.critical("Unable to make 'GBFF' subdirectory in database directory")
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

                # launch new genome
                try:
                    GlobalSettings.mainWindow.launch_newGenome()
                    logger.debug("Successfully initialized New Genome in startup window.")
                except Exception as e:
                    logger.critical("Unable to initialize New Genome from startup window.")
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

                self.close()
            else:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Not a directory")
                msgBox.setText(
                    "The directory you selected does not exist.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

        except Exception as e:
            logger.critical("Error in launchNewGenome() in startup window.")
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

    #event handler for user clicking "Main Program" button - used to launch Main Window
    def launchMainWindow(self):
        try:
            #make sure database directory variable is up-to-date based on what the user has in the text edit
            self.databaseDirectory = str(self.currentDirText.text())

            # make sure the path is a valid path before launching New Genome
            if os.path.isdir(self.databaseDirectory) == True:

                #check if database directory has CSPR files in it
                foundCSPRFiles = False
                for file in os.listdir(self.databaseDirectory):
                    if (file.find(".cspr") != -1):
                        foundCSPRFiles = True
                        break

                #if CSPR files found in database directory
                if foundCSPRFiles == True:

                    #change directory to database directory
                    os.chdir(self.databaseDirectory)

                    #update database directory global variable
                    GlobalSettings.CSPR_DB = self.databaseDirectory

                    #save database directory to CASPERinfo
                    self.saveDatabaseDirectory()

                    # make sure FNA and GBFF subdirectories are present
                    subdirs = os.listdir(self.databaseDirectory)
                    if "FNA" not in subdirs and os.path.isdir("FNA") == False:
                        try:
                            os.mkdir("FNA")
                            logger.debug("Successfully created FNA subdirectory in database directory.")
                        except Exception as e:
                            logger.critical("Unable to make 'FNA' subdirectory in database directory")
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
                    if "GBFF" not in subdirs and os.path.isdir("GBFF") == False:
                        try:
                            os.mkdir("GBFF")
                            logger.debug("Successfully created GBFF subdirectory in database directory.")
                        except Exception as e:
                            logger.critical("Unable to make 'GBFF' subdirectory in database directory")
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

                    ### fill in organism/endo/annotation dropdown information for main, mulit-targeting, and populatin analysis
                    try:
                        GlobalSettings.mainWindow.getData()
                        GlobalSettings.mainWindow.fill_annotation_dropdown()
                        logger.debug("Successfully loaded organism/endo/annotation drop down information in Main.")
                    except Exception as e:
                        logger.critical("Unable to load organism/endo/annotation drop down information in Main.")
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

                    try:
                        GlobalSettings.MTWin.launch()
                        logger.debug("Successfully loaded organism/endo drop down information in Multi-targeting.")
                    except Exception as e:
                        logger.critical("Unable to load organism/endo drop down information in Multi-targeting.")
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

                    try:
                        GlobalSettings.pop_Analysis.launch()
                        logger.debug("Successfully loaded organism/endo drop down information in Population Analysis.")
                    except Exception as e:
                        logger.critical("Unable to load organism/endo drop down information in Population Analysis.")
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

                    #show main window
                    if GlobalSettings.mainWindow.first_show == True:
                        GlobalSettings.mainWindow.show()
                        GlobalSettings.mainWindow.first_show = False
                    else:
                        GlobalSettings.mainWindow.show()
                    self.close()

                #no cspr file found
                else:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("Directory is invalid!")
                    msgBox.setText(
                        "You must select a directory with CSPR Files!")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    return
            #not a directory
            else:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Not a directory")
                msgBox.setText(
                    "The directory you selected does not exist.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return
        except Exception as e:
            logger.critical("Error in launchMain() in startup window.")
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


#initial function called during startup
def main():
    #log OS
    logger.info("System OS: %s" % (platform.system()))

    if hasattr(sys, 'frozen'):
        #log CASPER is in packaged format
        logger.info("Running a packaged version of CASPER.")

        GlobalSettings.appdir = sys.executable
        if platform.system() == 'Windows':
            GlobalSettings.appdir = GlobalSettings.appdir[:GlobalSettings.appdir.rfind("\\") + 1]
        else: # If platform is Mac
            GlobalSettings.appdir = GlobalSettings.appdir[:GlobalSettings.appdir.rfind("Contents/") + 9]
            GlobalSettings.appdir += "Resources/" 
    else:
        # log CASPER is not in packaged format
        logger.info("Running a non-packaged version of CASPER.")

        GlobalSettings.appdir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            GlobalSettings.appdir += '\\'
        else:
            GlobalSettings.appdir += '/'

    #QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)
    #QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_UseHighDpiPixmaps, True)

    app = QtWidgets.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")

    #setup logger
    fh = logging.FileHandler(GlobalSettings.appdir + 'CASPER.log', mode='w')
    fh_formatter = logging.Formatter('%(asctime)s %(levelname)s %(lineno)d:%(filename)s(%(process)d) - %(message)s')
    fh.setFormatter(fh_formatter)
    fh.setLevel(logging.DEBUG)
    GlobalSettings.logger.addHandler(fh)

    #log appdir
    logger.debug("App Directory: " + str(GlobalSettings.appdir))

    #load startup window
    try:
        startup = StartupWindow()
        logger.debug("Successfully initialized Startup Window.")
    except Exception as e:
        logger.critical("Can't start Startup window.")
        logger.critical(e)
        logger.critical(traceback.format_exc())
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet("font: " + str(startup.fontSize) + "pt 'Arial'")
        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
        msgBox.setWindowTitle("Fatal Error")
        msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
        msgBox.exec()

        exit(-1)

    #load main
    try:
        GlobalSettings.mainWindow = CMainWindow(os.getcwd())
        logger.debug("Successfully initialized Main Window.")
    except Exception as e:
        logger.critical("Can't start Main window.")
        logger.critical(e)
        logger.critical(traceback.format_exc())
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet("font: " + str(startup.fontSize) + "pt 'Arial'")
        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
        msgBox.setWindowTitle("Fatal Error")
        msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
        msgBox.exec()

        exit(-1)

    #load multi-targeting
    try:
        GlobalSettings.MTWin = multitargeting.Multitargeting()
        logger.debug("Successfully initialized Multi-targeting Window.")
    except Exception as e:
        logger.critical("Can't start Multi-targeting window.")
        logger.critical(e)
        logger.critical(traceback.format_exc())
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet("font: " + str(GlobalSettings.MTWin.fontSize) + "pt 'Arial'")
        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
        msgBox.setWindowTitle("Fatal Error")
        msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
        msgBox.exec()

        exit(-1)

    #load pop analysis
    try:
        GlobalSettings.pop_Analysis = populationAnalysis.Pop_Analysis()
        logger.debug("Successfully initialized Population Analysis Window.")
    except Exception as e:
        logger.critical("Can't start Population Analysis window.")
        logger.critical(e)
        logger.critical(traceback.format_exc())
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet("font: " + str(GlobalSettings.pop_Analysis.fontSize) + "pt 'Arial'")
        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
        msgBox.setWindowTitle("Fatal Error")
        msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
        msgBox.exec()

        exit(-1)

    startup.centerUI()
    startup.show()

    #GlobalSettings.MTWin.centerUI()
    #GlobalSettings.MTWin.show()

    sys.exit(app.exec_())


#call initial startup function main()
if __name__ == '__main__':
    main()
