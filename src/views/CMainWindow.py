import platform
# import controllers.ncbi as ncbi
import os
# from utils.Algorithms import get_table_headers
# from models.CSPRparser import CSPRparser
import glob
import models.GlobalSettings as GlobalSettings
from PyQt6 import QtWidgets, QtGui, QtCore, uic, QtGui
from utils.ui import scale_ui, center_ui, show_message, show_error
# from views.annotation_functions import *
# from views.AnnotationParser import Annotation_Parser
# from views.AnnotationWindow import AnnotationWindow
# import views.genomeBrowser as genomeBrowser
# from views.NewGenome import NewGenome
# from views.NewEndonuclease import NewEndonuclease
# from controllers.CoTargeting import CoTargeting
# from views.generateLib import genLibrary
# from controllers.Results import Results
# from views.export_tool import export_tool
# from views.closingWin import closingWindow
# from utils.web import ncbi_page, repo_page, ncbi_blast_page
# from controllers.populate_fna_files import PopulateFNAFiles

# logger = GlobalSettings.logger

fontSize = 12

class CMainWindow(QtWidgets.QMainWindow):
    def __init__(self, settings):
        try:
            super(CMainWindow, self).__init__()
            # uic.loadUi(os.path.join(self.settings.get_ui_dir(), 'startupCASPER.ui'), self)
            # uic.loadUi(GlobalSettings.appdir + 'ui/CASPER_main.ui', self)
            # print("path: ", GlobalSettings.appdir + 'ui/CASPER_main_copy_2.ui')
            # uic.loadUi(GlobalSettings.appdir + 'ui/CASPER_main_copy_2.ui', self)
            self.settings = settings
            print("path: ", os.path.join(self.settings.get_ui_dir(), 'CASPER_main.ui'))
            uic.loadUi(os.path.join(self.settings.get_ui_dir(), 'CASPER_main.ui'), self)
            self.setWindowTitle("CASPER")
            self.setWindowIcon(QtGui.QIcon(os.path.join(self.settings.get_assets_dir(), "cas9image.ico")))

            # self.dbpath = ""
            # self.inputstring = "" # This is the search string
            # # self.info_path = settings.get_app_dir()
            # # info_path = settings.get_app_dir()
            # self.anno_name = ""
            # self.endo_name = ""
            # self.fontSize = 12
            # self.org = ""
            # self.TNumbers = {}  # the T numbers from a kegg search
            # self.orgcodes = {}  # Stores the Kegg organism code by the format {full name : organism code}
            # self.gene_list = {}  # list of genes (no ides what they pertain to
            # self.searches = {}
            # self.checkBoxes = []
            # self.genlib_list = [] # This list stores selected SeqFeatures from annotation window
            # self.checked_info = {}
            # self.check_ntseq_info = {}  # the ntsequences that go along with the checked_info
            # self.annotation_parser = Annotation_Parser()
            # self.link_list = list()  # the list of the downloadable links from the NCBI search
            # self.organismDict = dict()  # the dictionary for the links to download. Key is the description of the organism, value is the ID that can be found in link_list
            # self.results_list = list()
            # self.organismData = list()
            # self.ncbi = ncbi.NCBI_search_tool()

            # groupbox_style = """
            # QGroupBox:title{subcontrol-origin: margin;
            #                 left: 10px;
            #                 padding: 0 5px 0 5px;}
            # QGroupBox#Step1{border: 2px solid rgb(111,181,110);
            #                 border-radius: 9px;
            #                 margin-top: 10px;
            #                 font: bold 14pt 'Arial';}
            #                 """

            # self.Step1.setStyleSheet(groupbox_style)
            # self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2"))
            # self.Step3.setStyleSheet(groupbox_style.replace("Step1", "Step3"))
            # self.CASPER_Navigation.setStyleSheet(groupbox_style.replace("Step1", "CASPER_Navigation").replace("solid","dashed").replace("rgb(111,181,110)","rgb(88,89,91)"))

            # self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            # self.pushButton_FindTargets.clicked.connect(self.gather_settings)
            # self.pushButton_ViewTargets.clicked.connect(self.view_results)
            # self.pushButton_ViewTargets.setEnabled(False)
            # self.GenerateLibrary.setEnabled(False)
            # self.radioButton_Gene.clicked.connect(self.toggle_annotation)
            # self.radioButton_Position.clicked.connect(self.toggle_annotation)
        
            """ Connect functions to buttons """ 
            # self.newGenome_button.clicked.connect(self.launch_newGenome) # Connect launch function to New Genome
            # self.newEndo_button.clicked.connect(self.launch_newEndonuclease) # Connect launch function to New Endonuclease
            # self.multitargeting_button.clicked.connect(self.changeto_multitargeting) # Connect launch function to Multitargeting
            # self.populationAnalysis_button.clicked.connect(self.changeto_population_Analysis) # Connect launch function to PA
            # self.GenerateLibrary.clicked.connect(self.prep_genlib)
            # self.combineFiles_button.clicked.connect(self.launch_populate_fna_files)

            """ Connect functions to actions (menu bar) """ 
            # self.actionOpen_Genome_Browser.triggered.connect(self.launch_newGenomeBrowser)
            # self.actionExit.triggered.connect(self.close_app)
            # self.visit_repo.triggered.connect(repo_page)
            # self.actionChange_Directory.triggered.connect(self.change_directory)
            # self.actionNCBI.triggered.connect(ncbi_page)
            # self.actionCasper2.triggered.connect(self.open_casper2_web_page)
            # self.actionNCBI_BLAST.triggered.connect(ncbi_blast_page)



            # self.progressBar.setMinimum(0)
            # self.progressBar.setMaximum(100)
            # self.progressBar.reset()
            # self.Annotation_Window = AnnotationWindow(info_path)
            # self.geneEntryField.setPlaceholderText("Example Inputs: \n\n"
            #                                 "Option 1: Feature (ID, Locus Tag, or Name)\n"
            #                                 "Example: 854068/YOL086C/ADH1 for S. cerevisiae alcohol dehydrogenase 1\n\n"
            #                                 "Option 2: Position (chromosome,start,stop)\n"
            #                                 "Example: 1,1,1000 for targeting chromosome 1, base pairs 1 to 1000\n\n"
            #                                 "Option 3: Sequence (must be within the selected organism)\n"
            #                                 "Example: Any nucleotide sequence between 100 and 10,000 base pairs.\n\n"
            #                                 "*Note: to multiplex, separate multiple queries by new lines*\n"
            #                                 "Example:\n"
            #                                 "1,1,1000\n"
            #                                 "5,1,500\n"
            #                                 "etc.")

            # show functionalities on window
            self.populate_fna_files = None 
            self._new_genome = None
            # self.newEndonuclease = NewEndonuclease()
            # self.CoTargeting = CoTargeting(info_path)
            # self.Results = Results()
            # self.export_tool_window = export_tool()
            # self.genLib = genLibrary()
            # self.myClosingWindow = closingWindow()
            # self.genomebrowser = genomeBrowser.genomebrowser()
            # self.launch_ncbi_button.clicked.connect(self.launch_ncbi)

            # self.first_show = True
            scale_ui(self, custom_scale_width=1150, custom_scale_height=650)
            # self.show()
            # self.load_dropdown_data()
            print("MainWindow initialized")
        except Exception as e:
            show_error("Error in __init__() in main", e)

    # def get_populate_fna_files(self):
    #     if self.populate_fna_files is None:
    #         self.populate_fna_files = PopulateFNAFiles(GlobalSettings.GlobalSettings1(GlobalSettings.appdir))
    #     return self.populate_fna_files
    
    # def launch_populate_fna_files(self):
    #     self.get_populate_fna_files().show()  # Ensure the window is shown

    # this function prepares everything for the generate library function
    # it is very similar to the gather settings, how ever it stores the data instead of calling the Annotation Window class
    # it moves the data onto the generateLib function, and then opens that window
    # def prep_genlib(self):
    #     # make sure the user actually inputs something
    #     try:
    #         inputstring = str(self.geneEntryField.toPlainText())
    #         if (inputstring.startswith("Example Inputs:") or inputstring == ""):
    #             show_message(
    #                 fontSize=12,
    #                 icon=QtWidgets.QMessageBox.Icon.Critical,
    #                 title="Error",
    #                 message="No gene has been entered.  Please enter a gene.",
    #                 button=QtWidgets.QMessageBox.StandardButton.Ok
    #             )
    #             return
    #         else:
    #             # standardize the input
    #             inputstring = inputstring.lower()
    #             found_matches_bool = True
    #             # call the respective function
    #             self.progressBar.setValue(10)
    #             if self.radioButton_Gene.isChecked():
    #                 if len(self.genlib_list) > 0:
    #                     found_matches_bool = True
    #                 else:
    #                     found_matches_bool = False
    #             elif self.radioButton_Position.isChecked() or self.radioButton_Sequence.isChecked():
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Critical,
    #                     title="Error",
    #                     message="Generate Library can only work with feature searches.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Ok
    #                 )
    #                 return
    #             """
    #             elif self.radioButton_Position.isChecked():
    #                 pinput = inputstring.split(';')
    #                 found_matches_bool = self.run_results("position", pinput,openAnnoWindow=False)
    #             elif self.radioButton_Sequence.isChecked():
    #                 sinput = inputstring
    #                 found_matches_bool = self.run_results("sequence", sinput, openAnnoWindow=False)
    #             """
    #             # if matches are found
    #             if found_matches_bool == True:
    #                 # get the cspr file name
    #                 cspr_file = self.organisms_to_files[self.orgChoice.currentText()][self.endoChoice.currentText()][0]
    #                 if platform.system() == 'Windows':
    #                     cspr_file = GlobalSettings.CSPR_DB + '\\' + cspr_file
    #                 else:
    #                     cspr_file = GlobalSettings.CSPR_DB + '/' + cspr_file
    #                 kegg_non = 'non_kegg'

    #                 # launch generateLib
    #                 self.progressBar.setValue(100)

    #                 # calculate the total number of matches found
    #                 tempSum = len(self.genlib_list)

    #                 # warn the user if the number is greater than 50
    #                 if tempSum > 50:
    #                     msgBox = QtWidgets.QMessageBox()
    #                     msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
    #                     msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
    #                     msgBox.setWindowTitle("Many Matches Found")
    #                     msgBox.setText("More than 50 matches have been found. Continuing could cause a slow down...\n\n Do you wish to continue?")
    #                     msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
    #                     msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
    #                     msgBox.exec()

    #                     if (msgBox.result() == QtWidgets.QMessageBox.No):
    #                         self.searches.clear()
    #                         self.progressBar.setValue(0)
    #                         return -2

    #                 self.genLib.launch(self.genlib_list,cspr_file, kegg_non)
    #             else:
    #                 self.progressBar.setValue(0)
    #     except Exception as e:
    #         show_error("Error in prep_genlib() in main", e)
         
    # # Function for collecting the settings from the input field and transferring them to run_results
    # def gather_settings(self):
    #     try:
    #         ### If user searches multiple times for the same thing, this avoids re-searching the entire annotation file
    #         check_org = self.orgChoice.currentText().lower()
    #         check_endo = self.endoChoice.currentText().lower()
    #         check_anno_name = self.annotation_files.currentText().lower()
    #         check_input = str(self.geneEntryField.toPlainText()).lower()
    #         if (check_input == self.inputstring and check_org == self.org and check_anno_name == self.anno_name and check_endo == self.endo_name):
    #             same_search = True
    #         else:
    #             self.org = check_org
    #             self.anno_name = check_anno_name 
    #             self.inputstring = check_input
    #             self.endo_name = check_endo
    #             same_search = False

    #         # Error check: make sure the user actually inputs something
    #         if (self.inputstring.startswith("Example Inputs:") or self.inputstring == ""):
    #             show_message(
    #                 fontSize=12,
    #                 icon=QtWidgets.QMessageBox.Icon.Critical,
    #                 title="Error",
    #                 message="No feature has been searched for. Please enter a search.",
    #                 button=QtWidgets.QMessageBox.StandardButton.Ok
    #             )
    #             return
    #         else:

    #             ### Remove additional scoring columns if necessary
    #             header = get_table_headers(self.Results.targetTable) # Returns headers of the target table in View Targets window 
    #             col_indices = [header.index(x) for x in GlobalSettings.algorithms if x in header] # Returns the index(es) of the alternative scoring column(s) in the target table of View Targets window
    #             if len(col_indices) > 0: # If alternative scoring has been done
    #                 for i in col_indices:
    #                     self.Results.targetTable.removeColumn(i)
    #             self.Results.targetTable.resizeColumnsToContents()

    #             self.progressBar.setValue(10)
    #             if self.radioButton_Gene.isChecked():
    #                 ginput = [x.strip() for x in self.inputstring.split('\n')] # Split search based on newline character and remove deadspace
    #                 self.run_results("feature", ginput, same_search)
    #             elif self.radioButton_Position.isChecked():
    #                 pinput = [x.strip() for x in self.inputstring.split('\n')] # Split search based on newline character and remove deadspace
    #                 self.run_results("position", pinput, same_search)
    #             elif self.radioButton_Sequence.isChecked():
    #                 sinput = self.inputstring
    #                 self.run_results("sequence", sinput, same_search)
    #     except Exception as e:
    #         show_error("Error in gather_settings() in main", e)
            
    # # ---- Following functions are for running the auxillary algorithms and windows ---- #
    # # this function is parses the annotation file given, and then goes through and goes onto results
    # # it will call other versions of collect_table_data and fill_table that work with these file types
    # # this function should work with the any type of annotation file, besides kegg.
    # # this assumes that the parsers all store the data the same way, which gff and feature table do
    # # please make sure the genbank parser stores the data in the same way
    # # so far the gff files seems to all be different. Need to think about how we want to parse it
    # def run_results_own_ncbi_file(self, inputstring, fileName, same_search, openAnnoWindow=True):
    #     try:
    #         self.set_progress(35)
    #         self.results_list = self.annotation_parser.genbank_search(inputstring, same_search)

    #         cspr_file = self.organisms_to_files[self.orgChoice.currentText()][self.endoChoice.currentText()][0]
    #         cspr_file = os.path.join(GlobalSettings.CSPR_DB, cspr_file)
            
    #         own_cspr_parser = CSPRparser(cspr_file)
    #         own_cspr_parser.read_first_lines()
    #         if len(own_cspr_parser.karystatsList) != self.annotation_parser.max_chrom:
    #             show_message(
    #                 fontSize=12,
    #                 icon=QtWidgets.QMessageBox.Icon.Warning,
    #                 title="Warning:",
    #                 message="The number of chromosomes do not match. This could cause errors.",
    #                 button=QtWidgets.QMessageBox.StandardButton.Ok
    #             )
    #         self.set_progress(60)

    #         self.searches.clear()

    #         self.set_progress(75)
    #         if not self.results_list:
    #             show_message(
    #                 fontSize=12,
    #                 icon=QtWidgets.QMessageBox.Icon.Critical,
    #                 title="No Matches Found",
    #                 message="No matches found with that search, please try again.",
    #                 button=QtWidgets.QMessageBox.StandardButton.Ok
    #             )
    #             self.set_progress(0) 
    #             return False if not openAnnoWindow else None

    #         self.set_progress(80)

    #         return self.Annotation_Window.fill_table_nonKegg(self, self.results_list) if openAnnoWindow else True
    #     except Exception as e:
    #         show_error(f"Error in run_results_own_ncbi_file() in main.", e)

    # def set_progress(self, value):
    #     self.progressBar.setValue(value)

    # def run_results(self, inputtype, inputstring, same_search, openAnnoWindow=True):
    #     try:
    #         file_name = self.annotation_files.currentText()
    #         for file in glob.glob(GlobalSettings.CSPR_DB + "/**/*.gb*", recursive=True):
    #             if file_name in file:
    #                 self.annotation_parser.annotationFileName = file
    #                 break
    #         self.Results.annotation_path = self.annotation_parser.annotationFileName

    #         progvalue = 15
    #         self.searches = {}
    #         self.gene_list = {}
    #         self.progressBar.setValue(progvalue)

    #         try:
    #             self.Results.endonucleaseBox.currentIndexChanged.disconnect()
    #         except Exception as e:
    #             pass
    #         # set Results endo combo box
    #         self.Results.endonucleaseBox.clear()

    #         # set the results window endoChoice box menu
    #         # set the mainWindow's endoChoice first, and then loop through and set the rest of them
    #         self.Results.endonucleaseBox.addItem(self.endoChoice.currentText())
    #         for item in self.organisms_to_endos[str(self.orgChoice.currentText())]:
    #             if item != self.Results.endonucleaseBox.currentText():
    #                 self.Results.endonucleaseBox.addItem(item)

    #         self.Results.endonucleaseBox.currentIndexChanged.connect(self.Results.changeEndonuclease)
    #         self.Results.get_endo_data()

    # #        self.Results.change_start_end_button.setEnabled(False)
    #         self.Results.displayGeneViewer.setChecked(0)

    #         if inputtype == "feature":
    #             fileType = self.annotation_parser.find_which_file_version()

    #             # if the parser retuns the 'wrong file type' error
    #             if fileType == -1:
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Critical,
    #                     title="Error",
    #                     message="Feature search requires a GenBank formatted annotation file. Please select a file from the dropdown menu or search by position",
    #                     button=QtWidgets.QMessageBox.StandardButton.Ok
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return

    #             # make sure an annotation file has been selected
    #             if self.annotation_files.currentText() == "None":
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Critical,
    #                     title="Error",
    #                     message="Search by feature requires a GenBank annotation file. Please select one from the dropdown menu or search by position.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Ok
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return

    #             # this now just goes onto the other version of run_results
    #             myBool = self.run_results_own_ncbi_file(inputstring, self.annotation_files.currentText(), same_search, openAnnoWindow=openAnnoWindow)
    #             if not openAnnoWindow:
    #                 return myBool
    #             else:
    #                 self.progressBar.setValue(0)
    #                 return

    #         if inputtype == "position":
    #             full_org = str(self.orgChoice.currentText())
    #             self.checked_info.clear()
    #             self.check_ntseq_info.clear()

    #             for item in inputstring: # Loop through each search
    #                 searchIndices = [x.strip() for x in item.split(',')] # Parse input query

    #                 if len(searchIndices) != 3:
    #                     show_message(
    #                         fontSize=12,
    #                         icon=QtWidgets.QMessageBox.Icon.Critical,
    #                         title="Error",
    #                         message="There are 3 arguments required for this function: chromosome, start position, and end position.",
    #                         button=QtWidgets.QMessageBox.StandardButton.Ok
    #                     )
    #                     self.progressBar.setValue(0)
    #                     return

    #                 if not searchIndices[0].isdigit() or not searchIndices[1].isdigit() or not searchIndices[2].isdigit():
    #                     show_message(
    #                         fontSize=12,
    #                         icon=QtWidgets.QMessageBox.Icon.Critical,
    #                         title="Error",
    #                         message="The positions given must be integers. Please try again.",
    #                         button=QtWidgets.QMessageBox.StandardButton.Ok
    #                     )
    #                     self.progressBar.setValue(0)
    #                     return
    #                 elif int(searchIndices[1]) >= int(searchIndices[2]):
    #                     show_message(
    #                         fontSize=12,
    #                         icon=QtWidgets.QMessageBox.Icon.Critical,
    #                         title="Error",
    #                         message="The start index must be less than the end index.",
    #                         button=QtWidgets.QMessageBox.StandardButton.Ok
    #                     )
    #                     self.progressBar.setValue(0)
    #                     return
    #                 elif abs(int(searchIndices[2])-int(searchIndices[1])) > 50000:
    #                     show_message(
    #                         fontSize=12,
    #                         icon=QtWidgets.QMessageBox.Icon.Critical,
    #                         title="Error",
    #                         message="The search range must be less than 50,000 nt.",
    #                         button=QtWidgets.QMessageBox.StandardButton.Ok
    #                     )
    #                     self.progressBar.setValue(0)
    #                     return
    #                 elif int(searchIndices[0]) > self.annotation_parser.get_max_chrom():
    #                     show_message(
    #                         fontSize=12,
    #                         icon=QtWidgets.QMessageBox.Icon.Critical,
    #                         title="Error",
    #                         message="Chromosome %s does not exist in the selected annotation file." % searchIndices[0],
    #                         button=QtWidgets.QMessageBox.StandardButton.Ok
    #                     )
    #                     self.progressBar.setValue(0)
    #                     return
    #                 # append the data into the checked_info
    #                 tempString = 'chrom: ' + str(searchIndices[0]) + ',start: ' + str(searchIndices[1]) + ',end: ' + str(searchIndices[2])
    #                 self.checked_info[tempString] = (int(searchIndices[0]), int(searchIndices[1])-1, int(searchIndices[2]))

    #             self.progressBar.setValue(50)
    #             self.Results.transfer_data(full_org, self.organisms_to_files[full_org], [str(self.endoChoice.currentText())], os.getcwd(), self.checked_info, self.check_ntseq_info,inputtype)
    #             self.Results.load_gene_viewer()
    #             self.progressBar.setValue(100)
    #             self.pushButton_ViewTargets.setEnabled(True)
    #             self.GenerateLibrary.setEnabled(True)

    #         if inputtype == "sequence":
    #             fileType = self.annotation_parser.find_which_file_version()

    #             if fileType == -1:
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Critical,
    #                     title="Error",
    #                     message="Search by sequence requires a GenBank annotation file. Please select one from the dropdown menu or search by position.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Ok
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return
    #             if self.annotation_files.currentText() == "None":
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Critical,
    #                     title="Error",
    #                     message="Search by sequence requires a GenBank annotation file. Please select one from the dropdown menu or search by position.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Ok
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return

    #             checkString = 'AGTCN'
    #             full_org = str(self.orgChoice.currentText())
    #             self.checked_info.clear()
    #             self.progressBar.setValue(10)
    #             inputstring = inputstring.replace('\n','').upper().strip()

    #             for letter in inputstring:
    #                 if letter not in checkString:
    #                     show_message(
    #                         fontSize=12,
    #                         icon=QtWidgets.QMessageBox.Icon.Critical,
    #                         title="Error",
    #                         message="The sequence must consist of A, G, T, C, or N. No other characters are allowed.",
    #                         button=QtWidgets.QMessageBox.StandardButton.Ok
    #                     )
    #                     self.progressBar.setValue(0)
    #                     return

    #             if len(inputstring) < 100:
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Critical,
    #                     title="Error",
    #                     message="The sequence given is too small. At least 100 characters are required.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Ok
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return

    #             if len(inputstring) > 10000:
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Question,
    #                     title="Large Sequence Detected",
    #                     message="The sequence given is too large one.\n\nPlease input a sequence less than 10kb in length.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Yes
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return

    #             self.progressBar.setValue(30)

    #             # Check the GBFF file for the sequence
    #             my_check = self.annotation_parser.get_sequence_info(inputstring)

    #             self.progressBar.setValue(55) # Update progress bar

    #             if type(my_check) == bool: 
    #                 show_message(
    #                     fontSize=12,
    #                     icon=QtWidgets.QMessageBox.Icon.Question,
    #                     title="Sequence Not Found",
    #                     message="The sequence entered was not found.\n\nPlease input a sequence that is in the selected organism.",
    #                     button=QtWidgets.QMessageBox.StandardButton.Yes
    #                 )
    #                 self.progressBar.setValue(0)
    #                 return

    #             else: 
    #                 tempString = 'chrom: ' + str(my_check[0]) + ',start: ' + str(my_check[1]) + ',end: ' + str(my_check[2])
    #                 self.checked_info[tempString] = (int(my_check[0]), int(my_check[1])-1, int(my_check[2]))

    #             self.progressBar.setValue(75) 

    #             self.Results.transfer_data(full_org, self.organisms_to_files[full_org], [str(self.endoChoice.currentText())], os.getcwd(), self.checked_info, self.check_ntseq_info, inputtype)
    #             self.Results.load_gene_viewer()
    #             self.progressBar.setValue(100)
    #             self.pushButton_ViewTargets.setEnabled(True)
    #             self.GenerateLibrary.setEnabled(True)
    #     except Exception as e:
    #         show_error("Error in run_results() in main", e)
            
    # def handle_feature_search(self, input_string, open_anno_window):
    #     file_type = self.annotation_parser.find_which_file_version()
    #     if file_type == -1 or self.annotation_files.currentText() == "None":
    #         self.show_error_message("Feature search requires a GenBank formatted annotation file.")
    #         return False

    #     return self.run_results_own_ncbi_file(input_string, self.annotation_files.currentText(), same_search, open_anno_window)

    # def launch_newGenome(self):
    #     try:
    #         # Update endo list
    #         self.get_new_genome().fillEndo()
    #         if self.get_new_genome().first_show:
    #             center_ui(self.get_new_genome())
    #             self.get_new_genome().first_show = False
    #         self.hide()
    #         self.get_new_genome().show()
    #     except Exception as e:
    #         show_error("Error in launch_newGenome() in main", e)

    # def launch_newEndonuclease(self):
    #     try:
    #         center_ui(self.newEndonuclease)
    #         self.newEndonuclease.show()
    #         self.newEndonuclease.activateWindow()
    #     except Exception as e:
    #         show_error("Error in launch_newEndonuclease() in main", e)

    # #launch genome browser tool
    # def launch_newGenomeBrowser(self):
    #     try:
    #         self.genomebrowser.createGraph(self)
    #     except Exception as e:
    #         show_error("Error in launch_newGenomeBrowser() in main", e)

    # def launch_ncbi(self):
    #     try:
    #         show_message(
    #             fontSize=12,
    #             icon=QtWidgets.QMessageBox.Icon.Information,
    #             title="Note:",
    #             message="NCBI Annotation Guidelines:\n\nDownload annotation files of the exact species and strain used in Analyze New Genome.\n\nMismatched annotation files will inhibit downstream analyses.",
    #             button=QtWidgets.QMessageBox.StandardButton.Ok
    #         )
    #         if self.ncbi.first_show:
    #             self.ncbi.first_show = False
    #             center_ui(self.ncbi)

    #         self.ncbi.show()
    #         self.ncbi.activateWindow()
    #     except Exception as e:
    #         show_error("launch_ncbi() in main", e)

    # # this function does the same stuff that the other collect_table_data does, but works with the other types of files
    # def collect_table_data_nonkegg(self):
    #     try:
    #         # start out the same as the other collect_table_data
    #         self.checked_info.clear()
    #         self.genlib_list.clear()
    #         self.check_ntseq_info.clear()
    #         full_org = str(self.orgChoice.currentText())
    #         holder = ()
    #         selected_indices = []
    #         selected_rows = self.Annotation_Window.tableWidget.selectionModel().selectedRows()
    #         for ind in sorted(selected_rows):
    #             selected_indices.append(ind.row())

    #         for item in self.checkBoxes:
    #             feature = item[1]
    #             # If inidices of checkBoxes list and selected rows in table match...
    #             if item[2] in selected_indices:
    #                 holder = (item[0],int(feature.location.start),int(feature.location.end)) # Tuple order: Feature chromosome/scaffold number, feature start, feature end
    #                 ### If locus tag available, combine with gene name to create dict key
    #                 if 'locus_tag' in feature.qualifiers:
    #                     tag = feature.qualifiers['locus_tag'][0]
    #                     key = tag + ": " + get_name(feature)
    #                 else:
    #                     key = get_name(feature)
    #                 self.checked_info[key] = holder
    #                 self.genlib_list.append((item[0],feature)) # Tuple order: Feature chromosome/scaffold number, SeqFeature object
    #             else:
    #                 # If item was not selected in the table, go to the next item
    #                 continue

    #         # now call transfer data
    #         self.progressBar.setValue(95)
    #         self.Results.transfer_data(full_org, self.organisms_to_files[full_org], [str(self.endoChoice.currentText())], os.getcwd(),
    #                                    self.checked_info, self.check_ntseq_info,inputtype="feature")
    #         self.Results.load_gene_viewer()

    #         self.progressBar.setValue(100)
    #         self.pushButton_ViewTargets.setEnabled(True)
    #         self.GenerateLibrary.setEnabled(True)
    #     except Exception as e:
    #         show_error("Error in collect_table_data_nonkegg() in main", e)
            
    # def separate_line(self, input_string):
    #     try:
    #         export_array = []
    #         while True:
    #             index = input_string.find('\n')
    #             if index == -1:
    #                 if len(input_string) == 0:
    #                     return export_array
    #                 else:
    #                     export_array.append(input_string)
    #                     return export_array
    #             export_array.append(input_string[:index])
    #             input_string = input_string[index + 1:]
    #     except Exception as e:
    #         show_error("Error in seperate_line() in main", e)
           
    # def removeWhiteSpace(self, strng):
    #     try:
    #         while True:
    #             if len(strng) == 0 or (strng[0] != " " and strng[0] != "\n"):
    #                 break
    #             strng = strng[1:]
    #         while True:
    #             if len(strng) == 0 or (strng[len(strng) - 1] != " " and strng[0] != "\n"):
    #                 return strng
    #             strng = strng[:len(strng) - 1]
    #     except Exception as e:
    #         show_error("Error in removeWhiteSpace() in main", e)

    # # Function to enable and disable the Annotation function if searching by position or sequence
    # def toggle_annotation(self):
    #     try:
    #         if self.radioButton_Gene.isChecked():
    #             self.Step2.setEnabled(True)
    #         else:
    #             self.Step2.setEnabled(True)
    #     except Exception as e:
    #         show_error("Error in toggle_annotation() in main", e)

    # def fill_annotation_dropdown(self):
    #     try:
    #         #recursive search for all GenBank files in casper db folder
    #         self.annotation_files.clear()
    #         annotation_files = glob.glob(GlobalSettings.CSPR_DB + "/**/*.gb*", recursive=True)
    #         if platform.system() == "Windows":
    #             for i in range(len(annotation_files)):
    #                 annotation_files[i] = annotation_files[i].replace("/","\\")
    #                 annotation_files[i] = annotation_files[i][annotation_files[i].rfind("\\") + 1:]
    #         else:
    #             for i in range(len(annotation_files)):
    #                 annotation_files[i] = annotation_files[i].replace("\\","/")
    #                 annotation_files[i] = annotation_files[i][annotation_files[i].rfind("/") + 1:]

    #         annotation_files.sort(key=str.lower)
    #         self.annotation_files.addItems(annotation_files)
    #         self.annotation_files.addItems(["None"])
    #     except Exception as e:
    #         show_error("Error in fill_annotation_dropdown() in main", e)
            
    # def make_dictonary(self):
    #     try:
    #         url = "https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:" + self.TNumbers[
    #             self.Annotations_Organism.currentText()]
    #         source_code = requests.get(url, verify=False)
    #         plain_text = source_code.text
    #         buf = io.StringIO(plain_text)

    #         while True:
    #             line = buf.readline()
    #             if line[0] == "-":
    #                 break
    #         while True:
    #             line = buf.readline()
    #             if line[1] != "a":
    #                 return
    #             line = line[line.find(">") + 1:]
    #             seq = line[line.find(":") + 1:line.find("<")]
    #             line = line[line.find(">") + 1:]

    #             i = 0
    #             while True:
    #                 if line[i] == " ":
    #                     i = i + 1
    #                 else:
    #                     break
    #             key = line[i:line.find("\n") - 1]
    #             if key in self.gene_list:
    #                 if seq not in self.gene_list[key]:
    #                     self.gene_list[key].append(seq)
    #             else:
    #                 self.gene_list[key] = [seq]
    #             z = 5
    #     except Exception as e:
    #         show_error("Error in make_dictionary() in main", e)

    # def organism_finder(self, long_str):
    #     try:
    #         semi = long_str.find(";")
    #         index = 1
    #         while True:
    #             if long_str[semi - index] == " ":
    #                 break
    #             index = index + 1
    #         return long_str[:semi - index]
    #     except Exception as e:
    #         show_error("Error in organism_finder() in main", e)

    # # This method is for testing the execution of a button call to make sure the button is linked properly
    # def testexe(self):
    #     try:
    #         msgBox = QtWidgets.QMessageBox()
    #         msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
    #         msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
    #         msgBox.setWindowTitle("Extract!")
    #         msgBox.setText(
    #             "Are you sure you want to quit?")
    #         msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
    #         msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
    #         msgBox.exec()

    #         if msgBox.result() == QtWidgets.QMessageBox.Yes:
    #             # print(self.orgChoice.currentText())
    #             sys.exit()
    #         else:
    #             pass
    #     except Exception as e:
    #         show_error("Error in testexe() in main", e)
            
    # def getData(self):
    #     try:
    #         try:
    #             self.orgChoice.currentIndexChanged.disconnect()
    #         except Exception as e:
    #              pass

    #         self.orgChoice.clear()
    #         self.endoChoice.clear()
    #         mypath = os.getcwd()
    #         found = False
    #         self.dbpath = mypath
    #         onlyfiles = [str(f) for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    #         onlyfiles.sort(key=str.lower)
    #         self.organisms_to_files = {}
    #         self.organisms_to_endos = {}
    #         first = True
    #         for file in onlyfiles:
    #             if file.find('.cspr') != -1:
    #                 if first == True:
    #                     first = False
    #                 found = True
    #                 newname = file[0:-4]
    #                 endo = newname[newname.rfind("_")+1:-1]
    #                 hold = open(file, 'r')
    #                 buf = (hold.readline())
    #                 buf = str(buf)
    #                 buf = buf.strip()
    #                 species = buf.replace("GENOME: ",'')

    #                 if species in self.organisms_to_files:
    #                     self.organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]
    #                 else:
    #                     self.organisms_to_files[species] = {}
    #                     self.organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]

    #                 if species in self.organisms_to_endos:
    #                     self.organisms_to_endos[species].append(endo)
    #                 else:
    #                     self.organisms_to_endos[species] = [endo]
    #                     if self.orgChoice.findText(species) == -1:
    #                         self.orgChoice.addItem(species)

    #         if found == False:
    #             return False

    #         self.endoChoice.clear()
    #         self.endoChoice.addItems(self.organisms_to_endos[str(self.orgChoice.currentText())])
    #         self.orgChoice.currentIndexChanged.connect(self.changeEndos)
    #     except Exception as e:
    #         show_error("Error in getData() in main.", e)

    # def changeEndos(self):
    #     try:
    #         if self.orgChoice.currentText() != "Custom Input Sequences":
    #             self.Step2.setEnabled(True)
    #             self.endoChoice.setEnabled(True)
    #             self.radioButton_Gene.show()
    #             self.radioButton_Position.show()
    #             self.endoChoice.clear()
    #             self.endoChoice.addItems(self.organisms_to_endos[str(self.orgChoice.currentText())])
    #         else:
    #             self.Step2.setEnabled(False)
    #             self.endoChoice.clear()
    #             self.endoChoice.setEnabled(False)
    #             self.radioButton_Gene.hide()
    #             self.radioButton_Position.hide()
    #     except Exception as e:
    #         show_error("Error in changeEndos() in main", e)

    # def change_directory(self):
    #     try:
    #         mydir = QtWidgets.QFileDialog.getExistingDirectory(
    #             None, "Open a folder...", self.dbpath, QtWidgets.QFileDialog.Option.ShowDirsOnly)

    #         if not os.path.isdir(mydir):
    #             show_message(
    #                 fontSize=12,
    #                 icon=QtWidgets.QMessageBox.Icon.Critical,
    #                 title="Not a directory",
    #                 message="The directory you selected does not exist."
    #             )
    #             return

    #         if not any(file.endswith(".cspr") for file in os.listdir(mydir)):
    #             show_message(
    #                 fontSize=12,
    #                 icon=QtWidgets.QMessageBox.Icon.Critical,
    #                 title="Directory is invalid!",
    #                 message="You must select a directory with CSPR Files!"
    #             )
    #             return

    #         os.chdir(mydir)
    #         mydir = mydir.replace("/", "\\") if platform.system() == "Windows" else mydir
    #         GlobalSettings.CSPR_DB = mydir

    #         GlobalSettings.MTWin.directory = mydir
    #         GlobalSettings.MTWin.get_data()
    #         GlobalSettings.pop_Analysis.get_data()
    #         self.getData()
    #         self.fill_annotation_dropdown()
    #     except Exception as e:
    #         show_error("Error in change_directory() in main.", e)

    # def changeto_multitargeting(self):
    #     try:
    #         os.chdir(os.getcwd())
    #         if GlobalSettings.MTWin.first_show == True:
    #             GlobalSettings.MTWin.show()
    #             GlobalSettings.MTWin.first_show = False
    #         else:
    #             GlobalSettings.MTWin.show()
    #             GlobalSettings.mainWindow.hide()

    #     except Exception as e:
    #         show_error("Error in changeto_multitargeting() in main.", e)

    # #change to population analysis window
    # def changeto_population_Analysis(self):
    #     try:
    #         GlobalSettings.pop_Analysis.launch()
    #         if GlobalSettings.pop_Analysis.first_show == True:
    #             center_ui(GlobalSettings.pop_Analysis)
    #             GlobalSettings.pop_Analysis.first_show = False
    #         GlobalSettings.pop_Analysis.show()
    #         GlobalSettings.mainWindow.hide()
    #     except Exception as e:
    #         show_error("Error in changeto_population_Analysis() in main.", e)

    # def annotation_information(self):
    #     try:
    #         show_message(
    #             fontSize=12,
    #             icon=QtWidgets.QMessageBox.Icon.Critical,
    #             title="Annotation Information",
    #             message="Annotation files are used for searching for spacers on a gene/locus basis and can be selected here using either " \
    #                     "NCBI databases or a local file."
    #         )
    #     except Exception as e:
    #         show_error("Error in annotation_information() in main.", e)
   
    # @QtCore.pyqtSlot()
    # def view_results(self):
    #     try:
    #         #center results window on current screen
    #         if self.Results.first_show == True:
    #             self.Results.first_show = False
    #             self.Results.centerUI()

    #         self.Results.show()
    #         self.hide()
    #     except Exception as e:
    #         show_error("Error in view_results() in main", e)

    # def closeFunction(self):
    #     try:
    #         # Attempt to close the NCBI window if it exists
    #         try:
    #             self.ncbi.close()
    #         except AttributeError:
    #             print("No NCBI window to close.")
            
    #         self.myClosingWindow.get_files()
    #         center_ui(self.myClosingWindow)
    #         self.myClosingWindow.show()
    #     except Exception as e:
    #         show_error("Error in closeFunction() in main", e)

    # def close_app(self):
    #     try:
    #         # Attempt to close the NCBI window if it exists
    #         try:
    #             self.ncbi.close()
    #         except Exception as e:
    #             print("No NCBI window to close.")

    #         self.closeFunction()
    #         self.close()
    #     except Exception as e:
    #         show_error("Error in close_app() in main", e)

    # def load_dropdown_data(self):
    #     """Fill in organism/endo/annotation dropdown information."""
    #     try:
    #         self.getData()
    #         self.fill_annotation_dropdown()
    #         # self.logger.debug("Successfully loaded organism/endo/annotation drop down information in Main.")
    #     except Exception as e:
    #         show_error("Error in load_dropdown_data() in Main", e)

    # # Call methods for other windows if needed
    # # self.load_mt_data()
    # # self.load_pop_analysis_data()
    