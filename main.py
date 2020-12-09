import sys
import os, platform
import io
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from APIs import Kegg
from CoTargeting import CoTargeting
from closingWin import closingWindow
from Results import Results, geneViewerSettings
from NewGenome import NewGenome
from NewEndonuclease import NewEndonuclease
import genomeBrowser
import gzip
import webbrowser
import requests
import GlobalSettings
import multitargeting
from AnnotationParser import Annotation_Parser
from export_to_csv import export_csv_window
from cspr_chromosome_selection import cspr_chromosome_selection
from generateLib import genLibrary
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
import populationAnalysis
import platform
import ncbi

# =========================================================================================
# CLASS NAME: AnnotationsWindow
# Inputs: Greg: fill in this information
# Outputs: Greg: fill in this information
# =========================================================================================

class AnnotationsWindow(QtWidgets.QMainWindow):

    def __init__(self, info_path):
        super(AnnotationsWindow, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'Annotation Details.ui', self)
        self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
        self.Submit_button.clicked.connect(self.submit)
        self.Go_Back_Button.clicked.connect(self.go_Back)
        self.select_all_checkbox.stateChanged.connect(self.select_all_genes)
        self.mainWindow = ""
        self.type = ""
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window

    def submit(self):
        if self.type == "kegg":
            self.mainWindow.collect_table_data()
            self.hide()
            self.mainWindow.mwfg.moveCenter(self.mainWindow.cp)  ##Center window
            self.mainWindow.move(self.mainWindow.mwfg.topLeft())  ##Center window
            self.mainWindow.show()
        elif self.type == "nonkegg":
            self.mainWindow.collect_table_data_nonkegg()
            self.hide()
            self.mainWindow.mwfg.moveCenter(self.mainWindow.cp)  ##Center window
            self.mainWindow.move(self.mainWindow.mwfg.topLeft())  ##Center window
            self.mainWindow.show()


    def go_Back(self):
        self.tableWidget.clear()
        self.mainWindow.checkBoxes.clear()
        self.mainWindow.searches.clear()
        self.tableWidget.setColumnCount(0)
        self.mainWindow.mwfg.moveCenter(self.mainWindow.cp)  ##Center window
        self.mainWindow.move(self.mainWindow.mwfg.topLeft())  ##Center window
        self.mainWindow.show()
        self.mainWindow.progressBar.setValue(0)
        self.hide()


    # this function is very similar to the other fill_table, it just works with the other types of annotation files
    def fill_table_nonKegg(self, mainWindow):
        self.tableWidget.clearContents()
        self.mainWindow = mainWindow
        index = 0
        self.tableWidget.setColumnCount(5)
        self.mainWindow.progressBar.setValue(85)
        self.tableWidget.setHorizontalHeaderLabels("Description;Chromosome #;Type;Gene ID;Select".split(";"))
        mainWindow.checkBoxes = []
        self.type = "nonkegg"

        # below chain of loops goes through and figures out how many rows are needed
        # for searchValue in mainWindow.searches:
        #    for definition in mainWindow.searches[searchValue]:
        #        for gene in mainWindow.searches[searchValue][definition]:
        #            index += 1
        # self.tableWidget.setRowCount(index)

        index = 0
        for searchValue in mainWindow.searches:
            for definition in mainWindow.searches[searchValue]:
                for gene in mainWindow.searches[searchValue][definition]:
                    if (gene[2] == 'gene' or gene[2] == 'tRNA'):
                        self.tableWidget.setRowCount(index + 1)
                        # set the checkbox
                        ckbox = QtWidgets.QCheckBox()
                        self.tableWidget.setCellWidget(index, 4, ckbox)

                        # set the description part of the window as well as set the correct data for the checkbox
                        if definition != gene[0]:
                            defin_obj = QtWidgets.QTableWidgetItem(definition)
                            self.tableWidget.setItem(index, 0, defin_obj)
                            mainWindow.checkBoxes.append([definition])
                        else:
                            checkValue = searchValue.upper()
                            defin_obj = QtWidgets.QTableWidgetItem(checkValue)
                            self.tableWidget.setItem(index, 0, defin_obj)
                            mainWindow.checkBoxes.append([checkValue])
                        mainWindow.checkBoxes[len(mainWindow.checkBoxes) - 1].append(gene)
                        mainWindow.checkBoxes[len(mainWindow.checkBoxes) - 1].append(ckbox)

                        # set the type in the window
                        type_obj = QtWidgets.QTableWidgetItem(gene[2])
                        self.tableWidget.setItem(index, 2, type_obj)

                        # set the gene id in the window
                        gene_id_obj = QtWidgets.QTableWidgetItem(gene[0])
                        self.tableWidget.setItem(index, 3, gene_id_obj)

                        chrom_number = QtWidgets.QTableWidgetItem(str(gene[1]))
                        self.tableWidget.setItem(index, 1, chrom_number)

                        index += 1
                    if index >= 1000:
                        break
                if index >= 1000:
                    break
            if index >= 1000:
                break





        index = 0
        self.tableWidget.resizeColumnsToContents()

        # if show all is checked, show the window so the user can select the genes they want
        if mainWindow.Show_All_Results.isChecked():
            mainWindow.hide()
            self.mwfg.moveCenter(self.cp)  ##Center window
            self.move(self.mwfg.topLeft())  ##Center window
            self.show()
        else:  # show all not checked
            if (len(mainWindow.checkBoxes) > 15):  # check the size, throw an error if it is too large
                error = QtWidgets.QMessageBox.question(self, "Large File Found",
                                                       "This annotation file and search parameter yieled many matches and could cause a slow down.\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                if (error == QtWidgets.QMessageBox.No):
                    return -2
            for obj in mainWindow.checkBoxes:  # check every match
                obj[2].setChecked(True)
            self.mainWindow.collect_table_data_nonkegg()
        return 0


    # this is the connection for the select all checkbox
    # selects/deselects all the genes in the table
    def select_all_genes(self):
        # check to see if we're selecting all of them or not
        if self.select_all_checkbox.isChecked():
            select_all = True
        else:
            select_all = False

        # go through and do the selection
        if self.type == 'kegg':
            for i in range(self.tableWidget.rowCount()):
                self.tableWidget.cellWidget(i, 2).setChecked(select_all)
        elif self.type == 'nonkegg':
            for i in range(self.tableWidget.rowCount()):
                self.tableWidget.cellWidget(i, 4).setChecked(select_all)


    def fill_Table(self, mainWindow):
        self.tableWidget.clearContents()
        self.mainWindow = mainWindow
        index = 0
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setColumnCount(3)
        self.mainWindow.progressBar.setValue(25)
        self.tableWidget.setHorizontalHeaderLabels("Description;Gene ID;Select".split(";"))
        self.type = "kegg"

        mainWindow.checkBoxes = []

        index = 1000
        self.tableWidget.setRowCount(index)
        if index == 0:
            return -1
        index = 0


        for sValues in mainWindow.searches:
            for definition in mainWindow.searches[sValues]:
                defin_obj = QtWidgets.QTableWidgetItem(definition)
                self.tableWidget.setItem(index, 0, defin_obj)
                for gene in mainWindow.searches[sValues][definition]:
                    ckbox = QtWidgets.QCheckBox()
                    mainWindow.checkBoxes.append([definition + " " + gene])
                    mainWindow.checkBoxes[len(mainWindow.checkBoxes) - 1].append(ckbox)
                    gene_obj = QtWidgets.QTableWidgetItem(gene)
                    self.tableWidget.setItem(index, 1, gene_obj)
                    self.tableWidget.setCellWidget(index, 2, ckbox)
                    index = index + 1
                    if index >= 1000:
                        break
                if index >= 1000:
                    break
            if index >= 1000:
                break
        self.tableWidget.resizeColumnsToContents()
        self.mainWindow.progressBar.setValue(50)
        if mainWindow.Show_All_Results.isChecked():
            mainWindow.hide()
            self.mwfg.moveCenter(self.cp)  ##Center window
            self.move(self.mwfg.topLeft())  ##Center window
            self.show()
        else:  # Show all not checked
            if (len(mainWindow.checkBoxes) > 15):  # check the size, throw an error if it is too large
                error = QtWidgets.QMessageBox.question(self, "Large File Found",
                                                       "This annotation file and search parameter yieled many matches and could cause a slow down.\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                if (error == QtWidgets.QMessageBox.No):
                    return -2
            self.mainWindow.progressBar.setValue(65)
            for obj in mainWindow.checkBoxes:  # check every match
                obj[1].setChecked(True)
            self.mainWindow.collect_table_data()  # collect the data
        return 0


    # this function calls the closingWindow class.
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()


# =========================================================================================
# CLASS NAME: CMainWindow
# Inputs: Takes in the path information from the startup window and also all input parameters
# that define the search for targets e.g. endonuclease, organism genome, gene target etc.
# Outputs: The results of the target search process by generating a new Results window
# =========================================================================================


class CMainWindow(QtWidgets.QMainWindow):

    def __init__(self, info_path):
        super(CMainWindow, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'CASPER_main.ui', self)
        self.dbpath = ""
        self.info_path = info_path
        self.data = {}  # each org genome name and the endonucleases along with it
        self.TNumbers = {}  # the T numbers from a kegg search
        self.shortHand = {}  # each org's short name IE bacillus subtillis is bsu
        self.orgcodes = {}  # Stores the Kegg organism code by the format {full name : organism code}
        self.gene_list = {}  # list of genes (no ides what they pertain to
        self.searches = {}
        self.checkBoxes = []
        self.add_orgo = []
        self.checked_info = {}
        self.check_ntseq_info = {}  # the ntsequences that go along with the checked_info
        self.annotation_parser = Annotation_Parser()
        self.link_list = list()  # the list of the downloadable links from the NCBI search
        self.organismDict = dict()  # the dictionary for the links to download. Key is the description of the organism, value is the ID that can be found in link_list
        self.organismData = list()
        #self.ncbi = ncbi.NCBI_search_tool()




        # --- Style Modifications --- #
        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        font: 11pt "Sans Serif";
                        font: bold;
                        margin-top: 10px;}"""

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2").replace("rgb(111,181,110)", "rgb(77,158,89)"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1", "Step3").replace("rgb(111,181,110)", "rgb(53,121,93)"))

        # --- Button Modifications --- #


        #self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir.encode()))
        self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + 'cas9image.png'))
        self.pushButton_FindTargets.clicked.connect(self.gather_settings)
        self.pushButton_ViewTargets.clicked.connect(self.view_results)
        self.pushButton_ViewTargets.setEnabled(False)
        self.radioButton_Gene.clicked.connect(self.toggle_annotation)
        self.radioButton_Position.clicked.connect(self.toggle_annotation)

        self.seq_label.hide()

        self.actionUpload_New_Genome.triggered.connect(self.launch_newGenome)
        self.actionUpload_New_Endonuclease.triggered.connect(self.launch_newEndonuclease)
        self.actionOpen_Genome_Browser.triggered.connect(self.launch_newGenomeBrowser)
        self.Add_Orgo_Button.clicked.connect(self.add_Orgo)
        self.Remove_Organism_Button.clicked.connect(self.remove_Orgo)
        self.endoChoice.currentIndexChanged.connect(self.endo_Changed)
        self.GenerateLibrary.clicked.connect(self.prep_genlib)
        self.actionExit.triggered.connect(self.close_app)
        self.visit_repo.triggered.connect(self.visit_repo_func)
        self.refresh_button.clicked.connect(self.fill_annotation_dropdown)
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()
        self.Annotation_Window = AnnotationsWindow(info_path)

        # Hide Added orgo boxes
        self.Added_Org_Combo.hide()
        self.Remove_Organism_Button.hide()
        self.Added_Org_Label.hide()
        # --- Menubar commands --- #
        self.actionChange_Directory.triggered.connect(self.change_directory)
        self.actionMultitargeting.triggered.connect(self.changeto_multitargeting)
        self.actionPopulation_Analysis.triggered.connect(self.changeto_population_Analysis)
        self.actionNCBI.triggered.connect(self.open_ncbi_web_page)
        self.actionCasper2.triggered.connect(self.open_casper2_web_page)
        self.actionNCBI_BLAST.triggered.connect(self.open_ncbi_blast_web_page)

        self.Question_Button_add_org.clicked.connect(self.add_org_popup)


        # --- Setup for Gene Entry Field --- #
        self.geneEntryField.setPlainText("Example Inputs: \n"
                                         "Gene (LocusID): YOL086C  *for Saccharomyces Cerevisiae ADH1 gene* \n"
                                         "Position: chromosome,start,stop\n chromosome,start,stop...\n"
                                         "Sequence: *Pure sequence. CASPER will search for targets and report off"
                                         "targets based on the genome selected, if any*")

        # show functionalities on window
        self.newGenome = NewGenome(info_path)
        self.newEndonuclease = NewEndonuclease()
        self.CoTargeting = CoTargeting(info_path)
        self.Results = Results()
        self.gene_viewer_settings = geneViewerSettings()
        self.export_csv_window = export_csv_window()
        self.cspr_selector = cspr_chromosome_selection()
        self.genLib = genLibrary()
        self.myClosingWindow = closingWindow()
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.actionUpload_New_Genome.setEnabled(False)
        self.genomebrowser = genomeBrowser.genomebrowser()
        #GlobalSettings.mainWindow.ncbi = ncbi.NCBI_search_tool()
        self.launch_ncbi_button.clicked.connect(self.launch_ncbi)

        self.checkBox.setEnabled(False)

    def endo_Changed(self):
        self.add_orgo.clear()
        self.Add_Orgo_Combo.clear()
        self.Added_Org_Combo.clear()
        self.addOrgoCombo()
        self.Added_Org_Combo.hide()
        self.Added_Org_Label.hide()
        self.Remove_Organism_Button.hide()


    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def remove_Orgo(self):
        self.add_orgo.remove(self.Added_Org_Combo.currentText())
        self.Add_Orgo_Combo.addItem(self.Added_Org_Combo.currentText())
        self.Added_Org_Combo.removeItem(self.Added_Org_Combo.currentIndex())
        if len(self.add_orgo) == 0:
            self.Added_Org_Combo.hide()
            self.Added_Org_Label.hide()
            self.Remove_Organism_Button.hide()


    def add_Orgo(self):
        if self.Add_Orgo_Combo.currentText() == "Select Organism":
            QtWidgets.QMessageBox.question(self, "Must Select Organism",
                                           "You must select an organism to add.",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.add_orgo.append(self.Add_Orgo_Combo.currentText())
        self.Added_Org_Combo.addItem(self.Add_Orgo_Combo.currentText())
        self.Add_Orgo_Combo.removeItem(self.Add_Orgo_Combo.currentIndex())
        self.Added_Org_Combo.show()
        self.Added_Org_Label.show()
        self.Remove_Organism_Button.show()


    # this function prepares everything for the generate library function
    # it is very similar to the gather settings, how ever it stores the data instead of calling the Annotation Window class
    # it moves the data onto the generateLib function, and then opens that window
    def prep_genlib(self):
        # make sure the user actually inputs something
        inputstring = str(self.geneEntryField.toPlainText())
        if (inputstring.startswith("Example Inputs:") or inputstring == ""):
            QtWidgets.QMessageBox.question(self, "Error",
                                           "No gene has been entered. Please enter a gene.",
                                           QtWidgets.QMessageBox.Ok)
        else:
            # standardize the input
            inputstring = inputstring.lower()
            found_matches_bool = True

            # call the respective function
            self.progressBar.setValue(10)
            if self.radioButton_Gene.isChecked():
                ginput = inputstring.split(',')
                found_matches_bool = self.run_results("gene", ginput, openAnnoWindow=False)
            elif self.radioButton_Position.isChecked() or self.radioButton_Sequence.isChecked():
                QtWidgets.QMessageBox.question(self, "Error", "Generate Library can only work with gene (Locus ID).",
                                               QtWidgets.QMessageBox.Ok)
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
                fileName = self.shortHand[self.orgChoice.currentText()]
                fileName = GlobalSettings.CSPR_DB + '/' + fileName + '_' + self.endoChoice.currentText() + '.cspr'

                kegg_non = 'non_kegg'

                # launch generateLib
                self.progressBar.setValue(100)

                # calculate the total number of matches found
                tempSum = 0
                for item in self.searches:
                    tempSum += len(self.searches[item])

                # warn the user if the number is greater than 50
                if tempSum > 50:
                    error = QtWidgets.QMessageBox.question(self, "Many Matches Found",
                                                           "More than 50 matches have been found. Continuing could cause a slow down. .\n\n"
                                                           "Do you wish to continue?",
                                                           QtWidgets.QMessageBox.Yes |
                                                           QtWidgets.QMessageBox.No,
                                                           QtWidgets.QMessageBox.No)
                    if (error == QtWidgets.QMessageBox.No):
                        self.searches.clear()
                        self.progressBar.setValue(0)
                        return -2

                self.genLib.launch(self.searches, fileName, kegg_non)
            else:
                self.progressBar.setValue(0)


    # Function for collecting the settings from the input field and transferring them to run_results
    def gather_settings(self):
        inputstring = str(self.geneEntryField.toPlainText())

        # Error check: make sure the user actually inputs something
        if (inputstring.startswith("Example Inputs:") or inputstring == ""):
            QtWidgets.QMessageBox.question(self, "Error",
                                           "No gene has been entered. Please enter a gene.",
                                           QtWidgets.QMessageBox.Ok)
        else:
            # standardize the input
            inputstring = inputstring.lower()

            self.progressBar.setValue(10)
            if self.radioButton_Gene.isChecked():
                ginput = inputstring.split(',')
                self.run_results("gene", ginput)
            elif self.radioButton_Position.isChecked():
                self.run_results("position", inputstring)
            elif self.radioButton_Sequence.isChecked():
                sinput = inputstring
                self.run_results("sequence", sinput)


    # ---- Following functions are for running the auxillary algorithms and windows ---- #
    # this function is parses the annotation file given, and then goes through and goes onto results
    # it will call other versions of collect_table_data and fill_table that work with these file types
    # this function should work with the any type of annotation file, besides kegg.
    # this assumes that the parsers all store the data the same way, which gff and feature table do
    # please make sure the gbff parser stores the data in the same way
    # so far the gff files seems to all be different. Need to think about how we want to parse it
    def run_results_own_ncbi_file(self, inputstring, fileName, openAnnoWindow=True):
        self.annotation_parser = Annotation_Parser()
        self.annotation_parser.annotationFileName = fileName
        fileType = self.annotation_parser.find_which_file_version()

        # if the parser retuns the 'wrong file type' error
        if fileType == -1:
            QtWidgets.QMessageBox.question(self, "Error:",
                                           "We cannot parse the file type given. Please make sure to choose a GBFF, GFF, or Feature Table file."
                                           , QtWidgets.QMessageBox.Ok)
            self.progressBar.setValue(0)
            return

        self.progressBar.setValue(60)

        # this bit may not be needed here. Just a quick error check to make sure the chromosome numbers match
        full_org = str(self.orgChoice.currentText())
        cspr_file = (GlobalSettings.CSPR_DB + "/" + self.shortHand[full_org] + "_" + str(
            self.endoChoice.currentText()) + ".cspr")
        own_cspr_parser = CSPRparser(cspr_file)
        own_cspr_parser.read_first_lines()

        if len(own_cspr_parser.karystatsList) != self.annotation_parser.max_chrom:
            QtWidgets.QMessageBox.question(self, "Warning:",
                                           "The number of chromosomes do not match. This could cause errors."
                                           , QtWidgets.QMessageBox.Ok)

        ##may need to open a message window here
        ###################

        # now go through and search for the actual locus tag, in the case the user input that
        searchValues = self.separate_line(inputstring[0])
        self.searches.clear()
        for search in searchValues:
            search = self.removeWhiteSpace(search)
            if len(search) == 0:
                continue

            # set the searche's dict of search equal to a new dictionary
            self.searches[search] = {}

            # upper case it, because the files seem to be all uppercase for this part
            checkNormalDict = search.upper()
            if checkNormalDict in self.annotation_parser.reg_dict:  # if it is in the normal dictionary
                for item in self.annotation_parser.reg_dict[checkNormalDict]:  # for each list item in that position
                    if item[0] not in self.searches[search]:  # if its not in the search's position yet
                        self.searches[search][item[0]] = self.annotation_parser.reg_dict[checkNormalDict]
                    elif item not in self.searches[search][
                        item[0]]:  # assume it is in the searches position, but do not store duplicates
                        self.searches[search][item[0]].append(self.annotation_parser.reg_dict[checkNormalDict])
        if len(self.searches[searchValues[0]]) >= 1:  # if the previous search yielded results, do not continue
            if openAnnoWindow:
                self.Annotation_Window.fill_table_nonKegg(self)
                return
            else:
                return True

        self.progressBar.setValue(75)
        # reset, and search the parallel dictionary now
        self.searches = {}
        for search in searchValues:
            search = self.removeWhiteSpace(search)
            if len(search) == 0:
                continue

            self.searches[search] = {}
            for item in self.annotation_parser.para_dict:
                checkingItem = item.lower()  # lowercase now, to match the user's input

                if search in checkingItem:  # if what they are searching for is somewhere in that key
                    # loop through each list item now
                    for i in range(len(self.annotation_parser.para_dict[item])):
                        # now loop through the regular dict's position and store the ones we want
                        if self.annotation_parser.para_dict[item][i] != '':
                            for match in self.annotation_parser.reg_dict[self.annotation_parser.para_dict[item][i]]:
                                if item not in self.searches[search]:
                                    self.searches[search][item] = [match]
                                elif item not in self.searches[search][item]:
                                    self.searches[search][item].append(match)
        # if the search returns nothing, throw an error
        if len(self.searches[searchValues[0]]) <= 0:
            QtWidgets.QMessageBox.question(self, "No Matches Found",
                                           "No matches found with that search, please try again.",
                                           QtWidgets.QMessageBox.Ok)
            self.progressBar.setValue(0)
            if openAnnoWindow:
                return
            else:
                return False

        # if we get to this point, that means that the search yieleded results, so fill the table
        self.progressBar.setValue(80)
        # check whether this function call is for Annotation Window, or for generate Lib
        if openAnnoWindow:
            self.Annotation_Window.fill_table_nonKegg(self)
        else:
            return True


    def run_results(self, inputtype, inputstring, openAnnoWindow=True):
        if(str(self.annotation_files.currentText()).find('.gbff') == -1):
            QtWidgets.QMessageBox.information(self, "Genomebrowser Error", "Filetype must be GBFF.",
                                              QtWidgets.QMessageBox.Ok)
            self.progressBar.setValue(0)
            return


        progvalue = 15
        self.searches = {}
        self.gene_list = {}
        self.progressBar.setValue(progvalue)

        # set the gene viewer stuff to false so that the user can't open the gene viewer unless they have a
        # FNA/GBFF file
        # Note: if they are using kegg, the enabled's are set to True because it uses the kegg database by default
        self.Results.displayGeneViewer.setEnabled(False)
        self.Results.lineEditStart.setEnabled(False)
        self.Results.lineEditEnd.setEnabled(False)
        self.gene_viewer_settings.kegg_radio_button.setEnabled(False)
        self.Results.change_start_end_button.setEnabled(False)
        self.Results.displayGeneViewer.setChecked(0)

        if inputtype == "gene":
            # make sure an annotation file has been selected
            if self.annotation_files.currentText() == "":
                error = QtWidgets.QMessageBox.question(self, "No Annotation", "Please select an annotation from either KEGG, NCBI, or provide you own annotation file", QtWidgets.QMessageBox.Ok)
                self.progressBar.setValue(0)
                return
            # this now just goes onto the other version of run_results
            myBool = self.run_results_own_ncbi_file(inputstring, self.annotation_files.currentText(), openAnnoWindow=openAnnoWindow)
            if not openAnnoWindow:
                return myBool
            else:
                self.progressBar.setValue(0)
                return

        # position code below
        if inputtype == "position":
            inputstring = inputstring.replace(' ', '')
            searchInput = inputstring.split('\n')
            full_org = str(self.orgChoice.currentText())
            self.checked_info.clear()
            self.check_ntseq_info.clear()

            for item in searchInput:
                searchIndicies = item.split(',')
                # make sure the right amount of arguments were passed
                if len(searchIndicies) != 3:
                    QtWidgets.QMessageBox.question(self, "Position Error: Invalid Input",
                                                   "There are 3 arguments required for this function: chromosome, start position, and end position.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return

                # make sure user inputs digits
                if not searchIndicies[0].isdigit() or not searchIndicies[1].isdigit() or not searchIndicies[2].isdigit():
                    QtWidgets.QMessageBox.question(self, "Position Error: Invalid Input",
                                                   "The positions given must be integers. Please try again.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return
                # make sure start is less than end
                elif int(searchIndicies[1]) >= int(searchIndicies[2]):
                    QtWidgets.QMessageBox.question(self, "Position Error: Start Must Be Less Than End",
                                                   "The start index must be less than the end index.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return
                # append the data into the checked_info
                tempString = 'chrom: ' + str(searchIndicies[0]) + ' start: ' + str(searchIndicies[1]) + ' end: ' + str(searchIndicies[2])
                self.checked_info[tempString] = (int(searchIndicies[0]), int(searchIndicies[1]), int(searchIndicies[2]))

            self.progressBar.setValue(50)
            self.Results.transfer_data(self.shortHand[full_org], [str(self.endoChoice.currentText())], os.getcwd(), self.checked_info, self.check_ntseq_info, "")
            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)

        # sequence code below
        if inputtype == "sequence":
            checkString = 'AGTCN'
            self.checked_info.clear()
            self.progressBar.setValue(10)
            inputstring = inputstring.upper()

            # check to make sure that the use gave a long enough sequence
            if len(inputstring) < 100:
                QtWidgets.QMessageBox.question(self, "Error",
                                               "The sequence given is too small. At least 100 characters are required.",
                                               QtWidgets.QMessageBox.Ok)
                self.progressBar.setValue(0)
                return

            # give a warning if the length of the sequence is long
            if len(inputstring) > 30000:
                error = QtWidgets.QMessageBox.question(self, "Large Sequence Detected",
                                                       "The sequence given is a large one and could slow down the process.\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                if (error == QtWidgets.QMessageBox.No):
                    self.progressBar.setValue(0)
                    return

            # make sure all the chars are one of A, G, T, C, or N
            for letter in inputstring:
                # skip the end line character
                if letter == '\n':
                    continue
                if letter not in checkString:
                    QtWidgets.QMessageBox.question(self, "Sequence Error",
                                                   "The sequence must consist of A, G, T, C, or N. No other characters are allowed.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return
            self.progressBar.setValue(30)

            # build the CSPR file, and go into results
            fna_file_path = GlobalSettings.CSPR_DB + '/temp.fna'
            self.checked_info['Sequence Finder'] = (1, 0, len(inputstring))
            self.check_ntseq_info['Sequence Finder'] = inputstring.replace('\n', '')
            outFile = open(fna_file_path, 'w')
            outFile.write('>temp org here\n')
            outFile.write(inputstring)
            outFile.write('\n\n')
            outFile.close()
            self.progressBar.setValue(55)
            self.build_cspr_for_seq(fna_file_path)


    # this function actually builds the cspr file
    # I am not pulling out the data for the progress bar here, but it could be added easily. Figured it would be best without it since the progress should run quickly
    def build_cspr_for_seq(self, fna_file_path):
        seq_search_process = QtCore.QProcess()

        # this function decides what to do when the process is finished.
        # what it does:
        #       Deletes the temp FNA file, and calls transfer_data in results.
        def finished():
            # get the file name and process all of the targets
            cspr_file_name = GlobalSettings.CSPR_DB + '/' + code + '_' + myEndoChoice + '.cspr'
            os.remove(fna_file_path)
            seq_search_process.kill()
            self.Results.transfer_data('tempCode', [str(self.endoChoice.currentText())], os.getcwd(),
                                       self.checked_info, self.check_ntseq_info, "")

            self.Results.displayGeneViewer.setEnabled(True)
            self.Results.gene_viewer_settings_button.setEnabled(False)
            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)
            self.Results.seq_finder_cspr_file = cspr_file_name

        # --------------------getting all of the arguments-------------------
        myEndoChoice = str(self.endoChoice.currentText())
        my_seq = SeqTranslate()
        pam = my_seq.endo_info[myEndoChoice][0].split(',')[0]
        code = 'tempCode'
        if int(my_seq.endo_info[myEndoChoice][3]) == 3:
            pamdir = False
        else:
            pamdir = True

        output_location = GlobalSettings.CSPR_DB
        path_to_info = GlobalSettings.appdir + 'CASPERinfo'
        orgName = 'temp org'
        gRNA_length = my_seq.endo_info[myEndoChoice][2]
        seed_length = my_seq.endo_info[myEndoChoice][1]
        secondCode = 'test second code'
        # ----------------------------------------------------------------------

        # ------------compile the executable line------------------------------------
        program = '"' + GlobalSettings.appdir + '/Casper_Seq_Finder_Windows' + '" '
        args = '"' + myEndoChoice + '" '
        args = args + '"' + pam + '" '
        args = args + '"' + code + '" '
        args = args + str(pamdir) + ' '
        args = args + '"' + output_location + '" '
        args = args + '"' + path_to_info + '" '
        args = args + '"' + fna_file_path + '" '
        args = args + '"' + orgName + '" '
        args = args + gRNA_length + ' '
        args = args + seed_length + ' '
        args = args + '"' + secondCode + '"'
        program = program + args
        # ----------------------------------------------------------------------------

        # start the process
        seq_search_process.start(program)
        seq_search_process.finished.connect(finished)


    def launch_newGenome(self):
        self.hide()
        self.newGenome.mwfg.moveCenter(self.newGenome.cp)  ##Center window
        self.newGenome.move(self.mwfg.topLeft())  ##Center window
        self.newGenome.show()


    def launch_newEndonuclease(self):
        self.newEndonuclease.show()


    def launch_newGenomeBrowser(self):
        print("creating graph")
        self.genomebrowser.createGraph(self)


    def launch_ncbi(self):
        self.ncbi = ncbi.NCBI_search_tool()
        self.ncbi.show()

    # this function does the same stuff that the other collect_table_data does, but works with the other types of files
    def collect_table_data_nonkegg(self):
        # start out the same as the other collect_table_data
        self.checked_info.clear()
        self.check_ntseq_info.clear()
        full_org = str(self.orgChoice.currentText())
        holder = ()

        for item in self.checkBoxes:
            if item[2].isChecked():
                # if they searched base on Locus Tag
                if item[0] in self.annotation_parser.reg_dict:
                    # go through the dictionary, and if they match, store the item in holder
                    for match in self.annotation_parser.reg_dict[item[0]]:
                        if item[1] == match:
                            holder = (match[1], match[3], match[4])
                            self.checked_info[item[0]] = holder
                else:
                    # now we need to go through the para_dict
                    for i in range(len(self.annotation_parser.para_dict[item[0]])):
                        # now go through the matches in the normal dict's data
                        for match in self.annotation_parser.reg_dict[self.annotation_parser.para_dict[item[0]][i]]:
                            # if they match, store it in holder
                            if item[1] == match:
                                holder = (match[1], match[3], match[4])
                                self.checked_info[item[0]] = holder

        # now call transfer data
        self.progressBar.setValue(95)
        self.Results.transfer_data(self.shortHand[full_org], [str(self.endoChoice.currentText())], os.getcwd(),
                                   self.checked_info, self.check_ntseq_info, "")
        self.progressBar.setValue(100)
        self.pushButton_ViewTargets.setEnabled(True)


    def collect_table_data(self):
        # need to change this code so that it works with other types of files
        self.checked_info.clear()
        self.check_ntseq_info.clear()

        k = Kegg()
        full_org = str(self.orgChoice.currentText())
        organism = self.Annotations_Organism.currentText().split(" ")[1]
        nameFull = ""
        holder = ()
        for item in self.checkBoxes:
            if item[1].isChecked() == True:
                nameFull = item[0].split(" ")
                name = nameFull[len(nameFull) - 1]
                gene_info = k.gene_locator(organism + ":" + name)
                # print(nameFull)

                # get kegg's ntsequence and store it
                if gene_info != -1:
                    nt_sequence = k.get_nt_sequence(organism + ":" + name)

                    # print(item[0])
                    holder = (gene_info[0], gene_info[2], gene_info[3], gene_info[4])
                    self.checked_info[item[0]] = holder
                    self.check_ntseq_info[item[0]] = nt_sequence

        if len(self.checked_info) > 0:
            self.progressBar.setValue(80)
            self.Results.transfer_data(self.shortHand[full_org], [str(self.endoChoice.currentText())], os.getcwd(),
                                       self.checked_info, self.check_ntseq_info, "")
            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)
        else:
            print("No items were found. Please search again.")

    # ------------------------------------------------------------------------------------------------------ #

    # ----- Following Code is helper functions for processing input data ----- #
    def separate_line(self, input_string):
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


    def removeWhiteSpace(self, strng):
        while True:
            if len(strng) == 0 or (strng[0] != " " and strng[0] != "\n"):
                break
            strng = strng[1:]
        while True:
            if len(strng) == 0 or (strng[len(strng) - 1] != " " and strng[0] != "\n"):
                return strng
            strng = strng[:len(strng) - 1]


    # Function to enable and disable the Annotation function if searching by position or sequence
    def toggle_annotation(self):
        if self.radioButton_Gene.isChecked():
            self.Step2.setEnabled(True)
        else:
            self.Step2.setEnabled(True)

        # check to see if the sequence button is pressed, and act accordingly -- OLD code
        # elif self.radioButton_Sequence.isChecked():
        #     self.Step2.setEnabled(False)
        #     mySeq = SeqTranslate()
        #     seq_checker = False
        #     # time to reset the endo's
        #     self.endoChoice.clear()
        #     for item in mySeq.endo_info:
        #         self.endoChoice.addItem(item)
        # else:
        #     self.Step2.setEnabled(False)
        #     seq_checker = True
        #     self.changeEndos()

        #self.orgChoice.setEnabled(seq_checker)



    def fill_annotation_dropdown(self):
        self.annotation_files.clear()
        for file in os.listdir(GlobalSettings.CSPR_DB):
            if ".gbff" in file or ".gff" in file:
                self.annotation_files.addItem(file)


    def make_dictonary(self):
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


    def organism_finder(self, long_str):
        semi = long_str.find(";")
        index = 1
        while True:
            if long_str[semi - index] == " ":
                break
            index = index + 1
        return long_str[:semi - index]


    # This method is for testing the execution of a button call to make sure the button is linked properly
    def testexe(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to quit?",
                                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            # print(self.orgChoice.currentText())
            sys.exit()
        else:
            pass


    def addOrgoCombo(self):
        self.Add_Orgo_Combo.addItem("Select Organism")
        for item in self.data:
            if (self.endoChoice.currentText() in self.data[item]) and (item != str(self.orgChoice.currentText())):
                self.Add_Orgo_Combo.addItem(item)


    # ----- CALLED IN STARTUP WINDOW ------ #
    def getData(self):
        mypath = os.getcwd()
        found = False
        self.dbpath = mypath
        onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
        orgsandendos = {}
        shortName = {}
        first = True
        for file in onlyfiles:
            if file.find('.cspr') != -1:
                if first == True:
                    first = False
                found = True
                newname = file[0:-4]
                s = newname.split('_')
                hold = gzip.open(file, 'r')
                buf = (hold.readline())
                buf = str(buf)
                buf = buf.strip("'b")
                buf = buf[:len(buf) - 4]
                species = buf[8:]
                endo = str(s[1][:len(s[1]) - 1])
                if species not in shortName:
                    shortName[species] = s[0]
                if species in orgsandendos:
                    orgsandendos[species].append(endo)
                else:
                    orgsandendos[species] = [endo]
                    if self.orgChoice.findText(species) == -1:
                        self.orgChoice.addItem(species)

        #self.orgChoice.addItem("Custom Input Sequences")
        # auto fill the kegg search bar with the first choice in orgChoice
        if found == False:
            return False
        self.data = orgsandendos
        self.shortHand = shortName
        self.endoChoice.clear()
        self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])
        self.orgChoice.currentIndexChanged.connect(self.changeEndos)


    def changeEndos(self):
        if self.orgChoice.currentText() != "Custom Input Sequences":
            self.Step2.setEnabled(True)
            self.endoChoice.setEnabled(True)
            self.seq_label.hide()
            self.radioButton_Gene.show()
            self.radioButton_Position.show()
            self.endoChoice.clear()
            file = str(self.shortHand[str(self.orgChoice.currentText())]) + '_' + str(self.data[str(self.orgChoice.currentText())][0]) + '.cspr'
            self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])
        else:
            self.Step2.setEnabled(False)
            self.endoChoice.clear()
            self.endoChoice.setEnabled(False)
            self.radioButton_Gene.hide()
            self.radioButton_Position.hide()
            self.seq_label.show()


    def change_directory(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                           self.dbpath, QtWidgets.QFileDialog.ShowDirsOnly)

        if (os.path.isdir(mydir)):
            os.chdir(mydir)
        self.getData()


    # Tanner - added this function to allow the Tools->Multitargeting button to work
    # Function launches the multitargeting window and closing the current one
    def changeto_multitargeting(self):
        os.chdir(os.getcwd())
        GlobalSettings.MTWin.mwfg.moveCenter(GlobalSettings.MTWin.cp)  ##Center window
        GlobalSettings.MTWin.move(GlobalSettings.MTWin.mwfg.topLeft())  ##Center window
        GlobalSettings.MTWin.show()
        GlobalSettings.mainWindow.hide()


    def changeto_population_Analysis(self):
        GlobalSettings.pop_Analysis.launch(GlobalSettings.CSPR_DB)
        GlobalSettings.pop_Analysis.mwfg.moveCenter(GlobalSettings.pop_Analysis.cp)  ##Center window
        GlobalSettings.pop_Analysis.move(GlobalSettings.pop_Analysis.mwfg.topLeft())  ##Center window
        GlobalSettings.pop_Analysis.show()
        GlobalSettings.mainWindow.hide()


    def add_org_popup(self):
        info = "This functionality will allow users to use different organisms for off-target analysis in a future " \
               "version of the software. If you need to run analysis on multiple organisms, please use the Population " \
               "Analysis feature."
        QtWidgets.QMessageBox.information(self, "Add Organism Information", info, QtWidgets.QMessageBox.Ok)


    def annotation_information(self):
        info = "Annotation files for searching for targets on a gene/locus basis can be selected here using either KEGG " \
               "or NCBI databases, or uploading your own file. Note that KEGG searches are best done with exact matches " \
               "(e.g include strain designation), whereas NCBI will often return multiple assemblies of the same species. " \
               "If you have trouble deciding on the search returned annotation, go to the website, download the annotation " \
               "file to your local computer and choose ‘Choose Annotation File’"
        QtWidgets.QMessageBox.information(self, "Annotation Information", info, QtWidgets.QMessageBox.Ok)


    def open_ncbi_blast_web_page(self):
        webbrowser.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi', new=2)


    def open_ncbi_web_page(self):
        webbrowser.open('https://www.ncbi.nlm.nih.gov/', new=2)


    def open_casper2_web_page(self):
        webbrowser.open('http://casper2.org/', new=2)


    def visit_repo_func(self):
        webbrowser.open('https://github.com/TrinhLab/CASPERapp')


    @QtCore.pyqtSlot()
    def view_results(self):
        self.hide()

        # set Results endo combo box
        self.Results.endonucleaseBox.clear()

        # set the results window endoChoice box menu
        # set the mainWindow's endoChoice first, and then loop through and set the rest of them
        self.Results.endonucleaseBox.addItem(self.endoChoice.currentText())
        for item in self.data[str(self.orgChoice.currentText())]:
            if item != self.Results.endonucleaseBox.currentText():
                self.Results.endonucleaseBox.addItem(item)

        self.Results.mwfg.moveCenter(self.Results.cp)  ##Center window
        self.Results.move(self.Results.mwfg.topLeft())  ##Center window
        self.Results.show()


    # this function calls the closingWindow class.
    def closeEvent(self, event):
        self.closeFunction()
        event.accept()


    # this if the function that is called when the user closes the program entirely.
    # so far I only know of 4 spots that can do this
    #       1. mainWindow
    #       2. annotationsWindow
    #       3. Results
    #       4. Multitargetting
    def closeFunction(self):
        try:
            self.ncbi.close()
        except:
            print("no ncbi window to close")
        self.myClosingWindow.get_files()
        self.myClosingWindow.show()


    def close_app(self):
        try:
            self.ncbi.close()
        except:
            print("no ncbi window to close")

        self.closeFunction()
        self.close()


# ----------------------------------------------------------------------------------------------------- #
# =========================================================================================
# CLASS NAME: StartupWindow
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class StartupWindow(QtWidgets.QDialog):
    def __init__(self):
        super(StartupWindow, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'startupCASPER.ui', self)
        self.setWindowModality(2)  # sets the modality of the window to Application Modal
        self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
        pixmap = QtGui.QPixmap(GlobalSettings.appdir + 'CASPER-logo.jpg')
        self.labelforart.setPixmap(pixmap)
        self.pushButton_2.setDefault(True)
        # Check to see the operating system you are on and change this in Global Settings:
        GlobalSettings.OPERATING_SYSTEM_ID = platform.system()
        self.info_path = os.getcwd()
        self.gdirectory = self.check_dir()
        GlobalSettings.CSPR_DB = self.gdirectory
        self.lineEdit.setText(self.gdirectory)

        self.pushButton_3.clicked.connect(self.changeDir)
        self.pushButton_2.clicked.connect(self.show_window)
        #self.pushButton.clicked.connect(self.errormsgmulti)

        self.show()
        self.pushButton.setEnabled(False)


    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def changeDir(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                           self.gdirectory, QtWidgets.QFileDialog.ShowDirsOnly)


        if (os.path.isdir(mydir) == False):
            QtWidgets.QMessageBox.critical(self, "Invalid Directory!", "Invalid Directory!", QtWidgets.QMessageBox.Ok)
            return

        found = False
        for file in os.listdir(mydir):
            if(file.find(".cspr") != -1):
                found = True
                break

        if(found == False):
            QtWidgets.QMessageBox.critical(self, "Directory is invalid!", "You must select a directory with CSPR Files!",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.lineEdit.setText(mydir)
        self.gdirectory = mydir
        GlobalSettings.CSPR_DB = mydir


    def errormsgmulti(self):
        self.gdirectory = str(self.lineEdit.text())
        # print(self.gdirectory)
        if "Please select a directory that contains .cspr files" in self.gdirectory:
            QtWidgets.QMessageBox.question(self, "Must select directory", "You must select your directory.",
                                           QtWidgets.QMessageBox.Ok)

        elif (os.path.isdir(self.gdirectory)):
            os.chdir(self.gdirectory)
            # change dir, still load main window, still load MT data, and then open main window and newGenome window
            GlobalSettings.filedir = self.gdirectory
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            self.re_write_dir()
            GlobalSettings.mainWindow.launch_newGenome()
            self.close()
        else:
            QtWidgets.QMessageBox.question(self, "Not a directory", "The directory you selected does not exist.",
                                           QtWidgets.QMessageBox.Ok)


    def check_dir(self):
        cspr_info = open(GlobalSettings.appdir + "CASPERinfo", 'r+')
        cspr_info = cspr_info.read()
        lines = cspr_info.split('\n')
        line = ""
        for item in lines:
            if 'DIRECTORY:' in item:
                line = item
                break
        if len(line) < 11:
            return os.path.expanduser("~\Documents").replace('\\', '/')
        else:
            return line[10:]


    def re_write_dir(self):
        cspr_info = open(GlobalSettings.appdir + "CASPERinfo", 'r+')
        cspr_info_text = cspr_info.read()
        cspr_info_text = cspr_info_text.split('\n')
        full_doc = ""
        for item in cspr_info_text:
            if 'DIRECTORY:' in item:
                line = item
                break
        line_final = "DIRECTORY:" + self.gdirectory
        for item in cspr_info_text:
            if item == line:
                full_doc = full_doc + "\n" + line_final
            else:
                full_doc = full_doc + "\n" + item
        full_doc = full_doc[1:]
        cspr_info.close()
        cspr_info = open(GlobalSettings.appdir + "CASPERinfo", 'r+')
        cspr_info.write(full_doc)

        cspr_info.close()


    def show_window(self):
        self.gdirectory = str(self.lineEdit.text())
        # print(self.gdirectory)
        if "Please select a directory that contains .cspr files" in self.gdirectory:
            QtWidgets.QMessageBox.question(self, "Must select directory", "You must select your directory.",
                                           QtWidgets.QMessageBox.Ok)
        elif (os.path.isdir(self.gdirectory)):

            os.chdir(self.gdirectory)
            found = GlobalSettings.mainWindow.getData()
            GlobalSettings.mainWindow.fill_annotation_dropdown()
            if found == False:
                QtWidgets.QMessageBox.question(self, "No .cspr files",
                                               "Please select a directory that contains cspr files.",
                                               QtWidgets.QMessageBox.Ok)
                return
            GlobalSettings.filedir = self.gdirectory
            self.re_write_dir()
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            GlobalSettings.mainWindow.show()

            # Tanner - still setup data for MT
            GlobalSettings.MTWin.launch(self.gdirectory)

            self.close()
        else:
            QtWidgets.QMessageBox.question(self, "Not a directory", "The directory you selected does not exist.",
                                           QtWidgets.QMessageBox.Ok)


def main():
    if hasattr(sys, 'frozen'):
        GlobalSettings.appdir = sys.executable
        if platform.system() == 'Windows':
            GlobalSettings.appdir = GlobalSettings.appdir[:GlobalSettings.appdir.rfind("\\") + 1]
        else:
            GlobalSettings.appdir = GlobalSettings.appdir[:GlobalSettings.appdir.rfind("/") + 1]
    else:
        GlobalSettings.appdir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            GlobalSettings.appdir += '\\'
        else:
            GlobalSettings.appdir += '/'

    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = Qt.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")
    startup = StartupWindow()
    GlobalSettings.mainWindow = CMainWindow(os.getcwd())
    GlobalSettings.MTWin = multitargeting.Multitargeting()
    GlobalSettings.pop_Analysis = populationAnalysis.Pop_Analysis()
    sys.exit(app.exec_())

if __name__ == '__main__':

    try:
        main()
    except Exception as e:
        with open("debug.txt", "w") as f:
            f.write(e)