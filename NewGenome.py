import sys, os
import subprocess as sub
from threading import Thread
from queue import Queue, Empty
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
from bioservices import KEGG
from NCBI_API import Assembly, GBFF_Parse
import GlobalSettings


def iter_except(function, exception):
    """Works like builtin 2-argument `iter()`, but stops on `exception`."""
    try:
        while True:
            yield function()
    except exception:
        return


##############################################
# Class-Name: NCBI_Search_File
# What it does: this is the ncbi search window that opens up. Allows the user to search for a ncbi refseq or genbank file
##############################################
class NCBI_Search_File(QtWidgets.QDialog):
    def __init__(self):
        # Qt init stuff
        super(NCBI_Search_File, self).__init__()
        uic.loadUi("NCBI_File_Search.ui", self)
        self.setWindowTitle("Search NCBI For a File")
        self.searchProgressBar.setValue(0)
        self.setWindowIcon(Qt.QIcon("cas9image.png"))

        # selection table stuff
        self.selectionTableWidget.setColumnCount(1)  # hardcoded because there will always be 1 columns
        self.selectionTableWidget.setShowGrid(True)
        self.selectionTableWidget.setHorizontalHeaderLabels("Description;Selection".split(";"))
        self.selectionTableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.selectionTableWidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.selectionTableWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        # submission table stuff
        self.submissionTableWidget.setColumnCount(1) # hardcoded because it will always be 1
        self.submissionTableWidget.setShowGrid(True)
        self.submissionTableWidget.setHorizontalHeaderLabels("Description;Selection".split(";"))
        self.submissionTableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.submissionTableWidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.submissionTableWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        # variables below:
        self.searchType = "" #RefSeq or Genbank
        self.database_url_list = list()
        self.ncbi_searcher = Assembly()
        self.compressedFilePaths = list()
        self.decompressedFilePaths = list()

        # dictionary: key = the description of the organism
        # value = the GCA_ID that links it to the database link
        self.idList_dict = dict()

        # button modifications below:
        self.searchPushButton.clicked.connect(self.searchFunction)
        self.selectButton.clicked.connect(self.selectHighlighted)
        self.deselectButton.clicked.connect(self.deselectHighlighted)
        self.cancelPushButton.clicked.connect(self.cancelFunction)
        self.submitPushButton.clicked.connect(self.submissionFunction)


    # button functions below
    # this function searchs ncbi for the organism file, adds any matches to the selectionTableWidget
    def searchFunction(self):

        # if there is nothing set in the line edit, set the ret max to 20
        if self.ncbi_ret_max_line_edit.displayText() != "":
            # if the string is only digits (and non negative), set ret max to that, otherwise just set it to 20
            if self.ncbi_ret_max_line_edit.displayText().isdigit() and int(self.ncbi_ret_max_line_edit.displayText()) > 0:
                ret_max = int(self.ncbi_ret_max_line_edit.displayText())
            else:
                ret_max = 20
        else:
            ret_max = 20

        # clear any previous searches
        self.selectionTableWidget.clearContents()
        self.submissionTableWidget.clearContents()
        self.selectionTableWidget.setRowCount(0)
        self.submissionTableWidget.setRowCount(0)

        # check and make sure that the user enters an organism
        if self.organismLineEdit.displayText() == "":
            QtWidgets.QMessageBox.question(self, "No Organism Entered",
                                           "No organism has been entered to search for. "
                                           "Please type in an organism to search for.",
                                           QtWidgets.QMessageBox.Ok)
            return
        # make sure that the retmax value is not too large
        if ret_max > 100:
            QtWidgets.QMessageBox.question(self, "Error",
                                           "Return Max number is too high, please set it to something below 100",
                                           QtWidgets.QMessageBox.Ok)
            return

        else:
            # search the ncbi database
            self.searchProgressBar.setValue(15)
            searchOrganism = self.organismLineEdit.displayText()
            self.database_url_list, self.idList_dict = self.ncbi_searcher.getDataBaseURL(searchOrganism, self.searchType, ret_max)
            self.searchProgressBar.setValue(85)


            # check and make sure something was found
            if len(self.idList_dict) <= 0:
                QtWidgets.QMessageBox.question(self, "No matches found!", "No matches found for the search, please try "
                                                                          "again.",
                                               QtWidgets.QMessageBox.Ok)
                self.searchProgressBar.setValue(0)
                return
            # now update the table
            self.selectionTableWidget.setRowCount(len(self.idList_dict))
            loopCount = 0
            for item in self.idList_dict:
                tabWidget = QtWidgets.QTableWidgetItem(item)
                self.selectionTableWidget.setItem(loopCount, 0, tabWidget)
                loopCount += 1
            self.selectionTableWidget.resizeColumnsToContents()
            self.searchProgressBar.setValue(100)

    # this function takes the selected table widgets, and adds them to the submissionTableWidget
    def selectHighlighted(self):
        # get the selected items from selectionTableWidget
        selectedList = self.selectionTableWidget.selectedItems()

        # print an error if the user doesn't select anything
        if len(selectedList) <= 0 and self.selectionTableWidget.rowCount() > 0:
            QtWidgets.QMessageBox.question(self, "Nothing Seleted", "No items selected. Please selected items. "
                                                                      "again.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # break out if the size of the selectionTableWidget is 0
        elif self.selectionTableWidget.rowCount() <= 0:
            return

        # then set up the next table
        self.submissionTableWidget.setRowCount(len(selectedList))
        for i in range(len(selectedList)):
            tabWidget = QtWidgets.QTableWidgetItem(selectedList[i])
            self.submissionTableWidget.setItem(i, 0, tabWidget)
        self.submissionTableWidget.resizeColumnsToContents()

    # this function takes all of the selected widgets out of the submissionTableWidget
    def deselectHighlighted(self):
        # get the selected indicies
        selectedList = self.submissionTableWidget.selectedItems()

        # error checking
        if len(selectedList) <= 0 and self.submissionTableWidget.rowCount() > 0:
            QtWidgets.QMessageBox.question(self, "Nothing Seleted", "No items selected. Please selected items. "
                                                                    "again.",
                                           QtWidgets.QMessageBox.Ok)
            return

        elif self.submissionTableWidget.rowCount() <= 0:
            return

        # go through the delete the ones selected
        for i in range(self.submissionTableWidget.rowCount()):
            for j in range(len(selectedList)):
                if self.submissionTableWidget.item(i, 0) in selectedList:
                    self.submissionTableWidget.removeRow(i)


    # this function clears everything and hides the window
    def cancelFunction(self):
        # reset everything, and then hide the window
        self.organismLineEdit.setText("")
        self.selectionTableWidget.clearContents()
        self.submissionTableWidget.clearContents()
        self.selectionTableWidget.setRowCount(0)
        self.submissionTableWidget.setRowCount(0)
        self.searchProgressBar.setValue(0)
        self.hide()

    # this function downloads the selected compressed files, decompresses them, and then deletes the used
    # decompressed files. It also hides the window
    # decompressed file paths are currently stored in this window's class is a list. Can change when needed
    def submissionFunction(self):
        # clear previous searches
        self.compressedFilePaths.clear()
        self.decompressedFilePaths.clear()
        downloadedList = list()

        # basic error checking
        if self.submissionTableWidget.rowCount() <= 0 and self.selectionTableWidget.rowCount() <= 0:
            return
        elif self.submissionTableWidget.rowCount() <= 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected", "No items selected, please select some files to download",
                                           QtWidgets.QMessageBox.Ok)
            return

        # go through and download the compressed files
        for i in range(self.submissionTableWidget.rowCount()):
            tableWidgetData = self.submissionTableWidget.item(i, 0)
            for j in range(len(self.database_url_list)):
                # download ones that are selected and aren't already downloaded
                if self.idList_dict[tableWidgetData.text()] in self.database_url_list[j] and self.database_url_list[j] not in downloadedList:
                    self.compressedFilePaths.append(self.ncbi_searcher.download_compressed_file(self.database_url_list[j]))
                    downloadedList.append(self.database_url_list[j])

        # go through and decompress those files
        for i in range(len(self.compressedFilePaths)):
            self.decompressedFilePaths.append(self.ncbi_searcher.decompress_file(self.compressedFilePaths[i]))

        # go through and delete the compressed files
        file_names = os.listdir(GlobalSettings.CSPR_DB)
        for file in file_names:
            if ".gz" in file:
                print("Deleting: ", file)
                os.remove(file)

        self.hide()

#########################END OF the NCBI_Search_Window class

class NewGenome(QtWidgets.QMainWindow):
    def __init__(self, info_path):
        super(NewGenome, self).__init__()
        uic.loadUi('NewGenome.ui', self)
        self.setWindowTitle('New Genome')
        self.k = KEGG()
        self.info_path = info_path
        #---Button Modifications---#

        self.setWindowIcon(Qt.QIcon("cas9image.png"))
        self.whatsthisButton.clicked.connect(self.whatsthisclicked)
        self.KeggSearchButton.clicked.connect(self.updatekegglist)
        self.resetButton.clicked.connect(self.reset)
        self.submitButton.clicked.connect(self.submit)
        self.browseForFile.clicked.connect(self.selectFasta)
        self.NCBI_File_Search.clicked.connect(self.prep_ncbi_search)
        self.JobsQueueBox.setReadOnly(True)
        self.output_browser.setText("Waiting for program initiation...")
        self.CompletedJobs.setText(" ")

        self.comboBoxEndo.currentIndexChanged.connect(self.endo_settings)

        self.runButton.clicked.connect(self.run_jobs)
        self.clearButton.clicked.connect(self.clear_job_queue)

        self.viewStatButton.setEnabled(False)

        self.JobIndexInProgress = 0
        self.JobsQueue = []  # holds Job classes.
        self.Endos = dict()
        self.file = ""

        self.process = QtCore.QProcess()
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.finished.connect(self.upon_process_finishing)




        #show functionalities on window
        self.fillEndo()
        #self.show()




    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def selectFasta(self):

        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose a File")
        if (myFile[0] != ""):

            if not myFile[0].endswith(".fa") and not myFile[0].endswith(".fna") and not myFile[0].endswith(".gbff") and not myFile[0].endswith(".fasta"):
                QtWidgets.QMessageBox.question(self, "File Selection Error",
                                               "You have selected an incorrect type of file. "
                                               "Please choose a genbank, fasta, gbff, or a fna file.",
                                               QtWidgets.QMessageBox.Ok)
                return
            else:
                self.nameFile.setText(myFile[0])
                self.file = myFile[0]
        """cdir = self.lineEdit.text()
        os.chdir(mydir)
        self.gdirectory = mydir
        print(mydir)
        print(cdir)"""

    # this function figures out which type of file the user is searching for, and then shows the
    # ncbi_search_dialog window
    # connected to the button: self.NCBI_File_Search
    def prep_ncbi_search(self):
        fileSearchType = ""

        # check and see which type of file to search for
        if self.ref_seq_box.isChecked() and self.genbank_box.isChecked():
            fileSearchType = "GenBank"
        elif self.ref_seq_box.isChecked() and not self.genbank_box.isChecked():
            fileSearchType = "RefSeq"
        elif not self.ref_seq_box.isChecked() and self.genbank_box.isChecked():
            fileSearchType = "GenBank"
        elif not self.ref_seq_box.isChecked() and not self.genbank_box.isChecked():
            QtWidgets.QMessageBox.question(self, "Selection Error", "Please select Genbank or RefSeq", QtWidgets.QMessageBox.Ok)
            return

        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(0)
        GlobalSettings.mainWindow.ncbi_search_dialog.searchType = fileSearchType
        GlobalSettings.mainWindow.ncbi_search_dialog.organismLineEdit.setText(self.lineEdit_1.displayText())
        GlobalSettings.mainWindow.ncbi_search_dialog.show()


    def submit(self):
        warning = ""
        if len(self.lineEdit_1.text())==0:
            warning = warning + "\nYou need to include the organism's name."
        if len(self.file) == 0:
            warning = warning + "\nYou need to select a file."
        if len(warning) !=0:
            QtWidgets.QMessageBox.information(self, "Required Information", warning, QtWidgets.QMessageBox.Ok)
            return

        if len(self.lineEdit_2.text()) == 0:
            warning = warning + "\nIt is recommended to include the organism's subspecies/strain."
        if len(self.lineEdit_3.text()) == 0:
            warning = warning + "\nYou must include an organism code (KEGG code recommended)."
        if len(warning)!=0:
            hold = QtWidgets.QMessageBox.question(self, "Missing information", warning+
                                        "\n\nDo you wish to continue without including this information?"
                                       , QtWidgets.QMessageBox.Yes |
                                       QtWidgets.QMessageBox.No,
                                       QtWidgets.QMessageBox.No)
            if hold == QtWidgets.QMessageBox.No:
                return

        myjob = CasperJob(self.lineEdit_1.text() + self.lineEdit_2.text(), self.lineEdit_2.text(),
                          self.Endos[self.comboBoxEndo.currentText()], self.lineEdit_3.text(), self.file,
                          self.tot_len_box.text(), self.seed_len_box.text(), self.pamBox.isChecked())
        self.JobsQueue.append(myjob)
        nxtLine=""
        if len(self.JobsQueueBox.toPlainText())!=0:
            nxtLine = "\n"
        self.JobsQueueBox.setPlainText(self.JobsQueueBox.toPlainText()+nxtLine+myjob.name)

    def fillEndo(self):
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            f = open(self.info_path + "\\CASPERinfo")
        else:
            f = open(self.info_path + "/CASPERinfo")
        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if(line[0]=="-"):
                        break
                    line_tokened = line.split(";")
                    endo = line_tokened[0]
                    # Checking to see if there is more than one pam sequence in the list
                    if line_tokened[1].find(",") != -1:
                        p_pam = line_tokened[1].split(",")[0]
                    else:
                        p_pam = line_tokened[1]
                    default_seed_length = line_tokened[2].split(",")[0]
                    default_tot_length = line_tokened[2].split(",")[1]
                    self.Endos[endo + "PAM: " + p_pam] = (endo, p_pam, default_seed_length, default_tot_length)

                break
        f.close()
        self.comboBoxEndo.addItems(self.Endos.keys())

    def endo_settings(self):
        self.tot_len_box.setText(self.Endos[self.comboBoxEndo.currentText()][3])
        self.seed_len_box.setText(self.Endos[self.comboBoxEndo.currentText()][2])

    def findFasta(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to Quit?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            sys.exit()
        else:
            pass

    def testcheckandradio(self, state):
        if state == QtCore.Qt.Checked:
            pass

    def whatsthisclicked(self):
        QtWidgets.QMessageBox.information(self, "Organism Code", "The organism code is the manner in which CASPER will"
                                                                 "label it's data files and references for the organism"
                                                                 "you are importing here. It is HIGHLY RECOMMENDED that"
                                                                 "you use the 3-4 letter code used by KEGG as this will"
                                                                 "aid in automatic accession of annotations from the"
                                                                 "database.", QtWidgets.QMessageBox.Ok)

    def updatekegglist(self):

        self.keggSuggested.clear()
        kegg_orglist = self.k.lookfor_organism(self.lineEdit_1.text())
        holder = 0
        self.keggSuggested.setColumnCount(2)

        for item in kegg_orglist:

            second_space = item[item.find(" ") + 1:].find(" ")+item.find(" ")
            code = item[item.find(" ")+1:second_space+1]
            item = item[second_space+2:]
            semi = item.find(";")
            index = 1
            while True:
                if item[semi - index] == " ":
                    break
                index = index + 1
            organism = item[:semi - index]
            self.keggSuggested.setRowCount(holder+1)
            table_code = QtWidgets.QTableWidgetItem(code)
            table_organism = QtWidgets.QTableWidgetItem(organism)
            self.keggSuggested.setItem(holder, 0, table_organism)
            self.keggSuggested.setItem(holder, 1, table_code)
            # self.keggsearchresults.insertPlainText(item)
            holder+=1
        self.keggSuggested.resizeColumnsToContents()

    def run_jobs(self):
        # Top layer for loop to go through all of the jobs in the queue:
        job = self.JobsQueue[self.JobIndexInProgress]
        #program = GlobalSettings.CASPER_FOLDER_LOCATION + "/Casper_Seq_Finder_" + GlobalSettings.OPERATING_SYSTEM_ID
        print(GlobalSettings.appdir)
        program = GlobalSettings.appdir + '\\Casper_Seq_Finder_Windows.exe'
        self.JobInProgress.setText(job.name)
        self.process.start(program, job.get_arguments())
        self.process.readyReadStandardOutput.connect(self.output_stdout)

    def output_stdout(self):
        outputBytes = self.process.readAllStandardOutput().data()
        outputUnicode = outputBytes.decode('utf-8')
        self.output_browser.append(outputUnicode)

    def upon_process_finishing(self):
        self.CompletedJobs.append(self.JobsQueue[self.JobIndexInProgress].name)
        self.process.close()
        if self.JobIndexInProgress < len(self.JobsQueue)-1:
            self.JobIndexInProgress += 1
            self.run_jobs()

    def clear_job_queue(self):
        self.JobsQueue = []
        self.JobsQueueBox.setPlainText("")

    def reset(self):
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_3.clear()
        self.keggSuggested.clear()
        self.nameFile.setText("Name Of File")


class CasperJob:
    def __init__(self, org, suborg, endo, org_code, ref_file, tot_len, seed_len, pamdir):
        self.name = endo[0] + " targets in " + org
        self.organism_name = org
        self.substrain = suborg
        self.organism_code = org_code
        self.endo_name = endo[0]
        self.endo_pam = endo[1]
        self.reference_file = ref_file

        # These are endonuclease specific settings that should be pulled from CASPERinfo
        self.anti = pamdir
        self.sequence_length = tot_len
        self.seed_length = seed_len

    def get_arguments(self):
        ret_array = [self.endo_name, self.endo_pam, self.organism_code]
        # attach the 5' or 3' direction
        if self.anti:
            ret_array.append("TRUE")
        else:
            ret_array.append("FALSE")
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            ret_array.append(GlobalSettings.CSPR_DB + "\\")
        else:
            ret_array.append(GlobalSettings.CSPR_DB + "/")
        ret_array.append(GlobalSettings.CASPER_FOLDER_LOCATION + "/CASPERinfo")
        ret_array.append(self.reference_file)
        ret_array.append(self.organism_name)
        ret_array.append(self.sequence_length)
        ret_array.append(self.seed_length)
        ret_array.append(self.substrain)
        print(ret_array)
        return ret_array