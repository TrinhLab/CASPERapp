import GlobalSettings
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
from NCBI_API import Assembly, GBFF_Parse
import os

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
        # check to see what database to search in
        if self.refseq_checkbox.isChecked():
            self.searchType = 'RefSeq'
        elif self.genbank_checkbox.isChecked():
            self.searchType = 'GenBank'
        elif not self.refseq_checkbox.isChecked() and not self.genbank_checkbox.isChecked():
            QtWidgets.QMessageBox.question(self, "Selection Error", "Please select Genbank or RefSeq databases",
                                           QtWidgets.QMessageBox.Ok)
            return

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
        self.refseq_checkbox.setChecked(False)
        self.genbank_checkbox.setChecked(False)
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
        file_names = [f for f in os.listdir(GlobalSettings.CSPR_DB) if os.path.isfile(os.path.join(GlobalSettings.CSPR_DB, f))]
        for file in file_names:
            if ".gz" in file:
                #print("Deleting: ", file)
                os.remove(file)

        GlobalSettings.pop_Analysis.get_data()

        self.cancelFunction()
