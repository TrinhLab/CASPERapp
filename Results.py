from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
from Bio.Seq import Seq
from CSPRparser import CSPRparser
import GlobalSettings
import os
import OffTarget
import platform
import traceback

#global logger
logger = GlobalSettings.logger

# =========================================================================================
# CLASS NAME: Results
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================

class Results(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        try:
            super(Results, self).__init__(parent)
            uic.loadUi(GlobalSettings.appdir + 'resultsWindow.ui', self)
            self.setWindowTitle('Results')
            self.geneViewer.setReadOnly(True)
            self.curgene = ""
            self.dbpath = ""
            # Main Data container
            # Keys: Gene names
            # Values: #
            self.annotation_path = ""
            self.AllData = {}
            self.highlighted = {}
            self.co_target_endo_list = list()
            self.startpos = 0
            self.endpos = 0
            self.directory = ""
            self.geneDict = dict() # dictionary passed into transfer_data
            self.geneNTDict = dict() #dictionary passed into transfer_data, same key as geneDict, but hols the NTSEQ
            self.switcher = [1,1,1,1,1,1,1,1]  # for keeping track of where we are in the sorting clicking for each column

            # Initialize Filter Options Object
            self.filter_options = Filter_Options()

            # Target Table settings #
            self.targetTable.setColumnCount(8)  # hardcoded because there will always be 8 columns
            self.targetTable.setShowGrid(False)
            self.targetTable.setHorizontalHeaderLabels("Location;Endonuclease;Sequence;Strand;PAM;Score;Off-Target;Details".split(";"))
            self.targetTable.horizontalHeader().setSectionsClickable(True)
            self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.targetTable.horizontalHeader().setSectionResizeMode(7, QtWidgets.QHeaderView.Stretch) #Ensures last column goes to the edge of table

            self.back_button.clicked.connect(self.goBack)
            self.targetTable.horizontalHeader().sectionClicked.connect(self.table_sorting)
            self.off_target_button.clicked.connect(self.Off_Target_Analysis)
            self.cotargeting_button.clicked.connect(self.open_cotarget)
            self.displayGeneViewer.stateChanged.connect(self.checkGeneViewer)
            self.filter_options.cotarget_checkbox.stateChanged.connect(self.prep_cotarget_checkbox)
            self.highlight_gene_viewer_button.clicked.connect(self.highlight_gene_viewer)
            self.checkBoxSelectAll.stateChanged.connect(self.selectAll)
            self.filter_options_button.clicked.connect(self.show_filter_options)

            self.change_start_end_button.clicked.connect(self.change_indices)
            self.export_button.clicked.connect(self.open_export_to_csv)

            #self.targetTable.itemSelectionChanged.connect(self.item_select)
            self.filter_options.minScoreLine.setText("0")

            # Connecting the filters to the displayGeneData function
            self.filter_options.fivegseqCheckBox.stateChanged.connect(self.displayGeneData)
            self.filter_options.minScoreLine.textChanged.connect(self.displayGeneData)

            # Setting up the score filter:
            self.filter_options.scoreSlider.setMinimum(0)
            self.filter_options.scoreSlider.setMaximum(100)
            self.filter_options.scoreSlider.setTracking(False)

            self.filter_options.scoreSlider.valueChanged.connect(self.update_score_filter)

            #bool used to make sure only 1 instance of the OffTarget window is created
            self.first_boot = True
            #OTA is used to hold the row numbers of the items selected by user for OffTargetAnalysis
            #using this helps speed up updating the chart
            self.OTA = []

            self.clear_gene_viewer_button.clicked.connect(self.clear_gene_viewer)

            self.detail_output_list = []
            self.rows_and_seq_list = []
            self.seq_and_avg_list = []
            self.files_list = []
            self.seq_finder_cspr_file = ''
            self.mwfg = self.frameGeometry()  ##Center window
            self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window

            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#guide_viewer{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            margin-top: 10px;
                            font: bold;}"""

            self.guide_viewer.setStyleSheet(groupbox_style)
            self.guide_analysis.setStyleSheet(groupbox_style.replace("guide_viewer", "guide_analysis"))
            self.gene_viewer.setStyleSheet(groupbox_style.replace("guide_viewer", "gene_viewer"))

            #set pixel width for scroll bars
            self.targetTable.verticalScrollBar().setStyleSheet("width: 16px;")
            self.targetTable.horizontalScrollBar().setStyleSheet("height: 16px;")
            self.geneViewer.verticalScrollBar().setStyleSheet("width: 16px;")
            self.geneViewer.horizontalScrollBar().setStyleSheet("height: 16px;")

            self.get_endo_data()
        except Exception as e:
            logger.critical("Error initializing results class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def get_endo_data(self):
        try:
            f = open(GlobalSettings.appdir + "CASPERinfo")
            self.endo_data = {}
            while True:
                line = f.readline()
                if line.startswith('ENDONUCLEASES'):
                    while True:
                        line = f.readline()
                        line = line.replace("\n","")
                        if (line[0] == "-"):
                            break
                        line_tokened = line.split(";")
                        if len(line_tokened) == 10:
                            endo = line_tokened[0]
                            five_length = line_tokened[2]
                            seed_length = line_tokened[3]
                            three_length = line_tokened[4]
                            prime = line_tokened[5]
                            hsu = line_tokened[9]
                            self.endo_data[endo] = [int(five_length) + int(three_length) + int(seed_length), prime, "MATRIX:" + hsu]

                    break
            f.close()
        except Exception as e:
            logger.critical("Error in get_endo_data() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function opens the export_to_csv window
    # first it makes sure that the user actually has some highlighted targets that they want exported
    def open_export_to_csv(self):
        try:
            select_items = self.targetTable.selectedItems()
            if len(select_items) <= 0:
                QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                               "No targets were highlighted."
                                               "Please highlight the targets you want to be exported to a CSV File!",
                                               QtWidgets.QMessageBox.Ok)
                return
            # now launch the window
            GlobalSettings.mainWindow.export_csv_window.launch(select_items, 7)
        except Exception as e:
            logger.critical("Error in open_export_to_csv() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def change_indices(self):
        try:
            # make sure the gene viewer is on
            if not self.displayGeneViewer.isChecked():
                QtWidgets.QMessageBox.question(self, "Gene Viewer Error",
                                               "Gene Viewer display is off! "
                                               "Please turn the Gene Viewer on in order to highlight the sequences selected",
                                               QtWidgets.QMessageBox.Ok)
                return


            # change the start and end values
            prevTuple = self.geneDict[self.curgene]
            tempTuple = (self.geneDict[self.curgene][0], int(self.lineEditStart.displayText()), int(self.lineEditEnd.displayText()))

            # make sure that the difference between indicies is not too large
            if abs(tempTuple[1] - tempTuple[2]) > 50000:
                #print(abs(tempTuple[1] - tempTuple[2]))
                QtWidgets.QMessageBox.question(self, "Sequence Too Long",
                                               "The sequence is too long! "
                                               "Please choose indicies that will make the sequence less than 50,000!",
                                               QtWidgets.QMessageBox.Ok)
                self.lineEditStart.setText(str(self.geneDict[self.curgene][1]))
                self.lineEditEnd.setText(str(self.geneDict[self.curgene][2]))
                return

            # if the user is using gbff
            if self.annotation_path.endswith(".gbff"):
                self.geneDict[self.curgene] = tempTuple
                sequence = self.gbff_sequence_finder(self.geneDict[self.curgene])
                self.geneNTDict[self.curgene] = sequence
            # if the user is using fna (deprecated)
            """
            elif self.annotation_path.endswith(".fna"):
                self.geneDict[self.curgene] = tempTuple
                sequence = GlobalSettings.mainWindow.gene_viewer_settings.fna_sequence_finder(self.geneDict[self.curgene])
                self.geneNTDict[self.curgene] = sequence
            """
            # check and see if we need to add lowercase letters
            changeInStart = tempTuple[1] - prevTuple[1]
            changeInEnd = tempTuple[2] - prevTuple[2]

            # check and see if the sequence is extended at all
            # if it is, make the extended part lower-case as opposed to upper case
            if changeInStart != 0 and changeInStart < 0:
                tempString = self.geneNTDict[self.curgene][:abs(changeInStart)].lower()
                tempString = tempString + self.geneNTDict[self.curgene][abs(changeInStart):]
                self.geneNTDict[self.curgene] = tempString
            if changeInEnd != 0 and changeInEnd > 0:
                tempString = self.geneNTDict[self.curgene][len(self.geneNTDict[self.curgene]) - abs(changeInEnd):].lower()
                tempString = self.geneNTDict[self.curgene][:len(self.geneNTDict[self.curgene]) - abs(changeInEnd)] + tempString
                self.geneNTDict[self.curgene] = tempString

            # update the gene viewer
            self.checkGeneViewer()
        except Exception as e:
            logger.critical("Error in change_indices() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function listens for a stateChange in selectAllShown
    # if it is checked, it selects all shown
    # if it is unchecked, it deselects all shown
    # Note: it is a little buggy, possibly because when you change the minimum score it resets it all
    def selectAll(self):
        try:
            if self.checkBoxSelectAll.isChecked():
                self.targetTable.selectAll()
            else:
                self.targetTable.clearSelection()
        except Exception as e:
            logger.critical("Error in selectAll() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # hightlights the sequences found in the gene viewer
    # highlighting should stay the exact same with fasta and genbank files, as this function only edits what
    #   is currently in the gene viewer text table anyways
    def highlight_gene_viewer(self):
        try:
            # make sure gene viewer is enabled
            if not self.displayGeneViewer.isChecked():
                QtWidgets.QMessageBox.question(self, "Gene Viewer Error",
                                               "Gene Viewer display is off! "
                                               "Please turn the Gene Viewer on in order to highlight the sequences selected",
                                               QtWidgets.QMessageBox.Ok)
                return

            # variables needed
            cursor = self.geneViewer.textCursor()
            format = QtGui.QTextCharFormat()
            noMatchString = ""

            # reset the gene viewer text
            self.geneViewer.setText(self.geneNTDict[self.curgene])

            # check and make sure still is actually highlighted!
            selectedList = self.targetTable.selectedItems()
            if len(selectedList) <= 0:
                QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                               "No targets were highlighted."
                                               "Please highlight the targets you want to be highlighted in the gene viewer!",
                                               QtWidgets.QMessageBox.Ok)
                return
            # this is the loop that actually goes in and highlights all the things
            for i in range(self.targetTable.rowCount()):
                if self.targetTable.item(i, 0).isSelected():
                    # get the strand and sequence strings
                    locationString = self.targetTable.item(i,0).text()
                    strandString = self.targetTable.item(i, 3).text()
                    sequenceString = self.targetTable.item(i, 2).text()
                    printSequence = ""
                    movementIndex = 0
                    left_right = ""
                    #print("Length of geneViewer: ", len(self.geneViewer.toPlainText()))


                    # get the location
                    location = int(locationString) - self.geneDict[self.curgene][1]
                    try:
                        movementIndex = self.endo_data[self.endonucleaseBox.currentText()][0]
                    except:
                        QtWidgets.QMessageBox.critical(self, "Endo data not found.", "Could not find length of sequences in CASPERinfo file based on endo selected.", QtWidgets.QMessageBox.Ok)
                        return

                    # get which way it's moving, and the real location. This is for checking edge cases
                    if int(self.endo_data[self.endonucleaseBox.currentText()][1]) == 5:
                        # if the strand is positive, it moves to the right, if the strand is negative, it moves to the left
                        if strandString == "-":
                            left_right = "-"
                            location = (location - len(self.targetTable.item(i, 4).text())) + 1
                        elif strandString == "+":
                            location = (location + len(self.targetTable.item(i,4).text())) + 1
                            left_right = "+"
                    else:
                        # if the strand is negative, it moves to the right if the strand is positive it moves to the left
                        if strandString == "-":
                            left_right = "+"
                            #location = location + len(self.targetTable.item(i, 3).text())
                        elif strandString == "+":
                            left_right = "-"
                            #location = location - len(self.targetTable.item(i, 3).text())

                    # get the right color and the revcom
                    if strandString == "+":
                        format.setBackground(QtGui.QBrush(QtGui.QColor("green")))
                    elif strandString == "-":
                        format.setBackground(QtGui.QBrush(QtGui.QColor("red")))
                        temp = Seq(sequenceString)
                        sequenceString = temp.reverse_complement()
                        sequenceString = sequenceString.__str__()

                    testSequence = sequenceString.lower()
                    testGeneViewer = self.geneViewer.toPlainText().lower()

                    #print("Location is: ", location)
                    #print("Length of geneome viewer: ", len(self.geneViewer.toPlainText()))

                    # check to see if the sequence is in the gene viewer to behind with
                    if testSequence in testGeneViewer:
                        #print("In the if testSequence in testGeneViewer")
                        indexInViewer = testGeneViewer.find(testSequence)
                        cursor.setPosition(indexInViewer)
                        for i in range(len(sequenceString)):
                            cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                        cursor.mergeCharFormat(format)

                    # below elif's are edge cases, in case the sequence is not fully in the gene viewer
                    # if the start is too far to the left, but part of the sequence is in gene viewer
                    # and it's being built right-to-left
                    elif left_right == "-" and location - movementIndex < 0 and location > 0:
                        #print("in the first elif")
                        cursor.setPosition(0)
                        for i in range(location + 1):
                            cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                        cursor.mergeCharFormat(format)
                    # being built left-to-right
                    # if start is too far to the left, but start + total movement is in the geneviewer
                    elif left_right == "+" and location < 0 and location + movementIndex > 0:
                        #print("In the second elif")
                        cursor.setPosition(0)
                        for i in range(location + movementIndex):
                            cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                        cursor.mergeCharFormat(format)
                    # being build right-to-left
                    # if the location is too far to the right, but location-movement is in the gene viewer
                    elif left_right == "-" and location > len(self.geneViewer.toPlainText()) and location - movementIndex < len(self.geneViewer.toPlainText()):
                        #print("In the third elif statement")
                        cursor.setPosition((location - movementIndex) + 1)
                        cursor.movePosition(QtGui.QTextCursor.End, 1)
                        cursor.mergeCharFormat(format)
                    # being built left-to-right
                    # if the location + movement is too far to the right, but location is in the geneviewer
                    elif left_right == "+" and location + movementIndex > len(self.geneViewer.toPlainText()) and location < len(self.geneViewer.toPlainText()):
                        #print("In the fourth elif")
                        cursor.setPosition(location)
                        cursor.movePosition(QtGui.QTextCursor.End, 1)
                        cursor.mergeCharFormat(format)

                    # else, it is not able to be found
                    else:
                        if noMatchString == "":
                            noMatchString = sequenceString
                        else:
                            noMatchString = noMatchString + ";;" + sequenceString


            # if any of the sequences return 0 matches, show the user which ones were not found
            if len(noMatchString) >= 5:
                QtWidgets.QMessageBox.question(self, "Warning", "The following sequence(s) were not found in the Gene Viewer"
                                                                    "text:\n\t" + noMatchString, QtWidgets.QMessageBox.Ok)
        except Exception as e:
            logger.critical("Error in highlight_gene_viewer() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function updates the gene viewer based on the user clicking 'display on'
    # if it is check marked, it displays the correct data
    # if it is un-marked, it hides the data
    def checkGeneViewer(self):
        try:
            if self.displayGeneViewer.isChecked():
                self.lineEditStart.setText(str(self.geneDict[self.curgene][1]))
                self.lineEditEnd.setText(str(self.geneDict[self.curgene][2]))
                self.geneViewer.setText(self.geneNTDict[self.curgene])
            elif not self.displayGeneViewer.isChecked():
                self.lineEditStart.clear()
                self.lineEditEnd.clear()
                self.geneViewer.clear()
        except Exception as e:
            logger.critical("Error in checkGeneViewer() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function opens when the user clicks the CoTargeting button
    def open_cotarget(self):
        try:
            endo_list = list()
            if self.endonucleaseBox.count() <= 1:
                QtWidgets.QMessageBox.question(self, "Not Enough Endonucleases",
                                               "There are not enough endonucleases with this organism. "
                                               "At least 2 endonucleases are required for this function. "
                                               "Use Analyze New Genome to create CSPR files with other endonucleases.",
                                               QtWidgets.QMessageBox.Ok)
                return

            for i in range(self.endonucleaseBox.count()):
                endo_list.append(self.endonucleaseBox.itemText(i))

            GlobalSettings.mainWindow.CoTargeting.launch(endo_list, GlobalSettings.mainWindow.orgChoice.currentText())
        except Exception as e:
            logger.critical("Error in open_cotarget() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function goes through and calls transfer_data again.
    # Uses data from the mainWindow in Globalsettings, but that's because that info should not change
    # unless the user closes out of the Results window
    def changeEndonuclease(self):
        try:
            full_org = str(GlobalSettings.mainWindow.orgChoice.currentText())
            # organism = GlobalSettings.mainWindow.shortHand[full_org]

            endoChoice = self.endonucleaseBox.currentText().split("|")

            # make sure the user actually selects a new endonuclease
            if self.endo == endoChoice:
                QtWidgets.QMessageBox.question(self, "Select a different Endonuclease",
                                               "Please be sure to select a different endonuclease!",
                                               QtWidgets.QMessageBox.Ok)
                return

            # enable the cotarget checkbox if needed
            if len(endoChoice) > 1:
                self.filter_options.cotarget_checkbox.setEnabled(True)
                self.filter_options.cotarget_checkbox.setChecked(0)
            else:
                self.filter_options.cotarget_checkbox.setEnabled(False)
                self.filter_options.cotarget_checkbox.setChecked(0)
            self.transfer_data(full_org, GlobalSettings.mainWindow.organisms_to_files[full_org], endoChoice, GlobalSettings.CSPR_DB, self.geneDict,
                               self.geneNTDict, "")
        except Exception as e:
            logger.critical("Error in changeEndonuclease() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # Function that is used to set up the results page.
    # it calls get_targets, which in turn calls display data
    def transfer_data(self, org, org_files, endo, path, geneposdict, geneNTSeqDict, fasta):
        try:
            # set all of the classes variables
            self.org = org
            self.org_files = org_files
            self.endo = endo
            self.directory = path
            self.fasta_ref = fasta
            self.comboBoxGene.clear()
            self.AllData.clear()
            self.geneDict = geneposdict
            self.geneNTDict = geneNTSeqDict

            self.highlighted.clear()
            self.detail_output_list.clear()
            self.seq_and_avg_list.clear()
            self.rows_and_seq_list.clear()
            self.OTA.clear()

            for gene in geneposdict:
                detail_output1 = {}
                rows_and_seq2 = {}
                seq_and_avg3 = {}
                temp_split = gene.split(";")
                temp_len = len(temp_split)
                gene_name = temp_split[temp_len-2] + ": " + temp_split[-1]
                self.detail_output_list.append(detail_output1)
                self.seq_and_avg_list.append(seq_and_avg3)
                self.rows_and_seq_list.append(rows_and_seq2)
                self.comboBoxGene.addItem(gene_name)
                self.get_targets(gene, geneposdict[gene])

            # Enable the combobox to be toggled now that the data is in AllData
            self.comboBoxGene.currentTextChanged.connect(self.displayGeneData)
            self.first_boot = True
        except Exception as e:
            logger.critical("Error in transfer_data() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def goBack(self):
        try:
            # check and see if they searched on Sequence. If so, delete the temp CSPR file
            if len(self.seq_finder_cspr_file) > 0:
                os.remove(self.seq_finder_cspr_file)
                GlobalSettings.mainWindow.pushButton_ViewTargets.setEnabled(False)
                self.seq_finder_cspr_file = ''

            #center main on current screen
            frameGm = GlobalSettings.mainWindow.frameGeometry()
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            frameGm.moveCenter(centerPoint)
            GlobalSettings.mainWindow.move(frameGm.topLeft())

            GlobalSettings.mainWindow.show()
            self.filter_options.cotarget_checkbox.setChecked(0)
            self.hide()
        except Exception as e:
            logger.critical("Error in goBack() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # called when the user hits 'gene viewer settings'
    def changeGeneViewerSettings(self):
        try:
            GlobalSettings.mainWindow.gene_viewer_settings.show()
        except Exception as e:
            logger.critical("Error in changeGeneViewerSettings() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this is the function that sets up the cotargeting.
    # it is called from the Cotargeting class, when the user hits submit
    # myBool is whether or not to change the endoChoice comboBox
    def populate_cotarget_table(self, myBool = True):
        try:
            try:
                self.endonucleaseBox.currentIndexChanged.disconnect()
            except:
                pass

            # make a string of the combinitation, separated by commas's
            endoBoxString = ""
            for i in range(len(self.co_target_endo_list)):
                if endoBoxString == "":
                    endoBoxString = self.co_target_endo_list[i]
                else:
                    endoBoxString = endoBoxString + '|' + self.co_target_endo_list[i]

            # put the new endoChoice at the beginning. THis is the only way i could find to do it
            # get a list of all endo choices, and put the newest at the front
            endoBoxList = list()
            endoBoxList.append(endoBoxString)
            for i in range(self.endonucleaseBox.count()):
                endoBoxList.append(self.endonucleaseBox.itemText(i))

            # clear the current endo choices, and append the new order
            if myBool:
                self.endonucleaseBox.clear()
                for i in range(len(endoBoxList)):
                    self.endonucleaseBox.addItem(endoBoxList[i])

            # enable the cotarget checkbox
            self.filter_options.cotarget_checkbox.setEnabled(True)
            self.filter_options.cotarget_checkbox.setChecked(0)

            self.endonucleaseBox.currentIndexChanged.connect(self.changeEndonuclease)
            # add it to the endoBox choices, and then call transfer_data
            self.transfer_data(self.org, self.org_files, self.co_target_endo_list, GlobalSettings.CSPR_DB, self.geneDict, self.geneNTDict, "")
        except Exception as e:
            logger.critical("Error in populate_cotarget_table() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # prep function for the checkbox for cotargeting
    # if the checkbox is checked, just go ahead and displayGeneData
    # if not, call populate_cotarget_table, as a reset to get all of the data there
    def prep_cotarget_checkbox(self):
        try:
            if self.filter_options.cotarget_checkbox.isChecked():
                self.displayGeneData()
            elif not self.filter_options.cotarget_checkbox.isChecked():
                self.populate_cotarget_table(myBool=False)
        except Exception as e:
            logger.critical("Error in prep_cotarget_checkbox() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # Function grabs the information from the .cspr file and adds them to the AllData dictionary
    #changed to now call CSPRparser's function. Same function essentially, just cleaned up here
    def get_targets(self, genename, pos_tuple):
        try:
            #get the right files
            for endo in self.endo:
                if platform.system() == "Windows":
                    file = self.directory + "\\" + self.org_files[endo][0]
                else:
                    file = self.directory + "/" + self.org_files[endo][0]

                #create the parser, read the targets store it. then display the GeneData screen
                parser = CSPRparser(file)

                # if genename is not in the dict, make that spot into a list
                if genename not in self.AllData:
                    self.AllData[genename] = list()
                # now append parser's data to it
                self.AllData[genename].append(parser.read_targets(genename, pos_tuple, endo))

                # for each list item
                for item in self.AllData[genename]:
                    # for each tuple item
                    for i in range(len(item)):
                        self.highlighted[item[i][1]] = False
            self.displayGeneData()

            # if the endo choice is greater than 1, call the combine function
            if len(self.endo) > 1:
                self.combine_coTargets(genename)
        except Exception as e:
            logger.critical("Error in get_targets() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    ###############################################################################################################
    # Main Function for updating the Table.  Connected to all filter buttons and the Gene toggling of the combobox.
    ###############################################################################################################
    def displayGeneData(self):
        try:
            self.curgene = str(self.comboBoxGene.currentText())  # Gets the current gene
            # Creates the set object from the list of the current gene:
            if self.curgene=='' or len(self.AllData)<1:
                return

            # Loop through dictionary and link org dropdown to dictionary entry
            for entry in self.AllData.keys():
                if self.curgene.split(":")[0] in entry:
                    self.curgene = entry

            subset_display = set()
            # set the start and end numbers, as well as set the geneViewer text, if the displayGeneViewer is checked
            if self.displayGeneViewer.isChecked():
                self.lineEditStart.setText(str(self.geneDict[self.curgene][1]))
                self.lineEditEnd.setText(str(self.geneDict[self.curgene][2]))
                self.geneViewer.setText(self.geneNTDict[self.curgene])

            # if this checkBox is checked, remove the single endo
            if self.filter_options.cotarget_checkbox.isChecked():
                gene = self.curgene
                self.remove_single_endo(gene)

            # Removing all sequences below minimum score and creating the set:
            # for each list item
            for item in self.AllData[self.curgene]:
                # for each tuple item
                for i in range(len(item)):
                    if int(item[i][3]) > int(self.filter_options.minScoreLine.text()):
                        # Removing all non 5' G sequences:
                        if self.filter_options.fivegseqCheckBox.isChecked():
                            if item[i][1].startswith("G"):
                                subset_display.add(item[i])
                        else:
                            subset_display.add(item[i])

            self.targetTable.setRowCount(len(subset_display))

            index = 0
            #changed the number items to use setData so that sorting will work correctly
            #because before the numbers were interpretted as strings and not numbers
            for item in subset_display:
                num = int(item[0])
                loc = QtWidgets.QTableWidgetItem()
                loc.setData(QtCore.Qt.EditRole, abs(num))
                seq = QtWidgets.QTableWidgetItem(item[1])
                strand = QtWidgets.QTableWidgetItem(str(item[4]))
                PAM = QtWidgets.QTableWidgetItem(item[2])
                num1 = int(item[3])
                endonuclease = QtWidgets.QTableWidgetItem(item[5])
                score = QtWidgets.QTableWidgetItem()
                score.setData(QtCore.Qt.EditRole, num1)
                self.targetTable.setItem(index, 0, loc)
                self.targetTable.setItem(index, 1, endonuclease)
                self.targetTable.setItem(index, 2, seq)
                self.targetTable.setItem(index, 3, strand)
                self.targetTable.setItem(index, 4, PAM)
                self.targetTable.setItem(index, 5, score)
                self.targetTable.setItem(index, 6, QtWidgets.QTableWidgetItem("--.--"))
                self.targetTable.removeCellWidget(index, 7)
                if (item[1] in self.seq_and_avg_list[self.comboBoxGene.currentIndex()].keys()):
                    OT = QtWidgets.QTableWidgetItem()
                    OT.setData(QtCore.Qt.EditRole, self.seq_and_avg_list[self.comboBoxGene.currentIndex()][item[1]])
                    self.targetTable.setItem(index, 6, OT)
                if (item[1] in self.detail_output_list[self.comboBoxGene.currentIndex()].keys()):
                    details = QtWidgets.QPushButton()
                    details.setText("Details")
                    details.clicked.connect(self.show_details)
                    self.targetTable.setCellWidget(index, 7, details)
                index += 1

            self.targetTable.resizeColumnsToContents()
        except Exception as e:
            logger.critical("Error in displayGeneData() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function is only entered if the user checks the show only cotargeted sequence checkbox
    def remove_single_endo(self, gene):
        try:
            removalDict = dict()
            # go through and figure out which ones need to be shown
            for i in range(len(self.AllData[gene])):
                for j in range(len(self.AllData[gene][i])):
                    endoData = self.AllData[gene][i][j][5].split("|")
                    if len(endoData) == 1:
                        if i not in removalDict:
                            removalDict[i] = list()

                        removalDict[i].append(j)

            # now go through and delete them. But you have to go in reverse
            for item in removalDict:
                # for the reverse of that list. This is to keep the program from crashing
                # easier than building a new list honestly
                for index in reversed(removalDict[item]):
                    self.AllData[gene][item].pop(index)
        except Exception as e:
            logger.critical("Error in remove_single_endo() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function goes through and combines table rows that have the same location and PAM dir
    # it edits the dictionary data itself
    # currently it does not take the PAM direction into account
    # parameter genename:  the key to which part in the dictionary to look at
    def combine_coTargets(self, genename):
        try:
            deletingDict = dict()
            endoList = list()

            #get endo data
            self.Endos = {}
            f = open(GlobalSettings.appdir + 'CASPERinfo')
            while True:
                line = f.readline()
                if line.startswith('ENDONUCLEASES'):
                    while True:
                        line = f.readline()
                        if (line[0] == "-"):
                            break
                        line_tokened = line.split(";")
                        endo = line_tokened[0]
                        self.Endos[endo] = line_tokened[5]
                    break
            f.close()



            # get the endo list data
            for i in range(len(self.AllData[genename])):
                # if one of them is empty, just return because co-targeting is useless for that one
                if len(self.AllData[genename][i]) == 0:
                    return
                endoList.append(self.AllData[genename][i][0][5])

            # for each endoNuclease in the genename block
            for i in range(len(self.AllData[genename])):
                endoData1 = endoList[i]
                # for each target in that gene
                for j in range(len(self.AllData[genename][i])):
                    # get first locations endo
                    locationData1 = self.AllData[genename][i][j][0]
                    sequenceData1 = self.AllData[genename][i][j][1]
                    pamData1 = self.AllData[genename][i][j][2]
                    scoreData1 = self.AllData[genename][i][j][3]
                    strandData1 = self.AllData[genename][i][j][4]

                    # for each endoNuclease in the genename block
                    for k in range(len(self.AllData[genename])):
                        # if k == i then we are on the same endo target list, so break out because there can't be any combinations
                        if k == i:
                            break
                        endoData2 = endoList[k]
                        # for each target in that gene
                        for l in range(len(self.AllData[genename][k])):
                            locationData2 = self.AllData[genename][k][l][0]
                            strandData2 = self.AllData[genename][k][l][4]
                            pamData2 = self.AllData[genename][k][l][2]

                            # check which PAM is longer, and store the longer one. Otherwise, just store the first one
                            if len(pamData1) > len(pamData2):
                                storePam = pamData1
                            elif len(pamData1) < len(pamData2):
                                storePam = pamData2
                            else:
                                storePam = pamData1

                            # get the directions
                            # dir1 = self.S.endo_info[endoData1][3]
                            # dir2 = self.S.endo_info[endoData2][3]
                            dir1 = self.Endos[endoData1]
                            dir2 = self.Endos[endoData2]

                            # check if can be combined
                            if locationData1 == locationData2 and endoData1 != endoData2 and endoData2 not in endoData1:
                                if dir1 == dir2 and strandData1 == strandData2:
                                    storeEndo = self.AllData[genename][i][j][5]
                                    if endoData2 not in storeEndo:
                                        storeEndo = storeEndo + "|" + endoData2
                                    # combine the endo data
                                        self.AllData[genename][i][j] = (locationData1, sequenceData1, storePam, scoreData1, strandData1, storeEndo)

                                    # store which ones to delete
                                    if k not in deletingDict:
                                        deletingDict[k] = list()
                                    deletingDict[k].append(l)


            # delete the ones that need to be deleted
            for item in deletingDict:
                # go in the reverse of that list. This is to keep the program from crashing.
                # it is easier than building a new list honestly
                for index in reversed(deletingDict[item]):
                    self.AllData[genename][item].pop(index)
        except Exception as e:
            logger.critical("Error in combine_coTargets() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    ########################################## END UPDATING FUNCTION #############################################

    def search_gene(self):
        try:
            search_trms = []
            checkBox = self.sender()
            index = self.targetTable.indexAt(checkBox.pos())
            # print(index.column(), index.row(), checkBox.isChecked())
            seq = self.targetTable.item(index.row(),2).text()
            self.highlighted[str(seq)] = checkBox.isChecked()

            x=1
        except Exception as e:
            logger.critical("Error in search_gene() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def item_select(self):
        try:
            print(self.targetTable.selectedItems())
        except Exception as e:
            logger.critical("Error in item_select() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def table_sorting(self, logicalIndex):
        try:
            self.switcher[logicalIndex] *= -1
            if self.switcher[logicalIndex] == -1:
                self.targetTable.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.targetTable.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)
        except Exception as e:
            logger.critical("Error in table_sorting() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #linked to when the user pushes tools->off target analysis
    def Off_Target_Analysis(self):
        try:
            #build temp file for offtarget to read from
            if platform.system() == 'Windows':
                f = open(GlobalSettings.appdir + 'OffTargetFolder' + '\\temp.txt','w+')
            else:
                f = open(GlobalSettings.appdir + 'OffTargetFolder' + '/temp.txt', 'w+')
            self.OTA.clear()
            #get selected rows
            selected_rows = []
            indexes = self.targetTable.selectionModel().selectedRows()

            if len(indexes) == 0:
                QtWidgets.QMessageBox.critical(self, "No Rows Selected",
                                               "Please rows from the table to pass into the off target analysis!",
                                               QtWidgets.QMessageBox.Ok)
                return

            for index in indexes:
                selected_rows.append(index.row())

            for row in sorted(selected_rows):
                self.OTA.append(row)
                loc = self.targetTable.item(row, 0).text()
                seq = self.targetTable.item(row,2).text()
                strand = self.targetTable.item(row,3).text()
                pam = self.targetTable.item(row,4).text()
                score = self.targetTable.item(row,5).text()
                self.rows_and_seq_list[self.comboBoxGene.currentIndex()][seq] = row
                output = str(loc) + ';' + str(seq) + ";" + str(pam) + ";" + score + ";" + str(strand)
                f.write(output + '\n')
            f.close()
            #only make off target object if first time, otherwise just
            #reshow the object
            if(self.first_boot == True):
                self.first_boot = False
                self.off_tar_win = OffTarget.OffTarget()
                self.off_tar_win.submitButton.clicked.connect(self.refresh_data)
            self.off_tar_win.show()
            f.close()
        except Exception as e:
            logger.critical("Error in Off_Target_Analysis() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #refresh data is linked to the submit button on the off target analysis UI
    def refresh_data(self):
        try:
            #setup filename based on output name given in OffTarget
            filename = self.off_tar_win.output_path

            # if the user hits submit without running thr program, do nothing
            if filename == '':
                QtWidgets.QMessageBox.critical(self, "File Not Found",
                                                        "There was an error with the Off Target execution. No results file was found.",
                                                        QtWidgets.QMessageBox.Ok)
                return


            for rows in range(0,self.targetTable.rowCount()):
                self.targetTable.removeCellWidget(rows,7)
                self.targetTable.removeCellWidget(rows,5)
                self.targetTable.setItem(rows, 6, QtWidgets.QTableWidgetItem("--.--"))
            self.off_tar_win.hide()
            filename = filename[:len(filename)-1]
            filename = filename[1:]
            filename = filename.replace(r'\\', '\\')
            filename = filename.replace('"', '')
            self.files_list.append(filename)
            try:
                out_file = open(filename, "r")
            except:
                QtWidgets.QMessageBox.critical(self, "Unable to Open File",
                                                        "There was an error with the Off Target execution. The results file was either not found or not able to be opened.",
                                                        QtWidgets.QMessageBox.Ok)
                return

            #read the first line : either AVG or DETAILED OUTPUT
            output_type = out_file.readline()
            output_type = output_type.strip('\r\n')
            #parse based on whether avg or detailed output
            line_cnt = 0
            if(output_type == "AVG OUTPUT"):
                for line in out_file:
                    line = line.strip('\n')
                    if (line != ''):
                        values = line.split(":")
                        row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                        OT = QtWidgets.QTableWidgetItem()
                        OT.setData(QtCore.Qt.EditRole, values[1])
                        self.targetTable.setItem(row, 6, OT)
                        self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
                        line_cnt += 1
            else:
                details_bool = False
                temp_list = []
                values = []
                for line in out_file:
                    line = line.strip('\n')
                    if(line.find(':') != -1):
                        if(details_bool == True):
                            self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list
                            details_bool = False
                            temp_list = []
                        values = line.split(":")
                        self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
                        row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                        OT = QtWidgets.QTableWidgetItem()
                        OT.setData(QtCore.Qt.EditRole, values[1])
                        self.targetTable.setItem(row, 6, OT)
                        line_cnt += 1
                    elif line != "":
                        details_bool = True
                        temp_list.append(line)
                        details = QtWidgets.QPushButton()
                        details.setText("Details")
                        details.clicked.connect(self.show_details)
                        self.targetTable.setCellWidget(row, 7, details)
                        line_cnt += 1
                if(details_bool == True):
                    self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list


                #make sure OT output file had lines
                if line_cnt < 1:
                    QtWidgets.QMessageBox.critical(self, "File Empty",
                                                        "There was an error with the Off Target execution. No results were found in the results file.",
                                                        QtWidgets.QMessageBox.Ok)
                    return

                #print(self.detail_output_list)
                self.targetTable.resizeColumnsToContents()
                out_file.close()
        except Exception as e:
            logger.critical("Error in refresh_data() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def show_details(self):
        try:
            #create msg box popup for the details | used html to make it easier to style with bold
            button = self.sender()
            index = self.targetTable.indexAt(button.pos())
            msg = QtWidgets.QMessageBox()
            msg.setWindowTitle("Details")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            key = str(self.targetTable.item(index.row(),2).text())
            temp_str = ''
            for items in self.detail_output_list[self.comboBoxGene.currentIndex()][key]:
                temp_str += items + "<br>"

            chromo_str = "<html><b>Chromosome: Location, Sequence, Strand, PAM, Score:<br></b></html>"
            input_str = self.targetTable.item(index.row(),0).text() + ' , ' + key + ' , ' + \
                        self.targetTable.item(index.row(),3).text() + ' , ' + self.targetTable.item(index.row(),4).text() + \
                        ' , ' + self.targetTable.item(index.row(),5).text() + "<br><br>"
            detail_str = "<html><b>Deatailed Output: Score, Chromosome, Location, Sequence:<br></b></html>"
            msg.setText(chromo_str + input_str + detail_str + temp_str)
            msg.exec()
        except Exception as e:
            logger.critical("Error in show_details() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #function to clear highlights in gene viewer and un-select any rows highlighted in main table
    def clear_gene_viewer(self):
        try:
            #clear gene viewer highlights
            if self.displayGeneViewer.isChecked():
                self.geneViewer.setText(self.geneNTDict[self.curgene])

            #clear selected rows in table
            self.checkBoxSelectAll.setChecked(False)
            self.targetTable.selectionModel().clearSelection()
        except Exception as e:
            logger.critical("Error in clear_gene_viewer() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # -----------------------------------------------------------------------------------------------------#
    # ---- All Filter functions below ---------------------------------------------------------------------#
    # -----------------------------------------------------------------------------------------------------#
    def update_score_filter(self):
        try:
            self.filter_options.minScoreLine.setText(str(self.filter_options.scoreSlider.value()))
        except Exception as e:
            logger.critical("Error in update_score_filter() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    ###Save data and open data functions are currently deprecated
    """
    #allows user to save what is currently in the table
    def save_data(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,
                                      "Enter Text File Name", ".txt",
                                      "Text Document (*.txt)" )

        if (str(filename[0]) != ''):
            f = open(str(filename[0]), "w+")
            for genomes in self.AllData:
                f.write("***"+genomes + "\n")
                for items in self.AllData[genomes]:
                    for i in items:
                        for j in i:
                            f.write(str(j) + '|')

                        f.write(str(self.highlighted[items[1][1]]))
                        f.write("\n")
            f.close()

    #open any saved .txt of previous tables opened
    def open_data(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File')
        if (os.path.isfile(str(filename[0]))):
            f = open(str(filename[0]), "r+")
            temp_str = f.readline()
            if(temp_str.startswith("***") == False):
                msg = QtWidgets.QMessageBox()
                msg.setWindowTitle("Error")
                msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("<font size=4>" + "Not correct file format" + "</font>")
                msg.exec()
                f.close()
            else:
                f.close()
                self.AllData.clear()
                self.highlighted.clear()
                first = 1
                list1 = []
                list2 = []
                s = ""
                self.comboBoxGene.clear()
                if (os.path.isfile(str(filename[0]))):
                    f = open(str(filename[0]), "r+")
                    for line in f:
                        if(line.startswith("***")):
                            if(first == 0):
                                list2.append(list1)
                                self.AllData[s] = list2
                            else:
                                first = 0
                            s = line[3:]
                            s = s.strip('\n')
                            list1 = []
                            self.comboBoxGene.addItem(s)
                        else:
                            temp = line.split("|")
                            h = temp.pop()
                            h = h.strip('\n')
                            self.highlighted[temp[1]] = eval(h)
                            tup = tuple(temp)
                            list1.append(tup)
                    list2.append(list1)
                    self.AllData[s] = list2
                    f.close()
                    self.displayGeneData()

    """

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        try:
            # check and see if they searched on Sequence. If so, delete the temp CSPR file
            if len(self.seq_finder_cspr_file) > 0:
                os.remove(self.seq_finder_cspr_file)
                GlobalSettings.mainWindow.pushButton_ViewTargets.setEnabled(False)
                self.seq_finder_cspr_file = ''

            GlobalSettings.mainWindow.closeFunction()
            event.accept()
        except Exception as e:
            logger.critical("Error in closeEvent() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    ############################################
    # All Gene Viewer functions are below!
    ############################################
    # It will go based on the lengths stored in the comboGeneBox dictionary
    def load_gene_viewer(self):
        try:
            sequence = ""
            # for each gene selected from the results window
            for item in self.geneDict:
            ### FNA support deprecated currently
                """
                if self.file_type == "fna":
                    sequence = self.fna_sequence_finder(GlobalSettings.mainWindow.Results.geneDict[item])
                    GlobalSettings.mainWindow.Results.geneNTDict[item] = sequence
                """
                if self.annotation_path.endswith(".gbff"):
                    sequence = self.gbff_sequence_finder(self.geneDict[item])
                    self.geneNTDict[item] = sequence
            self.lineEditStart.setEnabled(True)
            self.lineEditEnd.setEnabled(True)
            self.change_start_end_button.setEnabled(True)
            self.displayGeneViewer.setChecked(0)
            self.checkGeneViewer()
        except Exception as e:
            logger.critical("Error in load_gene_viewer() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function gets the sequence out of the GBFF file
    # may have indexing issues
    def gbff_sequence_finder(self, location_data):
        try:
            # start up the function
            fileStream = open(self.annotation_path)
            buffer = fileStream.readline()
            index = 1
            pre_sequence = ""

            # get to the first chromosome's origin
            while True:
                if "ORIGIN" in buffer:
                    buffer = fileStream.readline()
                    break
                buffer = fileStream.readline()

            # skip all of the data until we are at the chromosome we care about
            while index != location_data[0]:
                if "ORIGIN" in buffer:
                    index += 1
                buffer = fileStream.readline()

            # get the entire chromesome into a string
            while "//" not in buffer:
                if "LOCUS" in buffer or buffer == "":
                    break
                # replace digits with spaces (if i replace them with nothing the program will crash)
                for i in range(len(buffer)):
                    if buffer[i].isdigit():
                        buffer = buffer.replace(buffer[i], " ")

                # replace all of the spaces
                buffer = buffer.replace(" ", "")

                # append it to the stored version of the entire string
                if pre_sequence == "":
                    pre_sequence = buffer
                else:
                    pre_sequence = pre_sequence + buffer

                buffer = fileStream.readline()

            # take out the endlines and uppercase the string
            pre_sequence = pre_sequence.replace("\n", "")
            pre_sequence = pre_sequence.upper()
            #print("Length of the pre-sequence: ", len(pre_sequence))

            ###Get gene sequence and padding sequences (for visualizing gRNAs that appear at extreme ends of gene)
            if location_data[1] - 30 >= 0: ### Check to make sure there is enough 5' end of gene to pull the padding from, so indexing error isn't raised
                five_prime_tail = str(pre_sequence[location_data[1]-31:location_data[1]-1])
            else:
                five_prime_tail = ""
            if len(pre_sequence) > (location_data[2] + 30): ### Check to make sure there is enough 3' end of gene to pull the padding from, so indexing error isn't raised
                three_prime_tail = str(pre_sequence[location_data[2]:location_data[2]+30])
            else:
                three_prime_tail = ""

            gene_sequence = str(pre_sequence[location_data[1]-1:location_data[2]])
            ret_sequence = five_prime_tail.lower() + gene_sequence.upper() + three_prime_tail.lower()
            return ret_sequence

        except Exception as e:
            logger.critical("Error in gbff_sequence_finder() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def show_filter_options(self):
        try:
            self.filter_options.show()
        except Exception as e:
            logger.critical("Error in show_filter_options() in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

class Filter_Options(QtWidgets.QDialog):
    def __init__(self, parent=None):
        try:
            super(Filter_Options, self).__init__(parent)
            uic.loadUi(GlobalSettings.appdir + 'filter_options.ui', self)
            self.minScoreLine.setText("0")
            self.hide()
            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#filterBox{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            margin-top: 10px;}"""
            self.filterBox.setStyleSheet(groupbox_style)
        except Exception as e:
            logger.critical("Error initializing Filter_Options class in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)