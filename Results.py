import sys


from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
from bs4 import BeautifulSoup
import requests
import webbrowser
from Scoring import OnTargetScore
import Algorithms
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
import GlobalSettings
import os
from APIs import Kegg
import OffTarget

# =========================================================================================
# CLASS NAME: Results
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class Results(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(Results, self).__init__(parent)
        uic.loadUi('resultsWindow.ui', self)

        self.setWindowTitle('Results')
        self.geneViewer.setReadOnly(True)
        # Scoring Class object #
        self.onscore = OnTargetScore()
        self.S = SeqTranslate()
        self.dbpath = ""
        # Main Data container
        # Keys: Gene names
        # Values: #
        self.AllData = {}
        self.highlighted = {}

        self.startpos = 0
        self.endpos = 0
        self.directory = ""
        self.geneDict = dict() # dictionary passed into transfer_data
        self.geneNTDict = dict() #dictionary passed into transfer_data, same key as geneDict, but hols the NTSEQ

        self.switcher = [1,1,1,1,1,1,1]  # for keeping track of where we are in the sorting clicking for each column

        # Target Table settings #
        self.targetTable.setColumnCount(8)  # hardcoded because there will always be seven columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Location;Sequence;Strand;PAM;Score;Off-Target;Off-Target;Details".split(";"))
        self.targetTable.horizontalHeader().setSectionsClickable(True)
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        self.back_button.clicked.connect(self.goBack)
        self.targetTable.horizontalHeader().sectionClicked.connect(self.table_sorting)
        self.actionSave.triggered.connect(self.save_data)
        self.actionOpen.triggered.connect(self.open_data)
        self.actionOff_Target_Analysis.triggered.connect(self.Off_Target_Analysis)
        self.actionCoTargeting.triggered.connect(self.open_cotarget)
        self.changeEndoButton.clicked.connect(self.changeEndonuclease)
        self.displayGeneViewer.stateChanged.connect(self.checkGeneViewer)
        self.highlight_gene_viewer_button.clicked.connect(self.highlight_gene_viewer)
        self.checkBoxSelectAll.stateChanged.connect(self.selectAll)
        self.pushButton_Deselect_All.clicked.connect(self.deselectAll)
        self.gene_viewer_settings_button.clicked.connect(self.changeGeneViewerSettings)
        self.change_start_end_button.clicked.connect(self.change_indicies)


        #self.targetTable.itemSelectionChanged.connect(self.item_select)
        self.minScoreLine.setText("0")

        # Connecting the filters to the displayGeneData function
        self.fivegseqCheckBox.stateChanged.connect(self.displayGeneData)
        self.minScoreLine.textChanged.connect(self.displayGeneData)

        # Setting up the score filter:
        self.scoreSlider.setMinimum(0)
        self.scoreSlider.setMaximum(100)
        self.scoreSlider.setTracking(False)

        self.scoreSlider.valueChanged.connect(self.update_score_filter)

        #bool used to make sure only 1 instance of the OffTarget window is created
        self.first_boot = True
        #OTA is used to hold the row numbers of the items selected by user for OffTargetAnalysis
        #using this helps speed up updating the chart
        self.OTA = []

        self.selectAllButton.clicked.connect(self.selectAll_OffT)


        self.detail_output_list = []
        self.rows_and_seq_list = []
        self.seq_and_avg_list = []
        self.files_list = []

    def change_indicies(self):

        # make sure the gene viewer is on
        if not self.displayGeneViewer.isChecked():
            QtWidgets.QMessageBox.question(self, "Gene Viewer Error",
                                           "Gene Viewer display is off! "
                                           "Please turn the Gene Viewer on in order to highlight the sequences selected",
                                           QtWidgets.QMessageBox.Ok)
            return


        # change the start and end values
        tempTuple = (self.geneDict[self.comboBoxGene.currentText()][0], int(self.lineEditStart.displayText()), int(self.lineEditEnd.displayText()))
        self.geneDict[self.comboBoxGene.currentText()] = tempTuple

        # if the user is using kegg
        if GlobalSettings.mainWindow.gene_viewer_settings.file_type == "kegg":
            # build the URL string
            url = "https://www.genome.jp/dbget-bin/cut_sequence_genes.pl?FROM="
            url = url + str(self.geneDict[self.comboBoxGene.currentText()][1])
            url = url + "&TO="
            url = url + str(self.geneDict[self.comboBoxGene.currentText()][2])
            url = url + "&VECTOR=1&ORG="
            url = url + GlobalSettings.mainWindow.Annotations_Organism.currentText().split(" ")[1]

            # soup time
            source = requests.get(url)
            plain_text = source.text
            soup = BeautifulSoup(plain_text, "html.parser")


            buffer = ""

            # use the try because if user gives bad input
            try:
                for item in soup.find('pre'):
                    # get the first line
                    if ">" in item.string:
                        temp = item.string.replace("\n", "")
                        temp = temp.replace(" ", "")
                        newLineIndex = temp.find(")") + 1

                        if buffer == "":
                            buffer = temp[newLineIndex:]
                        else:
                            buffer = buffer + temp[newLineIndex:]
                    # get every other lie
                    else:
                        if buffer == "":
                            buffer = item.string
                        else:
                            buffer = buffer + item.string

                buffer = buffer.replace("\n", "")
                self.geneNTDict[self.comboBoxGene.currentText()] = buffer
            except:
                print("indexing error")

        # if the user is using gbff
        elif GlobalSettings.mainWindow.gene_viewer_settings.file_type == "gbff":
            sequence = GlobalSettings.mainWindow.gene_viewer_settings.gbff_sequence_finder(self.geneDict[self.comboBoxGene.currentText()])
            self.geneNTDict[self.comboBoxGene.currentText()] = sequence
        # if the user is using fna
        elif GlobalSettings.mainWindow.gene_viewer_settings.file_type == "fna":
            sequence = GlobalSettings.mainWindow.gene_viewer_settings.fna_sequence_finder(self.geneDict[self.comboBoxGene.currentText()])
            self.geneNTDict[self.comboBoxGene.currentText()] = sequence

        self.checkGeneViewer()

    # this function goes through and deselects all of the Off-Target checkboxes
    # also sets the selectAllShown check box to unchecked as well
    def deselectAll(self):
        for i in range(self.targetTable.rowCount()):
            self.targetTable.cellWidget(i, 6).setCheckState(0)
        self.checkBoxSelectAll.setCheckState(0)

    # this function listens for a stateChange in selectAllShown
    # if it is checked, it selects all shown
    # if it is unchecked, it deselects all shown
    # Note: it is a little buggy, possibly because when you change the minimum score it resets it all
    def selectAll(self):
        if self.checkBoxSelectAll.isChecked():
            for i in range(self.targetTable.rowCount()):
                self.targetTable.cellWidget(i, 6).setCheckState(2)
        elif not self.checkBoxSelectAll.isChecked():
            for i in range(self.targetTable.rowCount()):
                self.targetTable.cellWidget(i, 6).setCheckState(0)

    # hightlights the sequences found in the gene viewer
    # highlighting should stay the exact same with fasta and genbank files, as this function only edits what
    #   is currently in the gene viewer text table anyways
    def highlight_gene_viewer(self):
        # make sure gene viewer is enabled
        if not self.displayGeneViewer.isChecked():
            QtWidgets.QMessageBox.question(self, "Gene Viewer Error",
                                           "Gene Viewer display is off! "
                                           "Please turn the Gene Viewer on in order to highlight the sequences selected",
                                           QtWidgets.QMessageBox.Ok)
            return

        # variables needed
        k = Kegg()
        cursor = self.geneViewer.textCursor()
        format = QtGui.QTextCharFormat()
        noMatchString = ""

        # reset the gene viewer text
        self.geneViewer.setText(self.geneNTDict[self.comboBoxGene.currentText()])

        # check and make sure still is actually highlighted!
        selectedList = self.targetTable.selectedItems()
        if len(selectedList) <= 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "No targets were highlighted "
                                           "Please highlight the targets you want to be highlighted in the gene viewer!",
                                           QtWidgets.QMessageBox.Ok)
            return
        # this is the loop that actually goes in and highlights all the things
        for i in range(self.targetTable.rowCount()):
            if self.targetTable.item(i, 0).isSelected():
                # get the strand and sequence strings
                locationString = self.targetTable.item(i,0).text()
                strandString = self.targetTable.item(i, 2).text()
                sequenceString = ""
                printSequence = ""
                movementIndex = 0
                left_right = ""


                location = int(locationString) - self.geneDict[self.comboBoxGene.currentText()][1]

                print("Location is: ", location)

                # if the endo is Cas12
                if "Cas12" in self.endonucleaseBox.currentText():
                    # movement is always 24
                    movementIndex = 24

                    # if the strand is positive, it moves to the right, if the strand is negative, it moves to the left
                    if strandString == "-":
                        left_right = "-"
                        location = location - len(self.targetTable.item(i, 3).text())
                    elif strandString == "+":
                        location = (location + len(self.targetTable.item(i,3).text())) + 1
                        left_right = "+"
                # if the endo is Cas9
                elif "Cas9" in self.endonucleaseBox.currentText():
                    # movement is always 20
                    movementIndex = 20

                    # if the strand is negative, it moves to the right if the strand is positive it moves to the left
                    if strandString == "-":
                        left_right = "+"
                    elif strandString == "+":
                        left_right = "-"

                if left_right == "+":
                    printSequence = self.geneViewer.toPlainText()[location:location + movementIndex]
                elif left_right == "-":
                    for i in range(movementIndex):
                        if printSequence == "":
                            printSequence = self.geneViewer.toPlainText()[location + 1]
                        else:
                            printSequence = printSequence + self.geneViewer.toPlainText()[(location - i) + 1]

                print(printSequence)

                sequenceString = printSequence

                # check whether to be red or green
                if strandString == "+":
                    format.setBackground(QtGui.QBrush(QtGui.QColor("green")))
                elif strandString == "-":
                    format.setBackground(QtGui.QBrush(QtGui.QColor("red")))
                    sequenceString = sequenceString[::-1]

                # go through and highlight if it's in the geneviewer
                if sequenceString in self.geneViewer.toPlainText():
                    regex = QtCore.QRegExp(sequenceString)
                    index = regex.indexIn(self.geneViewer.toPlainText(), 0)
                    cursor.setPosition(index)
                    for i in range(len(sequenceString)):
                        cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                    cursor.mergeCharFormat(format)
                # catch the error if it is not in gene viewer
                else:
                    if noMatchString == "":
                        noMatchString = sequenceString
                    else:
                        noMatchString = noMatchString + ";;" + sequenceString

        # if any of the sequences return 0 matches, show the user which ones were not found
        if len(noMatchString) >= 5:
            QtWidgets.QMessageBox.question(self, "Warning", "The following sequence(s) were not found in the Gene Viewer"
                                                                "text:\n\t" + noMatchString, QtWidgets.QMessageBox.Ok)

    # this function updates the gene viewer based on the user clicking 'display on'
    # if it is check marked, it displays the correct data
    # if it is un-marked, it hides the data
    def checkGeneViewer(self):
        if self.displayGeneViewer.isChecked():
            self.lineEditStart.setText(str(self.geneDict[self.comboBoxGene.currentText()][1]))
            self.lineEditEnd.setText(str(self.geneDict[self.comboBoxGene.currentText()][2]))
            self.geneViewer.setText(self.geneNTDict[self.comboBoxGene.currentText()])
        elif not self.displayGeneViewer.isChecked():
            self.lineEditStart.clear()
            self.lineEditEnd.clear()
            self.geneViewer.clear()

    # this function just opens CoTargeting when the user clicks the CoTargeting button
    # opened the same way that main opens it
    def open_cotarget(self):
        GlobalSettings.mainWindow.CoTargeting.launch(GlobalSettings.mainWindow.data, GlobalSettings.CSPR_DB, GlobalSettings.mainWindow.shortHand)

    # this function goes through and calls transfer_data again.
    # Uses data from the mainWindow in Globalsettings, but that's because that info should not change
    # unless the user closes out of the Results window
    def changeEndonuclease(self):
        full_org = str(GlobalSettings.mainWindow.orgChoice.currentText())
        organism = GlobalSettings.mainWindow.shortHand[full_org]
        self.transfer_data(organism, str(self.endonucleaseBox.currentText()), GlobalSettings.CSPR_DB, self.geneDict,
                           self.geneNTDict, "")

    # Function that is called in main in order to pass along the information the user inputted and the information
    # from the .cspr files that was discovered
    def transfer_data(self, org, endo, path, geneposdict, geneNTSeqDict, fasta):
        self.org = org
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
            self.detail_output_list.append(detail_output1)
            self.seq_and_avg_list.append(seq_and_avg3)
            self.rows_and_seq_list.append(rows_and_seq2)
            self.comboBoxGene.addItem(gene)
            self.get_targets(gene, geneposdict[gene])



        # Enable the combobox to be toggled now that the data is in AllData
        self.comboBoxGene.currentTextChanged.connect(self.displayGeneData)
        self.first_boot = True

    def goBack(self):
        GlobalSettings.mainWindow.show()
        self.hide()

    def changeGeneViewerSettings(self):
        GlobalSettings.mainWindow.gene_viewer_settings.show()

    # Function grabs the information from the .cspr file and adds them to the AllData dictionary
    #changed to now call CSPRparser's function. Same function essentially, just cleaned up here
    def get_targets(self, genename, pos_tuple):
        #get the right files
        if self.directory.find("/") != -1:
            file = (self.directory+"/" + self.org + "_" + self.endo + ".cspr")
        else:
            file = (self.directory + "\\" + self.org + "_" + self.endo + ".cspr")

        #create the parser, read the targets store it. then display the GeneData screen
        parser = CSPRparser(file)
        temp = parser.read_targets(genename, pos_tuple)
        self.AllData[genename] = parser.read_targets(genename, pos_tuple)
        for item in self.AllData[genename]:
            self.highlighted[item[1]] = False
        self.displayGeneData()


    ###############################################################################################################
    # Main Function for updating the Table.  Connected to all filter buttons and the Gene toggling of the combobox.
    ###############################################################################################################
    def displayGeneData(self):

        curgene = str(self.comboBoxGene.currentText())  # Gets the current gene
        # Creates the set object from the list of the current gene:
        if curgene=='' or len(self.AllData)<1:
            return

        subset_display = set()

        # set the start and end numbers, as well as set the geneViewer text, if the displayGeneViewer is checked
        if self.displayGeneViewer.isChecked():
            self.lineEditStart.setText(str(self.geneDict[self.comboBoxGene.currentText()][1]))
            self.lineEditEnd.setText(str(self.geneDict[self.comboBoxGene.currentText()][2]))
            self.geneViewer.setText(self.geneNTDict[self.comboBoxGene.currentText()])

        # Removing all sequences below minimum score and creating the set:

        for item in self.AllData[curgene]:
            if int(item[3]) > int(self.minScoreLine.text()):
                # Removing all non 5' G sequences:
                if self.fivegseqCheckBox.isChecked():
                    if item[1].startswith("G"):
                        subset_display.add(item)
                else:
                    subset_display.add(item)

        self.targetTable.setRowCount(len(subset_display))

        index = 0
        #changed the number items to use setData so that sorting will work correctly
        #because before the numbers were interpretted as strings and not numbers
        for item in subset_display:
            num = int(item[0])
            loc = QtWidgets.QTableWidgetItem()
            loc.setData(QtCore.Qt.EditRole, num)
            seq = QtWidgets.QTableWidgetItem(item[1])
            strand = QtWidgets.QTableWidgetItem(str(item[4]))
            PAM = QtWidgets.QTableWidgetItem(item[2])
            num1 = int(item[3])
            score = QtWidgets.QTableWidgetItem()
            score.setData(QtCore.Qt.EditRole, num1)
            self.targetTable.setItem(index, 0, loc)
            self.targetTable.setItem(index, 1, seq)
            self.targetTable.setItem(index, 2, strand)
            self.targetTable.setItem(index, 3, PAM)
            self.targetTable.setItem(index, 4, score)
            self.targetTable.setItem(index, 5, QtWidgets.QTableWidgetItem("--.--"))
            self.targetTable.removeCellWidget(index, 7)
            ckbox = QtWidgets.QCheckBox()
            ckbox.clicked.connect(self.search_gene)
            self.targetTable.setCellWidget(index,6,ckbox)
            if(self.highlighted[str(item[1])] == True):
                ckbox.click()
            if (item[1] in self.seq_and_avg_list[self.comboBoxGene.currentIndex()].keys()):
                OT = QtWidgets.QTableWidgetItem()
                OT.setData(QtCore.Qt.EditRole, self.seq_and_avg_list[self.comboBoxGene.currentIndex()][item[1]])
                self.targetTable.setItem(index, 5, OT)
            if (item[1] in self.detail_output_list[self.comboBoxGene.currentIndex()].keys()):
                details = QtWidgets.QPushButton()
                details.setText("Details")
                details.clicked.connect(self.show_details)
                self.targetTable.setCellWidget(index, 7, details)
            index += 1



        self.targetTable.resizeColumnsToContents()

    def search_gene(self):
        search_trms = []
        checkBox = self.sender()
        index = self.targetTable.indexAt(checkBox.pos())
        # print(index.column(), index.row(), checkBox.isChecked())
        seq = self.targetTable.item(index.row(),1).text()
        self.highlighted[str(seq)] = checkBox.isChecked()

        x=1

    def item_select(self):
        print(self.targetTable.selectedItems())

    def table_sorting(self, logicalIndex):
        self.switcher[logicalIndex] *= -1
        if self.switcher[logicalIndex] == -1:
            self.targetTable.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.targetTable.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)

        #update the OTA list in case user already has stuff selected, because their rows
        #will most likely change after the sort
        self.OTA.clear()
        for row in range(self.targetTable.rowCount()):
            if (self.targetTable.cellWidget(row, 6).isChecked()):
                self.OTA.append(row)

    #linked to when the user pushes tools->off target analysis
    def Off_Target_Analysis(self):
        #build temp file for offtarget to read from
        f = open(GlobalSettings.appdir + '\\OffTargetFolder' + '\\temp.txt','w+')
        self.OTA.clear()
        for row in range(self.targetTable.rowCount()):
            if(self.targetTable.cellWidget(row,6).isChecked()):
                self.OTA.append(row)
                loc = self.targetTable.item(row, 0).text()
                seq = self.targetTable.item(row,1).text()
                strand = self.targetTable.item(row,2).text()
                pam = self.targetTable.item(row,3).text()
                score = self.targetTable.item(row,4).text()
                loc = self.S.compress(int(loc), 64)
                self.rows_and_seq_list[self.comboBoxGene.currentIndex()][seq] = row
                seq = self.S.compress(seq, 64)
                pam = self.S.compress(pam, 64)
                score = self.S.compress(int(score), 64)
                output = str(loc) + ',' + str(seq) + str(strand) + str(pam) + ',' + score
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

    #refresh data is linked to the submit button on the off target analysis UI
    def refresh_data(self):
        self.off_tar_win.hide()
        #setup filename based on output name given in OffTarget
        filename = self.off_tar_win.output_path
        filename = filename[:len(filename)-1]
        filename = filename[1:]
        filename = filename.replace(r'\\', '\\')
        filename = filename.replace('"', '')
        self.files_list.append(filename)
        out_file = open(filename, "r")
        #read the first line : either AVG or DETAILED OUTPUT
        output_type = out_file.readline()
        output_type = output_type.strip('\r\n')
        length = len(self.OTA)
        #parse based on whether avg or detailed output
        if(output_type == "AVG OUTPUT"):
            for i in range(0,length):
                line = out_file.readline()
                line = line.strip('\n')
                if (line != ''):
                    values = line.split(":")
                    row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                    OT = QtWidgets.QTableWidgetItem()
                    OT.setData(QtCore.Qt.EditRole, values[1])
                    self.targetTable.setItem(row, 5, OT)
                    self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
        else:
            for i in range(0,length):
                temp_list = []
                line = out_file.readline()
                line = line.strip('\n')
                if(line != ''):
                    values = line.split(":")
                    row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                    OT = QtWidgets.QTableWidgetItem()
                    OT.setData(QtCore.Qt.EditRole, values[1])
                    self.targetTable.setItem(row, 5, OT)
                    line = out_file.readline()
                    line = line.strip('\n')
                    if(line != ''):
                        temp_list.append(line)
                        details = QtWidgets.QPushButton()
                        details.setText("Details")
                        details.clicked.connect(self.show_details)
                        self.targetTable.setCellWidget(row, 7, details)
                        while(True):
                            line = out_file.readline()
                            line = line.strip('\n')
                            if(line == ''):
                                self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list
                                break
                            temp_list.append(line)
                    self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
        self.targetTable.resizeColumnsToContents()
        out_file.close()


    def show_details(self):
        #create msg box popup for the details | used html to make it easier to style with bold
        button = self.sender()
        index = self.targetTable.indexAt(button.pos())
        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle("Details")
        msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
        key = str(self.targetTable.item(index.row(),1).text())
        temp_str = ''
        for items in self.detail_output_list[self.comboBoxGene.currentIndex()][key]:
            temp_str += items + "<br>"

        chromo_str = "<html><b>Chromsome: Location, Sequence, Strand, PAM, Score:<br></b></html>"
        input_str = self.targetTable.item(index.row(),0).text() + ' , ' + key + ' , ' + \
                    self.targetTable.item(index.row(),2).text() + ' , ' + self.targetTable.item(index.row(),3).text() + \
                    ' , ' + self.targetTable.item(index.row(),4).text() + "<br><br>"
        detail_str = "<html><b>Deatailed Output: Score, Chromsome, Location, Sequence:<br></b></html>"
        msg.setText(chromo_str + input_str + detail_str + temp_str)
        msg.exec()

    #select all off target checkboxes in results window
    def selectAll_OffT(self):
        for row in range(self.targetTable.rowCount()):
            if(self.targetTable.cellWidget(row,6).isChecked() == False):
                self.targetTable.cellWidget(row,6).click()



    # Function for displaying the target in the gene viewer
    """def displayGene(self,fastafile=None, Kegg=False, NCBI=False):
        organism_genome = list()  # list of chromosomes/scaffolds
        if fastafile:
            f = open(fastafile)
            chr_string = str()
            for line in f:
                if not line.startswith(">"):
                    chr_string += line[:-1]
                else:
                    organism_genome.append(chr_string)
                    chr_string = ""
            return organism_genome
        elif Kegg:
            # Get the gene from the Kegg database
        elif NCBI:
            # Get the gene from NCBI database (RefSeq)
        else:
            return "Error: Cannot find reference sequence.  Search Kegg, NCBI, or download a FASTA file to create a genome reference."""""

    # -----------------------------------------------------------------------------------------------------#
    # ---- All Filter functions below ---------------------------------------------------------------------#
    # -----------------------------------------------------------------------------------------------------#
    def update_score_filter(self):
        self.minScoreLine.setText(str(self.scoreSlider.value()))

    #allows user to save what is currently in the table
    def save_data(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,
                                      "Enter Text File Name", ".txt",
                                      "Text Document (*.txt)" )
        f = open(str(filename[0]), "w+")
        for genomes in self.AllData:
            f.write("***"+genomes + "\n")
            for items in self.AllData[genomes]:
                for i in items:
                    f.write(str(i) + '|')
                f.write(str(self.highlighted[items[1]]))
                f.write("\n")
        f.close()

    #open any saved .txt of previous tables opened
    def open_data(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File')
        self.AllData.clear()
        self.highlighted.clear()
        first = 1
        list1 = []
        s = ""
        self.comboBoxGene.clear()
        if (os.path.isfile(str(filename[0]))):
            f = open(str(filename[0]), "r+")
            for line in f:
                if(line.startswith("***")):
                    if(first == 0):
                        self.AllData[s] = list1
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
            self.AllData[s] = list1
            f.close()
            self.displayGeneData()

        else:
            # change to dialog box
            print('Could not open file')


class geneViewerSettings(QtWidgets.QDialog):
    def __init__(self):
        # Qt init stuff
        super(geneViewerSettings, self).__init__()
        uic.loadUi("geneViewerSettings.ui", self)
        self.setWindowTitle("Change Gene Viewer Settings")
        self.setWindowIcon(Qt.QIcon("cas9image.png"))

        # button connections
        self.kegg_radio_button.clicked.connect(self.change_file_type)
        self.gbff_radio_button.clicked.connect(self.change_file_type)
        self.fna_radio_button.clicked.connect(self.change_file_type)
        self.browse_button.clicked.connect(self.browseForFile)
        self.cancel_button.clicked.connect(self.cancelFunction)
        self.submit_button.clicked.connect(self.submitFunction)

        # class variables
        self.file_type = ""

    # this function is called when the user changes the file type
    # it just sets a class variable to the type of file selected
    def change_file_type(self):
        if self.kegg_radio_button.isChecked():
            self.file_type = "kegg"
        elif self.gbff_radio_button.isChecked():
            self.file_type = "gbff"
        elif self.fna_radio_button.isChecked():
            self.file_type = "fna"


    # this function is only called when the user selects browse for a file option
    # it opens a window such that the user can search for a file to use for the gene viewer sequence
    def browseForFile(self):
        # return out if the user has selected Kegg
        if self.kegg_radio_button.isChecked():
            return

        # make sure that either GBFF or FNA is checked
        if not self.gbff_radio_button.isChecked() and not self.fna_radio_button.isChecked():
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "Please select either the GBFF or the FNA radio button.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # open a window so that the user can select a file
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose an Annotation File")

        # make sure they choose the correct type of file
        if self.file_type not in myFile[0]:
            QtWidgets.QMessageBox.question(self, "Wrong type of file selected",
                                           "Please select the same type of file selected in the radio buttons.",
                                           QtWidgets.QMessageBox.Ok)
            self.file_name_edit.setText("")
            return

        # if the file is not empty, then set it
        if (myFile[0] != ""):
            self.file_name_edit.setText(myFile[0])

    # this function is only called when the user clicks on cancel.
    # it resets the text for the file chosen, while also unchecking all of the radio buttons
    # then it hides the window
    def cancelFunction(self):
        self.file_name_edit.setText("Please choose a file!")

        # setting them to False does not seem to work, no idea why
        self.kegg_radio_button.setChecked(False)
        self.gbff_radio_button.setChecked(False)
        self.fna_radio_button.setChecked(False)
        self.hide()

    # this function is only called when the user clicks on submit
    # it will find all of the sequences for all of the genes in the comboGeneBox in the results window
    # It will go based on the lengths stored in the comboGeneBox dictionary
    def submitFunction(self):
        if not self.kegg_radio_button.isChecked() and (self.file_name_edit.displayText() == "" or self.file_name_edit.displayText() == "Please choose a file!"):
            print("No file chosen")
            return

        sequence = ""

        # for each gene selected from the results window
        for item in GlobalSettings.mainWindow.Results.geneDict:
            if self.file_type == "kegg":
                k = Kegg()
                organism = GlobalSettings.mainWindow.Annotations_Organism.currentText().split(" ")[1]
                nameFull = item.split(" ")
                name = nameFull[len(nameFull) - 1]
                # get kegg's ntsequence and store it
                nt_sequence = k.get_nt_sequence(organism + ":" + name)
                GlobalSettings.mainWindow.Results.geneNTDict[item] = nt_sequence
            if self.file_type == "fna":
                sequence = self.fna_sequence_finder(GlobalSettings.mainWindow.Results.geneDict[item])
                GlobalSettings.mainWindow.Results.geneNTDict[item] = sequence
            if self.file_type == "gbff":
                sequence = self.gbff_sequence_finder(GlobalSettings.mainWindow.Results.geneDict[item])
                GlobalSettings.mainWindow.Results.geneNTDict[item] = sequence


        GlobalSettings.mainWindow.Results.displayGeneViewer.setEnabled(True)
        GlobalSettings.mainWindow.Results.lineEditStart.setEnabled(True)
        GlobalSettings.mainWindow.Results.lineEditEnd.setEnabled(True)
        GlobalSettings.mainWindow.Results.change_start_end_button.setEnabled(True)
        GlobalSettings.mainWindow.Results.checkGeneViewer()
        self.hide()

    # this function gets the sequence out of the GBFF file
    # may have indexing issues
    def gbff_sequence_finder(self, location_data):
        # start up the function
        fileStream = open(self.file_name_edit.displayText())
        buffer = fileStream.readline()
        print(buffer)
        index = 1
        pre_sequence = ""

        # skip all of the data until we are at the chromesome we care about
        while index != location_data[0]:
            if "ORIGIN" in buffer:
                index += 1
            buffer = fileStream.readline()

        # get the entire chromesome into a string
        while "//" not in buffer:
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

        # get the correct location and return
        ret_sequence = pre_sequence[location_data[1]:location_data[2]]
        return ret_sequence

    # this function is the function that actually finds the sequence
    # May have indexing issues here
    def fna_sequence_finder(self, location_data):
        # Open the file and set the index to 0
        fileStream = open(self.file_name_edit.displayText())
        index = 1

        # skip the file until we get to the chromesome we want
        buffer = fileStream.readline()
        while index != location_data[0]:
            buffer = fileStream.readline()
            if buffer.startswith(">"):
                index += 1

        buffer = fileStream.readline()

        # now go through and get the whole chromesome
        sequence = buffer
        while not buffer.startswith(">"):
            buffer = fileStream.readline()
            if not buffer.startswith(">"):
                sequence = sequence + buffer

        # uppercase that chromesome, and change all endlines with spaces
        sequence = sequence.upper()
        sequence = sequence.replace("\n", "")

        # now set the return sequence to the substring that we want
        NTSequence = sequence[location_data[1]:location_data[2]]

        return NTSequence






# Window opening and GUI launching code for debugging #
# ----------------------------------------------------------------------------------------------------- #
"""
app = Qt.QApplication(sys.argv)
app.setOrganizationName("TrinhLab-UTK")
app.setApplicationName("CASPER")
window = Results()
window.transfer_data("yli", "spCas9", "/Users/brianmendoza/Dropbox/CrisprDB/", {"phos.carboxylase":(2,3030460,3032157)}, "/Volumes/Seagate_Drive/FASTAs/yli.fna")
sys.exit(app.exec_())"""