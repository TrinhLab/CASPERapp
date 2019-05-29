import sys

from PyQt5 import Qt, QtWidgets, uic, QtCore
from Scoring import OnTargetScore
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
import GlobalSettings
import os
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

        self.switcher = [1,1,1,1,1,1,1]  # for keeping track of where we are in the sorting clicking for each column

        # Target Table settings #
        self.targetTable.setColumnCount(7)  # hardcoded because there will always be seven columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Location;Sequence;Strand;PAM;Score;Off-Target;Highlight".split(";"))
        self.targetTable.horizontalHeader().setSectionsClickable(True)
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        self.back_button.clicked.connect(self.goBack)
        self.targetTable.horizontalHeader().sectionClicked.connect(self.table_sorting)
        self.actionSave.triggered.connect(self.save_data)
        self.actionOpen.triggered.connect(self.open_data)
        self.changeEndoButton.clicked.connect(self.changeEndonuclease)

        self.targetTable.itemSelectionChanged.connect(self.item_select)
        self.minScoreLine.setText("0")

        # Connecting the filters to the displayGeneData function
        self.fivegseqCheckBox.stateChanged.connect(self.displayGeneData)
        self.minScoreLine.textChanged.connect(self.displayGeneData)

        # Setting up the score filter:
        self.scoreSlider.setMinimum(0)
        self.scoreSlider.setMaximum(100)
        self.scoreSlider.setTracking(False)

        self.scoreSlider.valueChanged.connect(self.update_score_filter)



    # this function goes through and calls transfer_data again.
    # Uses data from the mainWindow in Globalsettings, but that's because that info should not change
    # unless the user closes out of the Results window
    def changeEndonuclease(self):
        full_org = str(GlobalSettings.mainWindow.orgChoice.currentText())
        organism = GlobalSettings.mainWindow.shortHand[full_org]
        self.transfer_data(organism, str(self.endonucleaseBox.currentText()), GlobalSettings.CSPR_DB, GlobalSettings.mainWindow.checked_info, "")

    # Function that is called in main in order to pass along the information the user inputted and the information
    # from the .cspr files that was discovered
    def transfer_data(self, org, endo, path, geneposdict, fasta):
        self.org = org
        self.endo = endo
        self.directory = path
        self.fasta_ref = fasta
        self.comboBoxGene.clear()
        self.AllData.clear()

        # print("Inside transfer_data, here is the stuff")
        # print("Org: ", org)
        # print("Endo: ", endo)
        # print("Path: ", path)
        # print("GenePosDict: ", geneposdict)
        # print("Fasta: ", fasta)

        self.highlighted.clear()
        for gene in geneposdict:
            self.comboBoxGene.addItem(gene)
            # print('gene: ' + gene)
            self.get_targets(gene, geneposdict[gene])

        # Enable the combobox to be toggled now that the data is in AllData
        self.comboBoxGene.currentTextChanged.connect(self.displayGeneData)
    def goBack(self):

        GlobalSettings.mainWindow.show()
        self.hide()



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
        print('postupe: '+str(pos_tuple))
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
        #temp = list(self.AllData.keys())
        #curgene = temp[0]
        #print('curgene: ' + curgene)
        # Creates the set object from the list of the current gene:
        if curgene=='' or len(self.AllData)<1:
            return

        subset_display = set()
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
        for item in subset_display:
            loc = QtWidgets.QTableWidgetItem(str(item[0]))
            seq = QtWidgets.QTableWidgetItem(item[1])
            strand = QtWidgets.QTableWidgetItem(str(item[4]))
            PAM = QtWidgets.QTableWidgetItem(item[2])
            score = QtWidgets.QTableWidgetItem(str(item[3]))
            self.targetTable.setItem(index, 0, loc)
            self.targetTable.setItem(index, 1, seq)
            self.targetTable.setItem(index, 2, strand)
            self.targetTable.setItem(index, 3, PAM)
            self.targetTable.setItem(index, 4, score)
            self.targetTable.setItem(index, 5, QtWidgets.QTableWidgetItem("--.--"))
            ckbox = QtWidgets.QCheckBox()
            ckbox.clicked.connect(self.search_gene)
            self.targetTable.setCellWidget(index,6,ckbox)
            #print(self.highlighted)
            if(self.highlighted[str(item[1])] == True):
                ckbox.click()
            index += 1
        self.targetTable.resizeColumnsToContents()

    ########################################## END UPDATING FUNCTION #############################################

    def search_gene(self):
        search_trms = []
        checkBox = self.sender()
        index = self.targetTable.indexAt(checkBox.pos())
        print(index.column(), index.row(), checkBox.isChecked())
        #self.geneViewer.setPlainText(self.geneViewer.text())
        seq = self.targetTable.item(index.row(),1).text()
        self.highlighted[str(seq)] = checkBox.isChecked()
        print(seq)
        x=1

    def item_select(self):
        print(self.targetTable.selectedItems())

    def table_sorting(self, logicalIndex):
        self.switcher[logicalIndex] *= -1
        if self.switcher[logicalIndex] == -1:
            self.targetTable.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.targetTable.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)


    # currently unused, but could be used with the offTargetButton when implemented
    def offtargetButtonClicked(self):
        all_targets = []
        index = 0
        holdt = ()
        QtWidgets.QMessageBox.question(self, "Text File Created",
                                       "A File with the off target data has been created.",
                                       QtWidgets.QMessageBox.Ok)

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


    def save_data(self):
        #filed = QtWidgets.QFileDialog()
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File')
        #print(filename[0])
        if(os.path.isfile(str(filename[0]))):
            #print('success')
            f = open(str(filename[0]), "w+")
            for genomes in self.AllData:
                f.write("***"+genomes + "\n")
                for items in self.AllData[genomes]:
                    for i in items:
                        f.write(str(i) + '|')
                    f.write(str(self.highlighted[items[1]]))
                    f.write("\n")
            f.close()
        else:
            #change to dialog box
            print('Could not open file')


    def open_data(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File')
        #print(filename[0])
        #print(self.AllData)
        print(self.highlighted)
        self.AllData.clear()
        self.highlighted.clear()
        #print(self.AllData)
        first = 1
        list1 = []
        s = ""
        self.comboBoxGene.clear()
        if (os.path.isfile(str(filename[0]))):
            print('success')
            f = open(str(filename[0]), "r+")
            for line in f:
                if(line.startswith("***")):
                    if(first == 0):
                        self.AllData[s] = list1
                        #print(list1)
                    else:
                        first = 0
                    s = line[3:]
                    s = s.strip('\n')
                    list1 = []
                    self.comboBoxGene.addItem(s)
                    print('s: '+s)
                else:
                    temp = line.split("|")
                    h = temp.pop()
                    h = h.strip('\n')
                    self.highlighted[temp[1]] = eval(h)
                    print(temp[1])
                    print('h: '+str(h))


                    tup = tuple(temp)
                    list1.append(tup)
            self.AllData[s] = list1
            f.close()
            #print(self.AllData)
            print(self.highlighted)
            self.displayGeneData()

        else:
            # change to dialog box
            print('Could not open file')



# Window opening and GUI launching code for debugging #
# ----------------------------------------------------------------------------------------------------- #
"""
app = Qt.QApplication(sys.argv)
app.setOrganizationName("TrinhLab-UTK")
app.setApplicationName("CASPER")
window = Results()
window.transfer_data("yli", "spCas9", "/Users/brianmendoza/Dropbox/CrisprDB/", {"phos.carboxylase":(2,3030460,3032157)}, "/Volumes/Seagate_Drive/FASTAs/yli.fna")
sys.exit(app.exec_())"""