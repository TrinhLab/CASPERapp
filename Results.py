import sys

from PyQt5 import Qt, QtWidgets, uic
from Scoring import OnTargetScore
from Algorithms import SeqTranslate

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

        # Main Data container
        # Keys: Gene names
        # Values: #
        self.AllData = {}


        self.startpos = 0
        self.endpos = 0
        self.directory = ""

        # Target Table settings #
        self.targetTable.setColumnCount(7)  # hardcoded because there will always be seven columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Location;Sequence;Strand;PAM;Score;Off-Target;Highlight".split(";"))
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        self.show()

    # Function that is called in main in order to pass along the information the user inputted and the information
    # from the .cspr files that was discovered
    def transfer_data(self, org, endo, path, geneposdict, fasta):
        self.org = org
        self.endo = endo
        self.directory = path
        self.fasta_ref = fasta
        for gene in geneposdict:
            self.comboBoxGene.addItem(gene)
            self.get_targets(gene, geneposdict[gene])

    # Function grabs the information from the .cspr file and adds them to the AllData dictionary
    def get_targets(self, genename, pos_tuple):
        targets = []
        if self.directory.find("/") != -1:
            file = open(self.directory+"/" + self.org + "_" + self.endo + ".cspr")
        else:
            file = open(self.directory + "\\" + self.org + "_" + self.endo + ".cspr")
        header = file.readline()
        # Find the right chromosome:
        while True:
            # in the right chromosome/scaffold?
            if header.find("(" + str(pos_tuple[0]) + ")"):
                while True:
                    # Find the appropriate location by quickly decompressing the location at the front of the line
                    myline = file.readline()
                    if self.S.decompress64(myline.split(",")[0]) >= pos_tuple[1]:
                        while self.S.decompress64(myline.split(",")[0]) < pos_tuple[2]:
                            targets.append(self.S.decompress_csf_tuple(myline))
                            myline = file.readline()
                    else:
                        continue
                    break
                break
            else:
                header = file.readline()

        self.AllData[genename] = targets
        self.displayGeneData()

    def displayGeneData(self):
        curgene = str(self.comboBoxGene.currentText())
        #cg = self.allGeneSeqs[curgene]
        #self.geneViewer.setPlainText(cg)
        #  --- Shifting numbers over based on start and end ---  #

        self.targetTable.setRowCount(len(self.AllData[curgene]))
        print(self.AllData[curgene])
        index = 0
        for item in self.AllData[curgene]:
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
            """st = item[0]-self.startpos
            print(st)
            seq = QtWidgets.QTableWidgetItem(cg[st-20:st])
            self.targetTable.setItem(index, 0, seq)
            pam = QtWidgets.QTableWidgetItem(cg[st:st+3])  # this is only for a 3nt PAM need to import pam info
            self.targetTable.setItem(index, 1, pam)
            strand = QtWidgets.QTableWidgetItem(item[1])
            self.targetTable.setItem(index, 2, strand)
            scr = self.onscore.returnScore(cg[st-6:st+30])
            print(scr)
            score = QtWidgets.QTableWidgetItem(str(scr))
            self.targetTable.setItem(index, 6, score)"""
            self.btn_sell = QtWidgets.QPushButton('Find Off Targets')
            self.btn_sell.clicked.connect(self.handleButtonClicked)
            self.targetTable.setCellWidget(index, 5, self.btn_sell)
            ckbox = QtWidgets.QCheckBox()
            ckbox.clicked.connect(self.search_gene)
            self.targetTable.setCellWidget(index,6,ckbox)
            """x = self.targetTable.cellWidget(index,6)
            y = x.text()
            x.setTextAlignment(Qt.AlignHCenter())"""
            index += 1
        self.targetTable.resizeColumnsToContents()

    def search_gene(self):
        search_trms = []
        checkBox = self.sender()
        index = self.targetTable.indexAt(checkBox.pos())
        print(index.column(), index.row(), checkBox.isChecked())
        #self.geneViewer.setPlainText(self.geneViewer.text())
        seq = self.targetTable.item(index.row(),1).text()
        print(seq)
        x=1



    def handleButtonClicked(self):
        # button = QtGui.qApp.focusWidget()
        button = self.sender()
        index = self.targetTable.indexAt(button.pos())
        if index.isValid():
            print(index.row(), index.column())


# Window opening and GUI launching code #
# ----------------------------------------------------------------------------------------------------- #
app = Qt.QApplication(sys.argv)
app.setOrganizationName("TrinhLab-UTK")
app.setApplicationName("CASPER")
window = Results()
window.transfer_data("yli", "spCas9", "/Users/brianmendoza/Dropbox/CrisprDB", {"myfakegene":(1,1293,3496)}, "/Volumes/Seagate_Drive/FASTAs/yli.fna")
sys.exit(app.exec_())