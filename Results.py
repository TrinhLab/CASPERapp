import sys

from PyQt5 import Qt, QtWidgets, uic
from Scoring import OnTargetScore

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
        # Data containers #
        self.allGenes = []
        self.allTargets = {}
        self.allGeneSeqs = {}

        self.startpos = 0
        self.endpos = 0
        self.directory ='C:/Users/GregCantrall/Documents/Cspr files'



        # Target Table settings #
        self.targetTable.setColumnCount(7)  # hardcoded because there will always be seven columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Location;Sequence;Strand;PAM;Score;Off-Target;Highlight".split(";"))
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.targetTable.itemSelectionChanged.connect(self.selected)

        self.fill_table_TEST()
        self.show()


    def selected(self):
        for item in self.targetTable.selectedItems():
            print(item.text())
    def getTargets(self, fileName):
        file = ""
        targets = []
        if self.directory.find("/") != -1:
            file = open(self.directory+"/"+fileName+".cspr")
        else:
            file = open(self.directory + "\\" + fileName + ".cspr")
        hold = file.readline()
        hold  = hold[8:hold.find('\n')]
        file.readline()
        while True:
            string = file.readline()
            if "REPEATS" in string:
                break
            item = self.splitCsprFile(string)
            targets.append(item)
        self.allTargets[hold] = targets
        self.displayGeneData()



    def splitCsprFile(self, holder):
        item = []
        strand = "+"
        sep = holder.find(',')
        item.append(holder[:sep])
        holder = holder[sep+1:]
        sep = holder.find('+')
        if sep == -1:
            strand = "-"
            sep = holder.find("-")
        item.append(holder[:sep])
        item.append(strand)
        holder = holder[sep+1:]
        sep = holder.find(",")
        item.append(holder[:sep])
        item.append(holder[sep+1:])
        return item


    def displayGeneData(self):
        #curgene = str(self.comboBoxGene.currentText())
        curgene = 'Thermoanaerobacterium saccharolyticum'
        #cg = self.allGeneSeqs[curgene]
        #self.geneViewer.setPlainText(cg)
        #  --- Shifting numbers over based on start and end ---  #

        self.targetTable.setRowCount(len(self.allTargets[curgene]))
        print(self.allTargets[curgene])
        index = 0
        for item in self.allTargets[curgene]:
            loc = QtWidgets.QTableWidgetItem(item[0])
            seq = QtWidgets.QTableWidgetItem(item[1])
            strand = QtWidgets.QTableWidgetItem(item[2])
            PAM = QtWidgets.QTableWidgetItem(item[3])
            score = QtWidgets.QTableWidgetItem(item[4])
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

    #-----Testing Methods -----#
    def fill_table_TEST(self):
        y=2
        x = self.getTargets("tsh_spCas9-VRER")

        #self.loadGenesandTargets("testing_seq1",1,3,["target1","target2","target3"],"testo")
        #self.splitCsprFile("BY,Xc9d+CV,q")
        """self.targetTable.setRowCount(3)
        seq = QtWidgets.QTableWidgetItem("testing")
        self.targetTable.setItem(0, 0,seq )
        seq = QtWidgets.QTableWidgetItem("other")
        self.targetTable.setItem(1, 0, seq)
        seq = QtWidgets.QTableWidgetItem("third")
        self.targetTable.setItem(2, 0, seq)
        self.geneViewer.setPlainText("this is testing the third other thing")
        self.geneViewer.setFontItalic(True)
        self.geneViewer.find("testing")"""





# ----------------------------------------------------------------------------------------------------- #
app = Qt.QApplication(sys.argv)
app.setOrganizationName("TrinhLab-UTK")
app.setApplicationName("CASPER")
window = Results()
sys.exit(app.exec_())