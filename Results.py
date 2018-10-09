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
        self.targetTable.setColumnCount(6)  # hardcoded because there will always be five columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Sequence;PAM;Strand;Score;Off Targets;search".split(";"))
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        self.fill_table_TEST()
        #self.show()

    def loadGenesandTargets(self, genesequence, start, end, targets, genename):
        self.startpos = start
        self.endpos = end
        self.allTargets[genename] = targets
        self.allGeneSeqs[genename] = genesequence
        self.comboBoxGene.addItem(genename)
        self.displayGeneData()
        self.comboBoxGene.currentIndexChanged.connect(self.displayGeneData)


    def getTargets(self, fileName):
        file = ""
        targets = []
        if self.directory.find("/") != -1:
            file  = open(self.directory+"/"+fileName+".cspr")
        else:
            file = open(self.directory + "/" + fileName + ".cspr")
        file.readline()
        file.readline()
        for string in file.readline():
            item  = self.splitCsprFile(string)
            targets.append(item)
        return targets


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
        curgene = str(self.comboBoxGene.currentText())
        cg = self.allGeneSeqs[curgene]
        self.geneViewer.setPlainText(cg)
        #  --- Shifting numbers over based on start and end ---  #

        self.targetTable.setRowCount(len(self.allTargets[curgene]))
        print(self.allTargets[curgene])
        index = 0
        for item in self.allTargets[curgene]:
            st = item[0]-self.startpos
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
            self.targetTable.setItem(index, 3, score)
            self.btn_sell = QtWidgets.QPushButton('Find Off Targets')
            self.btn_sell.clicked.connect(self.handleButtonClicked)
            self.targetTable.setCellWidget(index, 4, self.btn_sell)
            index += 1
            ckbox = QtWidgets.QCheckBox()
            ckbox.clicked.connect(self.search_gene)
            self.targetTable.setCellWidget(index,5,ckbox)

        self.targetTable.resizeColumnsToContents()

    def search_gene(self):
        search_trms = []
        for item in self.targetTable:
            if item[5].isChecked():
                search_trms.append(item[0])
        self.geneViewer.setPlainText(self.geneViewer.text())
        for item in search_trms:
            x=0





    def handleButtonClicked(self):
        # button = QtGui.qApp.focusWidget()
        button = self.sender()
        index = self.targetTable.indexAt(button.pos())
        if index.isValid():
            print(index.row(), index.column())

    #-----Testing Methods -----#
    def fill_table_TEST(self):
        #self.getTargets("")
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
"""app = Qt.QApplication(sys.argv)
app.setOrganizationName("TrinhLab-UTK")
app.setApplicationName("CASPER")
window = Results()
sys.exit(app.exec_())"""