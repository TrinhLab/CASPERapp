import sys
from PyQt5 import QtGui, QtWidgets, QtCore, Qt

my_array = [['AAGGGTACCATCACTGG', 'AGG', 'Reverse', '55'],
                    ['GAGGATCCATAACATAGAT', 'GGG', 'Forward', '60'],
                    ['AGGGAGAGGAAGAACGT', 'TGG', 'Forward', '27']]

my_array_2 = [['TACCATAGGATACCAGT', 'GGG', 'Forward', '44'],
                    ['AGAAGAGATTCCCAGATG', 'CGG', 'Forward', '39'],
                    ['GGGATACGACGATAGCG', 'CGG', 'Forward', '54']]

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(5)
        self.setCentralWidget(self.table)
        self.data = []

        self.table.setRowCount(len(my_array))
        self.table.setShowGrid(False)
        self.table.setHorizontalHeaderLabels("Sequence;PAM;Strand;Score;Off Targets".split(";"))
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)


        """self.myLabel = QtGui.QLabel()
        mystring = ""
        f = open('CASPERinfo.txt')
        for line in f:
            mystring += line
        f.close()
        self.myLabel.setText(mystring)
        self.myLabel.show()"""

        index = 0
        for item in my_array:
            seq = QtWidgets.QTableWidgetItem(item[0])
            self.table.setItem(index,0,seq)
            pam = QtWidgets.QTableWidgetItem(item[1])
            self.table.setItem(index,1, pam)
            strand = QtWidgets.QTableWidgetItem(item[2])
            self.table.setItem(index, 2, strand)
            score = QtWidgets.QTableWidgetItem(item[3])
            self.table.setItem(index, 3, score)
            self.btn_sell = QtWidgets.QPushButton('Find Off Targets')
            self.btn_sell.clicked.connect(self.handleButtonClicked)
            self.table.setCellWidget(index,4,self.btn_sell)
            index += 1

        self.table.resizeColumnsToContents()

    def handleButtonClicked(self):
        # button = QtGui.qApp.focusWidget()
        button = self.sender()
        index = self.table.indexAt(button.pos())
        if index.isValid():
            print(index.row(), index.column())

    def importData(self, data):
        self.data = data

def main(args):
    app = Qt.QApplication(args)
    table = MainWindow()
    table.importData(my_array)
    table.show()
    sys.exit(app.exec_())

if __name__=="__main__":
    main(sys.argv)