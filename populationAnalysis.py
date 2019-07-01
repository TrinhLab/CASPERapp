from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings

class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi('populationanalysis.ui', self)
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.goBackButton.clicked.connect(self.go_back)

    def go_back(self):
        GlobalSettings.mainWindow.show()
        self.hide()