from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic

class HoverButton(QtGui.QToolButton):
    def __init__(self, parent=none):
        super(HoverButton, self).__init__(parent)
        self.setMouseTracking(True)

    def enterEvent(self, event):
        print("Enter")
        self.setStyleSheet("background-color:#45b545;")

    def leaveEvent(self, event):
        self.setStyleSheet("background-color:yellow;")
        print("Leave")