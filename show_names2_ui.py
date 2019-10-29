import os
from PyQt5 import QtWidgets, uic, QtCore
from functools import partial
import GlobalSettings

class show_names_table2(QtWidgets.QDialog):
    def __init__(self):
        super(show_names_table2, self).__init__()
        uic.loadUi('name_form2.ui', self)
        self.name_table2.setColumnCount(1)

    def fill_table(self, names):
        index = 0
        labels = []
        self.name_table2.setRowCount(0)
        for name in names:
            self.name_table2.setRowCount(index+1)
            n = QtWidgets.QTableWidgetItem()
            n.setData(QtCore.Qt.EditRole, str(name))
            self.name_table2.setItem(index, 0, n)
            labels.append(str(index))
            index += 1
        self.name_table2.resizeColumnsToContents()
        self.name_table2.setVerticalHeaderLabels(labels)


    def closeEvent(self, event):
        self.hide()
        event.accept()
