import models.GlobalSettings as GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic
import traceback
import math
from utils.ui import show_error, scale_ui

logger = GlobalSettings.logger

###########################################################
# closingWindow: this class is a little window where the user can select which files they want to delete
# Once they hit 'submit' it will delete all of the files selected, and close the program.
# If no files are selected, the program closes and no files are deleted
# Inputs are taking from the user (selecting files to delete and hitting submit), as well as GlobalSettings for the files in CSPR_DB
# Outputs are the files are deleting, and the program is closed
###########################################################
class closingWindow(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(closingWindow, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "ui/closing_window.ui", self)
            self.setWindowTitle("Delete Files")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))

            # Button 
            self.submit_button.clicked.connect(self.submit_and_close)

            # Table 
            self.files_table.setColumnCount(1)
            self.files_table.setShowGrid(True)
            self.files_table.setHorizontalHeaderLabels("File Name;".split(";"))
            self.files_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.files_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.files_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)


            scale_ui(self, custom_scale_width=400, custom_scale_height=300)


        except Exception as e:
            show_error("Error initializing closingWindow class.", e)

    # this function will delete selected files, and then close the program
    def submit_and_close(self):
        try:
            # loop through the whole table
            for i in range(self.files_table.rowCount()):
                tabWidget = self.files_table.item(i, 0)

                # if that specific tab is selected, delete it. otherwise do nothing
                if tabWidget.isSelected():
                    os.remove(tabWidget.text())
            self.close()
        except Exception as e:
            show_error("Error in sumbit_and_close() in closing window.", e)
            
    # this function gets all of the files from the CSPR_DB and puts them all into the table
    def get_files(self):
        try:
            loopCount = 0
            # get the file names from CSPR_DB
            files_names = os.listdir(GlobalSettings.CSPR_DB)
            files_names.sort(key=str.lower)
            self.files_table.setRowCount(len(files_names))

            # loop through and add them to the table
            for file in files_names:
                tabWidget = QtWidgets.QTableWidgetItem(file)
                self.files_table.setItem(loopCount, 0, tabWidget)
                loopCount += 1
            self.files_table.resizeColumnsToContents()
        except Exception as e:
            show_error("Error in get_files() in closing window.", e)