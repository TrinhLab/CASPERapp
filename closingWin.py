import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic
import traceback
import math

#global logger
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
            # qt stuff
            super(closingWindow, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "closing_window.ui", self)
            self.setWindowTitle("Delete Files")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))

            # button connections
            self.submit_button.clicked.connect(self.submit_and_close)

            # table stuff
            self.files_table.setColumnCount(1)
            self.files_table.setShowGrid(True)
            self.files_table.setHorizontalHeaderLabels("File Name;".split(";"))
            self.files_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.files_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.files_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing closingWindow class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.fontSize = fontSize
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 1150x650
            scaledWidth = int((width * 400) / 1920)
            scaledHeight = int((height * 300) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(scaledWidth / 2))
            y = y - (math.ceil(scaledHeight / 2))
            self.setGeometry(x, y, scaledWidth, scaledHeight)

            self.repaint()
            QtWidgets.QApplication.processEvents()

        except Exception as e:
            logger.critical("Error in scaleUI() in closing window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            #center window on current screen
            width = self.width()
            height = self.height()
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in centerUI() in closing window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # this function will delete selected files, and then close the program
    def submit_and_close(self):
        try:
            # loop through the whole table
            for i in range(self.files_table.rowCount()):
                tabWidget = self.files_table.item(i, 0)

                # if that specific tab is selected, delete it. otherwise do nothing
                if tabWidget.isSelected():
                    os.remove(tabWidget.text())

            # close the program now
            self.close()
        except Exception as e:
            logger.critical("Error in sumbit_and_close() in closing window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

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
            logger.critical("Error in get_files() in closing window.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)