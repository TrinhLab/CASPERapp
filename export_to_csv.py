import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore
import platform
import traceback
import math

#global logger
logger = GlobalSettings.logger

# Class: export_csv_window
# This class opens a window for the user to select where they want the CSV file exported to, and the name of the file
# It takes the highlighted data from the Results page, and creates a CSV file from that
class export_csv_window(QtWidgets.QMainWindow):
    # init function. Sets all of the buttons
    def __init__(self):
        try:
            # qt stuff
            super(export_csv_window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'export_to_csv_window.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))

            # button connections
            self.browse_button.clicked.connect(self.browseForFolder)
            self.cancel_button.clicked.connect(self.cancel_function)
            self.export_button.clicked.connect(self.export_function)

            # variables
            self.location = self.fileLocation_line_edit.text()
            self.selected_table_items = []
            self.num_columns = []

            self.setWindowTitle("Export to CSV")
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing export_to_csv_window class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            screen = QtWidgets.QApplication.screens()[screen]
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
            scaledWidth = int((width * 650) / 1920)
            scaledHeight = int((height * 200) / 1080)
            
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
            logger.critical("Error in scaleUI() in export to csv.")
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
            logger.critical("Error in centerUI() in export to csv.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # launch function. Called in Results.
    # parameter expect: a list of the items selected from the window.
    def launch(self, select_items, num_columns):
        try:
            if platform.system() == "Windows":
                self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "\\")
            else:
                self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "/")
            self.selected_table_items = select_items
            self.num_columns = num_columns
            self.centerUI()
            self.show()
            self.activateWindow()
        except Exception as e:
            logger.critical("Error in launch() in export to csv.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # export function
    # Takes the path and file name and combines them
    # Writes the header line, as well as ever line selected to that file
    # calls the cancel function when it's done
    def export_function(self):
        try:
            # get the full path ( path and file name)
            file_name = self.filename_line_edit.text()
            if file_name == "":
                file_name = "exported_gRNAs"
            self.location = self.fileLocation_line_edit.text() + '/'
            full_path = ""
            if '.csv' in file_name:
                full_path = self.location + file_name
            else:
                full_path = self.location + file_name + '.csv'

            # try to do it
            try:
                #write the table headers
                if self.num_columns == 8: ###Change headers for multitargeting table export
                    output_data = open(full_path, 'w')
                    output_data.write('Seed,Total Repeats,Avg. Repeats/Scaffold,Consensus Sequence,% Consensus,Score,PAM,Strand\n')
                elif self.num_columns == 9:
                    output_data = open(full_path, 'w')
                    output_data.write('Seed,% Coverage,Total Repeats,Avg. Repeats/Scaffold,Consensus Sequence,% Consensus,Score,PAM,Strand\n')
                else: ###Change headers for view results export
                    output_data = open(full_path, 'w')
                    output_data.write('Location,Endonuclease,Sequence,Strand,PAM,Score,Off_Target\n')

                # loop through and write the other data
                i = 0
                for item in self.selected_table_items:
                    if i == self.num_columns-1:
                        output_data.write(item.text())
                        output_data.write('\n')
                        i = 0
                    else:
                        output_data.write(item.text())
                        output_data.write(',')
                        i += 1
            # catch the permission exception
            except PermissionError:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("File Cannot Open")
                msgBox.setText("This file cannot be opened. Please make sure that the file is not opened elsewhere and try again.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()


                return
            # catch any other exception
            except Exception as e:
                print(e)
                return

            # close the window
            self.cancel_function()
        except Exception as e:
            logger.critical("Error in export_function() in export to csv.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # Resets everything to the init funciton
    # then closes the window
    def cancel_function(self):
        try:
            if platform.system() == "Windows":
                self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "\\")
            else:
                self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "/")
            self.filename_line_edit.setText("")
            self.location = ""
            self.hide()
        except Exception as e:
            logger.critical("Error in cancel_function() in export to csv.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    # browse for folder function
    # allows user to browse for a folder where to store the CSV file
    def browseForFolder(self):
        try:
            # get the folder
            filed = QtWidgets.QFileDialog()
            mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                           GlobalSettings.CSPR_DB, QtWidgets.QFileDialog.ShowDirsOnly)
            if(os.path.isdir(mydir) == False):
                return


            if platform.system() == "Windows":
                self.fileLocation_line_edit.setText(mydir + "\\")
                self.location = mydir + "\\"
            else:
                self.fileLocation_line_edit.setText(mydir + "/")
                self.location = mydir + "/"
        except Exception as e:
            logger.critical("Error in browseForFolder() in export to csv.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)