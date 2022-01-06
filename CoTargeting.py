from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import traceback
import math

#global logger
logger = GlobalSettings.logger

######################################################
# CoTargeting class: This class is a window that lets the user select which endonucleases to co-target with
# inputs are from the user and from results
# from results: the organism name and the list of endonucleases for that organism
# from user: which endonucleases to co-target
######################################################
class CoTargeting(QtWidgets.QMainWindow):

    def __init__(self, path):
        try:
            # pyqt stuff
            super(CoTargeting, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'cotargeting.ui', self)
            self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Co-targeting")

            # endo_table stuff
            self.endo_table.setColumnCount(1)  # hardcoded because there will always be 1 columns
            self.endo_table.setShowGrid(True)
            self.endo_table.setHorizontalHeaderLabels("Endonuclease;".split(";"))
            self.endo_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.endo_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.endo_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.endo_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch) #Ensures last column goes to the edge of table

            # variables
            self.info_path = path

            # button connections
            self.cancel_button.clicked.connect(self.cancel_function)
            self.submit_button.clicked.connect(self.submission_function)

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing CoTargeting class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

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

            # CASPER header scaling
            fontSize = 20
            self.title.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 850x750
            scaledWidth = int((width * 450) / 1920)
            scaledHeight = int((height * 375) / 1080)

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
            logger.critical("Error in scaleUI() in cotargeting in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            # center window on current screen
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
            logger.critical("Error in centerUI() in cotargeting in results.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # launches the window
    # it is expecting endo_choices in the form of a list, and the orgName in the form of a string
    # it sets the organism name, and sets the table as well
    def launch(self, endo_choices, orgName):
        try:
            self.orgName.setText(orgName)
            setTableList = list()
            # only get the endo choices that were original
            for item in endo_choices:
                checkList = item.split(",")
                if len(checkList) == 1 and "|" not in item: #Prevent cotarget endos from being added back in
                    setTableList.append(item)

            # go through and set each table item, but also set the row count
            self.endo_table.setRowCount(len(setTableList))
            loopCount = 0
            for item in setTableList:
                tabWidget = QtWidgets.QTableWidgetItem(item)
                self.endo_table.setItem(loopCount, 0, tabWidget)
                loopCount += 1
            self.endo_table.resizeColumnsToContents()

            # now show
            self.centerUI()
            self.show()
            self.activateWindow()
        except Exception as e:
            logger.critical("Error in launch() in CoTargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # this is the cancel function
    # it clears the table, sets the organism to nothing, and hides the window
    def cancel_function(self):
        try:
            self.endo_table.clearContents()
            self.endo_table.setRowCount(0)
            self.orgName.setText("")
            self.hide()
        except Exception as e:
            logger.critical("Error in cancel_function() in CoTargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # this is the submission function.
    # it makes sure the user selects at least 2 endonucleases
    # then it goes through and returns the endonucleases selected
    # once it gets those, it then calls a function in Results that repopulates the table correctly
    def submission_function(self):
        try:
            #get endo data from CASPERinfo
            self.Endos = {}
            f = open(GlobalSettings.appdir + 'CASPERinfo')
            while True:
                line = f.readline()
                if line.startswith('ENDONUCLEASES'):
                    while True:
                        line = f.readline()
                        if (line[0] == "-"):
                            break
                        line_tokened = line.split(";")
                        endo = line_tokened[0]
                        self.Endos[endo] = ([line_tokened[2], line_tokened[3], line_tokened[4]],line_tokened[5])
                    break
            f.close()

            # set the selected_list, and make sure they select at least 2 endonucleases
            selected_list = self.endo_table.selectedItems()
            if len(selected_list) <= 1:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Nothing Selected")
                msgBox.setText(
                    "No endonucleases selected. Please select at least 2 endonucleases")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            # go through and get which endonuclease's have been selected
            ret_endo_list = list()
            for i in range(self.endo_table.rowCount()):
                if self.endo_table.item(i, 0).isSelected():
                    ret_endo_list.append(self.endo_table.item(i, 0).text())

            #invalid_flag = False
            for endo1 in ret_endo_list:
                for endo2 in ret_endo_list:
                    if endo1 == endo2:
                        continue
                    endo1_len = sum([int(x) for x in self.Endos[endo1][0]])
                    endo2_len = sum([int(x) for x in self.Endos[endo2][0]])
                    if endo1_len != endo2_len or self.Endos[endo1][1] != self.Endos[endo2][1]: # If endonucleases don't have the same length gRNA or don't have the same directionality, throw an error
                        msgBox = QtWidgets.QMessageBox()
                        msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                        msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                        msgBox.setWindowTitle("Invalid Endonucleases")
                        msgBox.setText(
                            "The selected endonucleases are not compatible.")
                        msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                        msgBox.exec()

                        return

            GlobalSettings.mainWindow.Results.co_target_endo_list = ret_endo_list
            GlobalSettings.mainWindow.Results.populate_cotarget_table()
            self.cancel_function()
        except Exception as e:
            logger.critical("Error in submission_function() in CoTargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

