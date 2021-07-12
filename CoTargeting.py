
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import traceback

#global logger
logger = GlobalSettings.logger

######################################################
# CoTargeting class: This class is a window that lets the user select which endonucleases to co-target with
# inputs are from the user and from results
# from results: the organism name and the list of endonucleases for that organism
# from user: which endonucleases to co-target
######################################################
class CoTargeting(QtWidgets.QDialog):

    def __init__(self, path):
        try:
            # pyqt stuff
            super(CoTargeting, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'co_targeting.ui', self)
            self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))

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

            #set pixel width for scroll bars
            self.endo_table.verticalScrollBar().setStyleSheet("width: 16px;")
            self.endo_table.horizontalScrollBar().setStyleSheet("height: 16px;")
        except Exception as e:
            logger.critical("Error initializing CoTargeting class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
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
                if len(checkList) == 1:
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
            self.show()
        except Exception as e:
            logger.critical("Error in launch() in CoTargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
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
                        self.Endos[endo] = [line_tokened[2], line_tokened[3], line_tokened[4]]
                    break
            f.close()

            # set the selected_list, and make sure they select at least 2 endonucleases
            selected_list = self.endo_table.selectedItems()
            if len(selected_list) <= 1:
                QtWidgets.QMessageBox.question(self, "Nothing Selected", "No endonucleases selected. Please select at least 2 endonucleases",
                                               QtWidgets.QMessageBox.Ok)
                return

            # go through and get which endonuclease's have been selected
            ret_endo_list = list()
            for i in range(self.endo_table.rowCount()):
                if self.endo_table.item(i, 0).isSelected():
                    ret_endo_list.append(self.endo_table.item(i, 0).text())

            #invalid_flag = False
            for endo1 in ret_endo_list:
                for endo2 in ret_endo_list:
                    if self.Endos[endo1][0] != self.Endos[endo2][0] or self.Endos[endo1][1] != self.Endos[endo2][1] or self.Endos[endo1][2] != self.Endos[endo2][2]:
                        QtWidgets.QMessageBox.critical(self, "Invalid Endonucleases", "The selected endonucleases are not compatible.",QtWidgets.QMessageBox.Ok)
                        return

            GlobalSettings.mainWindow.Results.co_target_endo_list = ret_endo_list
            GlobalSettings.mainWindow.Results.populate_cotarget_table()
            self.cancel_function()
        except Exception as e:
            logger.critical("Error in submission_function() in CoTargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

