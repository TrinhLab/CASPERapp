from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import models.GlobalSettings as GlobalSettings
import traceback
import math
from utils.ui import show_message, show_error, scale_ui, center_ui

logger = GlobalSettings.logger

######################################################
# This class is a window that lets the user select which endonucleases to co-target with
# inputs are from the user and from results
# from results: the organism name and the list of endonucleases for that organism
# from user: which endonucleases to co-target
######################################################
class CoTargeting(QtWidgets.QMainWindow):
    def __init__(self, path):
        try:
            super(CoTargeting, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/cotargeting.ui', self)
            self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Co-targeting")

            self.endo_table.setColumnCount(1)  # hardcoded because there will always be 1 columns
            self.endo_table.setShowGrid(True)
            self.endo_table.setHorizontalHeaderLabels("Endonuclease;".split(";"))
            self.endo_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.endo_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.endo_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.endo_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch) #Ensures last column goes to the edge of table

            self.info_path = path

            self.cancel_button.clicked.connect(self.cancel_function)
            self.submit_button.clicked.connect(self.submission_function)

            scale_ui(self, custom_scale_width=450, custom_scale_height=375)

        except Exception as e:
            show_error("Error initializing CoTargeting class.", e)
            
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

            center_ui(self)

            self.show()
            self.activateWindow()
        except Exception as e:
            show_error("Error in launch() in CoTargeting.", e)
            
    # it clears the table, sets the organism to nothing, and hides the window
    def cancel_function(self):
        try:
            self.endo_table.clearContents()
            self.endo_table.setRowCount(0)
            self.orgName.setText("")
            self.hide()
        except Exception as e:
            show_error("Error in cancel_function() in CoTargeting.", e)
            
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
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical, 
                    title="Nothing Selected", 
                    message="No endonucleases selected. Please select at least 2 endonucleases"
                )
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
                        show_message(
                            fontSize=12,
                            icon=QtWidgets.QMessageBox.Icon.Critical,
                            title="Invalid Endonucleases",
                            message="The selected endonucleases are not compatible."
                        )
                        return

            GlobalSettings.mainWindow.Results.co_target_endo_list = ret_endo_list
            GlobalSettings.mainWindow.Results.populate_cotarget_table()
            self.cancel_function()
        except Exception as e:
            show_error("Error in submission_function() in CoTargeting.", e)