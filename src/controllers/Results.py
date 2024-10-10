from utils.sequence_utils import get_table_headers
from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
from Bio.Seq import Seq
from Bio import SeqIO
from models.CSPRparser import CSPRparser
import models.GlobalSettings as GlobalSettings
import controllers.OffTarget as OffTarget
import platform
import traceback
import math
from controllers.scoring_window import Scoring_Window
from utils.ui import scale_ui, center_ui


#global logger
logger = GlobalSettings.logger

# =========================================================================================
# CLASS NAME: Results
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================
class Results(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        try:
            super(Results, self).__init__(parent)
            uic.loadUi(GlobalSettings.appdir + 'ui/results.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle('Results')
            self.geneViewer.setReadOnly(True)
            self.curgene = ""
            self.dbpath = ""
            # Main Data container
            # Keys: Gene names
            # Values: #
            self.annotation_path = ""
            self.AllData = {}
            self.highlighted = {}
            self.co_target_endo_list = list()
            self.startpos = 0
            self.endpos = 0
            self.directory = ""
            self.inputtype = ""
            self.featureDict = dict() # dictionary passed into transfer_data
            self.featureNTDict = dict() #dictionary passed into transfer_data, same key as featureDict, but hols the NTSEQ
            self.chromDict = dict() # Initialize dictionary for storing chromosome lengths, same key as featureDict
            self.switcher = [1,1,1,1,1,1,1,1]  # for keeping track of where we are in the sorting clicking for each column

            # Initialize Filter Options Object
            self.filter_options = Filter_Options()

            # Initialize Scoring Window Object
            self.scoring_window = Scoring_Window()

            # Target Table settings #
            self.targetTable.setColumnCount(8)  # 
            self.targetTable.setShowGrid(False)
            self.targetTable.setHorizontalHeaderLabels("Location;Endonuclease;Sequence;Strand;PAM;Score;Off-Target;Details".split(";"))
            self.targetTable.horizontalHeader().setSectionsClickable(True)
            self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.targetTable.horizontalHeader().setSectionResizeMode(7, QtWidgets.QHeaderView.Stretch) #Ensures last column goes to the edge of table
            self.targetTable.horizontalHeader().setSectionResizeMode(7, QtWidgets.QHeaderView.Stretch) #Ensures last column goes to the edge of table

            self.back_button.clicked.connect(self.goBack)
            self.targetTable.horizontalHeader().sectionClicked.connect(self.table_sorting)
            self.off_target_button.clicked.connect(self.Off_Target_Analysis)
            self.cotargeting_button.clicked.connect(self.open_cotarget)
            self.displayGeneViewer.stateChanged.connect(self.checkGeneViewer)
            self.filter_options.cotarget_checkbox.stateChanged.connect(self.prep_cotarget_checkbox)
            self.highlight_gene_viewer_button.clicked.connect(self.highlight_gene_viewer)
            self.checkBoxSelectAll.stateChanged.connect(self.selectAll)
            self.filter_options_button.clicked.connect(self.show_filter_options)
            self.scoring_options_button.clicked.connect(self.show_scoring_window)


            self.change_start_end_button.clicked.connect(self.change_indices)
            self.reset_location_button.clicked.connect(self.reset_location)
            self.export_button.clicked.connect(self.open_export_tool)

            #self.targetTable.itemSelectionChanged.connect(self.item_select)
            self.filter_options.minScoreLine.setText("0")

            # Connecting the filters to the displayGeneData function
            self.filter_options.fivegseqCheckBox.stateChanged.connect(self.displayGeneData)
            self.filter_options.minScoreLine.textChanged.connect(self.displayGeneData)

            # Setting up the score filter:
            self.filter_options.scoreSlider.setMinimum(0)
            self.filter_options.scoreSlider.setMaximum(100)
            self.filter_options.scoreSlider.setTracking(False)
            self.filter_options.scoreSlider.valueChanged.connect(self.update_score_filter)

            #bool used to make sure only 1 instance of the OffTarget window is created
            self.first_boot = True
            #OTA is used to hold the row numbers of the items selected by user for OffTargetAnalysis
            #using this helps speed up updating the chart
            self.OTA = []

            self.clear_highlighted_guides_button.clicked.connect(self.clear_highlighted_guides)

            self.detail_output_list = []
            self.rows_and_seq_list = []
            self.seq_and_avg_list = []
            self.files_list = []
            self.mwfg = self.frameGeometry()  ##Center window
            self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window

            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#guide_viewer{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            margin-top: 10px;
                            font: bold 14pt 'Arial';}"""

            self.guide_viewer.setStyleSheet(groupbox_style)
            self.guide_analysis.setStyleSheet(groupbox_style.replace("guide_viewer", "guide_analysis"))
            self.gene_viewer.setStyleSheet(groupbox_style.replace("guide_viewer", "gene_viewer"))

            self.get_endo_data()

            ### Make line edits only accept integers
            self.lineEditStart.setValidator(QtGui.QIntValidator())
            self.lineEditEnd.setValidator(QtGui.QIntValidator())


            #scale UI
            self.first_show = True
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing results class.")
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
            fontSize = 30
            self.title.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 850x750
            scaledWidth = int((width * 1250) / 1920)
            scaledHeight = int((height * 750) / 1080)

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
            logger.critical("Error in scaleUI() in results.")
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
            logger.critical("Error in centerUI() in results.")
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

    def get_endo_data(self):
        try:
            f = open(GlobalSettings.appdir + "CASPERinfo")
            self.endo_data = {}
            while True:
                line = f.readline()
                if line.startswith('ENDONUCLEASES'):
                    while True:
                        line = f.readline()
                        line = line.replace("\n","")
                        if (line[0] == "-"):
                            break
                        line_tokened = line.split(";")
                        if len(line_tokened) == 10:
                            endo = line_tokened[0]
                            five_length = line_tokened[2]
                            seed_length = line_tokened[3]
                            three_length = line_tokened[4]
                            prime = line_tokened[5]
                            hsu = line_tokened[9]
                            self.endo_data[endo] = [int(five_length) + int(three_length) + int(seed_length), prime, "MATRIX:" + hsu]

                    break
            f.close()
        except Exception as e:
            logger.critical("Error in get_endo_data() in results.")
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

    # This function resets (turns off, then on) the selection of all the rows in the ViewTargets table.
    # For some reason this is necessary after adding values to the table through alternative scoring
    # (different on-target scoring or off-target), otherwise the export function will bug due to newly
    # added items being queued to the back of the list of self.targetTable.selectedItems(), instead of
    # where they belong at the end of each row.
    def reset_selection(self):
        rows = sorted(set(index.row() for index in self.targetTable.selectedIndexes())) # Find selected rows
        self.targetTable.clearSelection() # Clear the selection
        for row in rows: # For each selected row
            self.targetTable.selectRow(row) # Reselect each row

    # this function opens the export_tool window
    # first it makes sure that the user actually has some highlighted targets that they want exported
    def open_export_tool(self):
        try:
            self.reset_selection()
            select_items = self.targetTable.selectedItems()
            if len(select_items) <= 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Nothing Selected")
                msgBox.setText("No targets were highlighted. Please highlight the targets you want to be exported to a CSV File!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return
            # now launch the window
            GlobalSettings.mainWindow.export_tool_window.launch(select_items,"vt")
        except Exception as e:
            logger.critical("Error in open_export_tool() in results.")
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

    def change_indices(self):
        try:
            ### Make sure the gene viewer is on
            if not self.displayGeneViewer.isChecked():
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Gene Viewer Error")
                msgBox.setText("Gene Viewer display is off! Please turn the Gene Viewer on in order to highlight the sequences selected")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            ### Change the start and end values
            prevTuple = self.featureDict[self.curgene]
            tempTuple = (self.featureDict[self.curgene][0], int(self.lineEditStart.displayText())-1, int(self.lineEditEnd.displayText()))

            ### Make sure both indices are greater than 0
            if tempTuple[1]+1 <= 0 or tempTuple[2] <= 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Invalid location indices.")
                msgBox.setText("Location indices cannot be negative or zero! Please set values larger than 0.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1))
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2]))
                self.reset_location()
                return

            ### Make sure start is less than stop
            if tempTuple[1] >= tempTuple[2]:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Invalid location indices.")
                msgBox.setText("Start location must be less than stop location.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1))
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2]))
                self.reset_location()
                return

            ### Make sure that the difference between indicies is not too large
            if abs(tempTuple[1] - tempTuple[2]) > 50000:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Sequence Too Long")
                msgBox.setText("The sequence is too long! Please choose indicies that will make the sequence less than 50,000!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1))
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2]))
                self.reset_location()
                return

            ### Make sure search is within chromosome range 
            if int(tempTuple[2]) > self.chromDict[self.curgene]:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Position Error: Region not in Chromosome")
                msgBox.setText(
                    "The stop location is greater than the chromosome's length.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1))
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2]))
                self.reset_location()
                return

            sequence, chrom_len = self.sequence_finder(tempTuple) # Get the appropriate NT sequence
            self.geneViewer.setText(sequence) # Set the gene viewer to display the sequence
        except Exception as e:
            logger.critical("Error in change_indices() in results.")
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

    # this function listens for a stateChange in selectAllShown
    # if it is checked, it selects all shown
    # if it is unchecked, it deselects all shown
    # Note: it is a little buggy, possibly because when you change the minimum score it resets it all
    def selectAll(self):
        try:
            if self.checkBoxSelectAll.isChecked():
                self.targetTable.selectAll()
            else:
                self.targetTable.clearSelection()
        except Exception as e:
            logger.critical("Error in selectAll() in results.")
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
    ### This function resets gene viewer to the appropriate sequence
    def reset_location(self):
        try:
            if self.displayGeneViewer.isChecked():
                self.geneViewer.setText(self.featureNTDict[self.curgene])
                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1))
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2]))
            else:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Gene Viewer Error")
                msgBox.setText("Gene Viewer display is off! Please turn the Gene Viewer on in order to reset the locations")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()
        except Exception as e:
            logger.critical("Error in reset_location() in results.")
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




    # hightlights the sequences found in the gene viewer
    # highlighting should stay the exact same with fasta and genbank files, as this function only edits what
    #   is currently in the gene viewer text table anyways
    def highlight_gene_viewer(self):
        try:
            # make sure gene viewer is enabled
            if not self.displayGeneViewer.isChecked():
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Gene Viewer Error")
                msgBox.setText("Gene Viewer display is off! Please turn the Gene Viewer on in order to highlight the sequences selected")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            # variables needed
            cursor = self.geneViewer.textCursor()
            format = QtGui.QTextCharFormat()
            failed_guides = []

            # check and make sure still is actually highlighted!
            selectedList = self.targetTable.selectedItems()
            if len(selectedList) <= 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Nothing Selected")
                msgBox.setText("No targets were highlighted. Please highlight the targets you want to be highlighted in the gene viewer!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return
            # this is the loop that actually goes in and highlights all the things
            for i in range(self.targetTable.rowCount()):
                if self.targetTable.item(i, 0).isSelected():
                    # get the strand and sequence strings
                    locationString = self.targetTable.item(i,0).text()
                    strandString = self.targetTable.item(i, 3).text()
                    sequenceString = self.targetTable.item(i, 2).text()
                    printSequence = ""
                    movementIndex = len(sequenceString)
                    left_right = ""
                    #print("Length of geneViewer: ", len(self.geneViewer.toPlainText()))

                    if strandString == "+":
                        format.setBackground(QtGui.QBrush(QtGui.QColor("green")))
                        index = self.geneViewer.toPlainText().upper().find(str(sequenceString))
                        if index != -1: # If gRNA is found in GeneViewer
                            cursor.setPosition(index)
                            for i in range(movementIndex): # Actually highlight the gRNA now
                                cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                            cursor.mergeCharFormat(format)
                        else:
                            failed_guides.append(sequenceString)
                            continue
                    else: # gRNA is on negative strand
                        format.setBackground(QtGui.QBrush(QtGui.QColor("red")))
                        index = self.geneViewer.toPlainText().upper().find(str(Seq(sequenceString).reverse_complement()))
                        if index != -1: # If gRNA is found in GeneViewer
                            cursor.setPosition(index)
                            for i in range(movementIndex): # Actually highlight the gRNA now
                                cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                            cursor.mergeCharFormat(format)
                        else:
                            failed_guides.append(sequenceString)
                            continue

            # if any of the sequences return 0 matches, show the user which ones were not found
            if len(failed_guides) > 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Warning")
                msgBox.setText(
                    "The following sequence(s) were not found in the Gene Viewer text:\n\t" + "\n".join(failed_guides))
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

        except Exception as e:
            logger.critical("Error in highlight_gene_viewer() in results.")
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

    # this function updates the gene viewer based on the user clicking 'display on'
    # if it is check marked, it displays the correct data
    # if it is un-marked, it hides the data
    def checkGeneViewer(self):
        try:
            if self.displayGeneViewer.isChecked():
                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1))
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2]))
                self.geneViewer.setText(self.featureNTDict[self.curgene])
            elif not self.displayGeneViewer.isChecked():
                self.lineEditStart.clear()
                self.lineEditEnd.clear()
                self.geneViewer.clear()
        except Exception as e:
            logger.critical("Error in checkGeneViewer() in results.")
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

    # this function opens when the user clicks the CoTargeting button
    def open_cotarget(self):
        try:
            endo_list = list()
            if self.endonucleaseBox.count() <= 1:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Not Enough Endonucleases")
                msgBox.setText(
                    "There are not enough endonucleases with this organism. At least 2 endonucleases are required for this function. Use Analyze New Genome to create CSPR files with other endonucleases.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            for i in range(self.endonucleaseBox.count()):
                endo_list.append(self.endonucleaseBox.itemText(i))

            GlobalSettings.mainWindow.CoTargeting.launch(endo_list, GlobalSettings.mainWindow.orgChoice.currentText())
        except Exception as e:
            logger.critical("Error in open_cotarget() in results.")
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

    # this function goes through and calls transfer_data again.
    # Uses data from the mainWindow in Globalsettings, but that's because that info should not change
    # unless the user closes out of the Results window
    def changeEndonuclease(self):
        try:
            full_org = str(GlobalSettings.mainWindow.orgChoice.currentText())
            # organism = GlobalSettings.mainWindow.shortHand[full_org]

            endoChoice = self.endonucleaseBox.currentText().split("|")

            # make sure the user actually selects a new endonuclease
            if self.endo == endoChoice:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Select a different Endonuclease")
                msgBox.setText(
                    "Please be sure to select a different endonuclease!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            # enable the cotarget checkbox if needed
            if len(endoChoice) > 1:
                self.filter_options.cotarget_checkbox.setEnabled(True)
                self.filter_options.cotarget_checkbox.setChecked(0)
            else:
                self.filter_options.cotarget_checkbox.setEnabled(False)
                self.filter_options.cotarget_checkbox.setChecked(0)
            self.transfer_data(full_org, GlobalSettings.mainWindow.organisms_to_files[full_org], endoChoice, GlobalSettings.CSPR_DB, self.featureDict,
                               self.featureNTDict, self.inputtype)
        except Exception as e:
            logger.critical("Error in changeEndonuclease() in results.")
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

    # Function that is used to set up the results page.
    # it calls get_targets, which in turn calls display data
    def transfer_data(self, org, org_files, endo, path, feature_dict, featureNTSeqDict, inputtype):
        try:
            # set all of the classes variables
            self.org = org
            self.org_files = org_files
            self.endo = endo
            self.directory = path
            self.comboBoxGene.clear()
            self.AllData.clear()
            self.featureDict =feature_dict
            self.featureNTDict = featureNTSeqDict
            self.inputtype = inputtype

            self.highlighted.clear()
            self.detail_output_list.clear()
            self.seq_and_avg_list.clear()
            self.rows_and_seq_list.clear()
            self.OTA.clear()

            for feature in feature_dict:
                if self.inputtype == "feature":
                    detail_output1 = {}
                    rows_and_seq2 = {}
                    seq_and_avg3 = {}
                    self.detail_output_list.append(detail_output1)
                    self.seq_and_avg_list.append(seq_and_avg3)
                    self.rows_and_seq_list.append(rows_and_seq2)
                    self.comboBoxGene.addItem(feature)
                    self.get_targets(feature, feature_dict[feature])
                if self.inputtype == "position":
                    detail_output1 = {}
                    rows_and_seq2 = {}
                    seq_and_avg3 = {}
                    self.detail_output_list.append(detail_output1)
                    self.seq_and_avg_list.append(seq_and_avg3)
                    self.rows_and_seq_list.append(rows_and_seq2)
                    self.comboBoxGene.addItem(feature)
                    self.get_targets(feature, feature_dict[feature])
                if self.inputtype == "sequence":
                    detail_output1 = {}
                    rows_and_seq2 = {}
                    seq_and_avg3 = {}
                    self.detail_output_list.append(detail_output1)
                    self.seq_and_avg_list.append(seq_and_avg3)
                    self.rows_and_seq_list.append(rows_and_seq2)
                    self.comboBoxGene.addItem(feature)
                    self.get_targets(feature, feature_dict[feature])


            # Enable the combobox to be toggled now that the data is in AllData
            self.comboBoxGene.currentTextChanged.connect(self.displayGeneData)
            self.first_boot = True
        except Exception as e:
            logger.critical("Error in transfer_data() in results.")
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

    def goBack(self):
        try:
            GlobalSettings.mainWindow.show()
            self.filter_options.cotarget_checkbox.setChecked(0)
            self.filter_options.hide()
            try:
                self.off_tar_win.hide()
            except:
                pass
            GlobalSettings.mainWindow.CoTargeting.hide()
            self.scoring_window.hide()
            self.hide()

        except Exception as e:
            logger.critical("Error in goBack() in results.")
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

    # called when the user hits 'gene viewer settings'
    def changeGeneViewerSettings(self):
        try:
            GlobalSettings.mainWindow.gene_viewer_settings.show()
            GlobalSettings.mainWindow.gene_viewer_settings.activateWindow()
        except Exception as e:
            logger.critical("Error in changeGeneViewerSettings() in results.")
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

    # this is the function that sets up the cotargeting.
    # it is called from the Cotargeting class, when the user hits submit
    # myBool is whether or not to change the endoChoice comboBox
    def populate_cotarget_table(self, myBool = True):
        try:
            try:
                self.endonucleaseBox.currentIndexChanged.disconnect()
            except:
                pass

            # make a string of the combination, separated by commas's
            endoBoxString = ""
            for i in range(len(self.co_target_endo_list)):
                if endoBoxString == "":
                    endoBoxString = self.co_target_endo_list[i]
                else:
                    endoBoxString = endoBoxString + '|' + self.co_target_endo_list[i]

            # put the new endoChoice at the beginning. This is the only way i could find to do it
            # get a list of all endo choices, and put the newest at the front
            endoBoxList = list()
            endoBoxList.append(endoBoxString)
            for i in range(self.endonucleaseBox.count()):
                if self.endonucleaseBox.itemText(i) not in endoBoxList: # Prevent duplicate entries
                    endoBoxList.append(self.endonucleaseBox.itemText(i))

            # clear the current endo choices, and append the new order
            if myBool:
                self.endonucleaseBox.clear()
                for i in range(len(endoBoxList)):
                    self.endonucleaseBox.addItem(endoBoxList[i])

            # enable the cotarget checkbox
            self.filter_options.cotarget_checkbox.setEnabled(True)
            self.filter_options.cotarget_checkbox.setChecked(0)

            self.endonucleaseBox.currentIndexChanged.connect(self.changeEndonuclease)
            # add it to the endoBox choices, and then call transfer_data
            self.transfer_data(self.org, self.org_files, self.co_target_endo_list, GlobalSettings.CSPR_DB, self.featureDict, self.featureNTDict,self.inputtype)
        except Exception as e:
            logger.critical("Error in populate_cotarget_table() in results.")
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

    # prep function for the checkbox for cotargeting
    # if the checkbox is checked, just go ahead and displayGeneData
    # if not, call populate_cotarget_table, as a reset to get all of the data there
    def prep_cotarget_checkbox(self):
        try:
            if self.filter_options.cotarget_checkbox.isChecked():
                self.displayGeneData()
            elif not self.filter_options.cotarget_checkbox.isChecked():
                self.populate_cotarget_table(myBool=False)
        except Exception as e:
            logger.critical("Error in prep_cotarget_checkbox() in results.")
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

    # Function grabs the information from the .cspr file and adds them to the AllData dictionary
    #changed to now call CSPRparser's function. Same function essentially, just cleaned up here
    def get_targets(self, genename, pos_tuple):
        try:
            #get the right files
            for endo in self.endo:
                if platform.system() == "Windows":
                    file = self.directory + "\\" + self.org_files[endo][0]
                else:
                    file = self.directory + "/" + self.org_files[endo][0]

                #create the parser, read the targets store it. then display the GeneData screen
                parser = CSPRparser(file)

                # if genename is not in the dict, make that spot into a list
                if genename not in self.AllData:
                    self.AllData[genename] = list()
                # now append parser's data to it
                self.AllData[genename].append(parser.read_targets(genename, pos_tuple, endo))

                # for each list item
                for item in self.AllData[genename]:
                    # for each tuple item
                    for i in range(len(item)):
                        self.highlighted[item[i][1]] = False

            # if the endo choice is greater than 1, call the combine
            if len(self.endo) > 1:
                self.combine_coTargets(genename)
            self.displayGeneData()
        except Exception as e:
            logger.critical("Error in get_targets() in results.")
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

    ###############################################################################################################
    # Main Function for updating the Table.  Connected to all filter buttons and the Gene toggling of the combobox.
    ###############################################################################################################
    def displayGeneData(self):
        try:
            self.curgene = str(self.comboBoxGene.currentText())  # Gets the current gene
            # Creates the set object from the list of the current gene:
            if self.curgene=='' or len(self.AllData)<1:
                return

            subset_display = []
            # set the start and end numbers, as well as set the geneViewer text, if the displayGeneViewer is checked
            if self.displayGeneViewer.isChecked():
                self.lineEditStart.setText(str(self.featureDict[self.curgene][1]+1)) # Set start index for gene (Add 1 to account for Python indexing convention)
                self.lineEditEnd.setText(str(self.featureDict[self.curgene][2])) # Set end index for gene
                self.geneViewer.setText(self.featureNTDict[self.curgene]) # Get gene sequence

            # if this checkBox is checked, remove the single endo
            if self.filter_options.cotarget_checkbox.isChecked():
                gene = self.curgene
                self.remove_single_endo(gene)

            # Removing all sequences below minimum score and creating the set:
            # for each list item
            for item in self.AllData[self.curgene]:
                # for each tuple item
                for i in range(len(item)):
                    if int(item[i][3]) > int(self.filter_options.minScoreLine.text()):
                        # Removing all non 5' G sequences:
                        if self.filter_options.fivegseqCheckBox.isChecked():
                            if item[i][1].startswith("G"):
                                subset_display.append(item[i])
                        else:
                            subset_display.append(item[i])

            self.targetTable.setRowCount(len(subset_display))
            index = 0
            #changed the number items to use setData so that sorting will work correctly
            #because before the numbers were interpretted as strings and not numbers


            ### Remove alternative scoring columns (Azimuth, etc.) when switching between genes or loading new data...this prevents carry-over of the wrong scores from previously scored genes
            ### One possible solution would be to add alternate scores to self.AllData, but I won't do that for now.
            header = get_table_headers(self.targetTable) # Returns headers of the target table
            num_cols = len(header) # Get number of columns
            col_indices = [header.index(x) for x in GlobalSettings.algorithms if x in header] # Returns the index(es) of the alternative scoring column(s) in the target table of View Targets window
            if len(col_indices) > 0: # If alternative scoring has been done
                for i in col_indices:
                    self.targetTable.removeColumn(i)

            for item in subset_display:
                num = int(item[0])
                loc = QtWidgets.QTableWidgetItem()
                loc.setData(QtCore.Qt.EditRole, abs(num))
                seq = QtWidgets.QTableWidgetItem(item[1])
                strand = QtWidgets.QTableWidgetItem(str(item[4]))
                PAM = QtWidgets.QTableWidgetItem(item[2])
                num1 = int(item[3])
                endonuclease = QtWidgets.QTableWidgetItem(item[5])
                score = QtWidgets.QTableWidgetItem()
                score.setData(QtCore.Qt.EditRole, num1)
                self.targetTable.setItem(index, 0, loc)
                self.targetTable.setItem(index, 1, endonuclease)
                self.targetTable.setItem(index, 2, seq)
                self.targetTable.setItem(index, 3, strand)
                self.targetTable.setItem(index, 4, PAM)
                self.targetTable.setItem(index, 5, score)
                self.targetTable.setItem(index, 6, QtWidgets.QTableWidgetItem("--.--")) # Give "blank" value for Off-Target
                self.targetTable.removeCellWidget(index, num_cols-1) # Leave the "Details" column empty
                if (item[1] in self.seq_and_avg_list[self.comboBoxGene.currentIndex()].keys()):
                    OT = QtWidgets.QTableWidgetItem()
                    OT.setData(QtCore.Qt.EditRole, self.seq_and_avg_list[self.comboBoxGene.currentIndex()][item[1]])
                    self.targetTable.setItem(index, 6, OT)
                if (item[1] in self.detail_output_list[self.comboBoxGene.currentIndex()].keys()):
                    details = QtWidgets.QPushButton()
                    details.setText("Details")
                    details.clicked.connect(self.show_details)
                    self.targetTable.setCellWidget(index, 7, details)
                index += 1
            self.targetTable.resizeColumnsToContents()
        except Exception as e:
            logger.critical("Error in displayGeneData() in results.")
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

    # this function is only entered if the user checks the show only cotargeted sequence checkbox
    def remove_single_endo(self, gene):
        try:
            removalDict = dict()
            # go through and figure out which ones need to be shown
            for i in range(len(self.AllData[gene])):
                for j in range(len(self.AllData[gene][i])):
                    endoData = self.AllData[gene][i][j][5].split("|")
                    if len(endoData) == 1:
                        if i not in removalDict:
                            removalDict[i] = list()

                        removalDict[i].append(j)

            # now go through and delete them. But you have to go in reverse
            for item in removalDict:
                # for the reverse of that list. This is to keep the program from crashing
                # easier than building a new list honestly
                for index in reversed(removalDict[item]):
                    self.AllData[gene][item].pop(index)
        except Exception as e:
            logger.critical("Error in remove_single_endo() in results.")
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

    # this function goes through and combines table rows that have the same location and PAM dir
    # it edits the dictionary data itself
    # currently it does not take the PAM direction into account
    # parameter genename:  the key to which part in the dictionary to look at
    def combine_coTargets(self, genename):
        try:
            cotarget_data = list()
            delete_list = list()
            endoList = list()
            endo_lengths = list()

            #get endo data
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
                        self.Endos[endo] = line_tokened[5]
                    break
            f.close()

            ### Get the endo list data
            for endo_data in self.AllData[genename]: # Loop through each endonuclease
                ### If there is no data for an endo, just return because co-targeting is useless for that one
                if len(endo_data) == 0:
                    return
                if len(endoList) == 0: # If nothing has been added to the endo list yet
                    endoList.append(endo_data[0][5]) # Add endo name to list
                    endo_lengths.append(len(endo_data[0][2])) # Add endo length to list
                elif len(endo_data[0][2]) > endo_lengths[-1]: # Put the endo with the longest PAM first in the list
                    endoList.insert(0,endo_data[0][5])
                    endo_lengths.insert(0,endo_data[0][5])
                else:
                    endoList.append(endo_data[0][5]) # Add endo name to end of list
                    endo_lengths.append(len(endo_data[0][2])) # Add endo length to end of list
            cotarget_endo = "|".join(endoList)

            ### Get Co-Targets!
            for i, endo_data in enumerate(self.AllData[genename]): # For each endo in the genename block
                if i == 0: # If first endonuclease checked, just append data to lists
                    for target_tuple in endo_data: # For each target 
                        cotarget_data.append(target_tuple) # Store tuple information in list
                else:
                    continue
            for i, endo_data in enumerate(self.AllData[genename]): # For each endo in the genename block
                if i == 0: # If first endonuclease checked, just continue
                   continue
                else: # If not first endo, check and see if any sequences are cotargets
                    for j,cotarget in enumerate(cotarget_data): # For each potential cotarget
                        tmp_list = [x for x in endo_data if cotarget[1] == x[1]] # Check if cotarget is also in this endo's targets
                        if len(tmp_list) > 0: # If cotarget was found
                            new_tuple = tmp_list[0] # cotarget tuple
                            if len(new_tuple[2]) > len(cotarget_data[j][2]): # If PAM is longer for this endo...
                                cotarget_data[j] = (new_tuple[0],new_tuple[1],new_tuple[2],new_tuple[3],new_tuple[4],cotarget_endo) # Overwrite existing entry with the one containing longest PAM
                            else:
                                cotarget_data[j] = (cotarget[0],cotarget[1],cotarget[2],cotarget[3],cotarget[4],cotarget_endo) # Overwrite existing entry to have the right endonuclease
                        else: # If not cotarget
                            delete_list.append(j)

            ### Delete ones from cotarget_data that failed (can't do this in the above loop for some reason)
            final_targets = [x for i,x in enumerate(cotarget_data) if i not in delete_list]
            self.AllData[genename].append(final_targets)

        except Exception as e:
            logger.critical("Error in combine_coTargets() in results.")
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

    def search_gene(self):
        try:
            search_trms = []
            checkBox = self.sender()
            index = self.targetTable.indexAt(checkBox.pos())
            seq = self.targetTable.item(index.row(),2).text()
            self.highlighted[str(seq)] = checkBox.isChecked()

            x=1
        except Exception as e:
            logger.critical("Error in search_gene() in results.")
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

    def item_select(self):
        try:
            print(self.targetTable.selectedItems())
        except Exception as e:
            logger.critical("Error in item_select() in results.")
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

    def getColumnIndexByHeaderName(self, table, header_name):
        # This functin takes a QTableWidget object and header (column) name and returns the index in the table.
        header = table.horizontalHeader()
        for column in range(header.count()):
           logical_index = header.logicalIndex(column)
           if header.model().headerData(logical_index, header.orientation()) == header_name:
               return logical_index
        return None
               
    def table_sorting(self, logicalIndex):
        try:
            if logicalIndex == self.getColumnIndexByHeaderName(self.targetTable, "Details"): # Prevent sorting by Details column to prevent crashing
                return
            self.switcher[logicalIndex] *= -1
            if self.switcher[logicalIndex] == -1:
                self.targetTable.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.targetTable.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)
        except Exception as e:
            logger.critical("Error in table_sorting() in results.")
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

    def update_score_filter(self):
        try:
            self.filter_options.minScoreLine.setText(str(self.filter_options.scoreSlider.value()))
        except Exception as e:
            logger.critical("Error in update_score_filter() in results.")
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

    #linked to when the user pushes tools->off target analysis
    def Off_Target_Analysis(self):
        try:
            #build temp file for offtarget to read from
            if platform.system() == 'Windows':
                f = open(GlobalSettings.appdir + 'OffTargetFolder' + '\\temp.txt','w+')
            else:
                f = open(GlobalSettings.appdir + 'OffTargetFolder' + '/temp.txt', 'w+')
            self.OTA.clear()
            #get selected rows
            selected_rows = []
            indexes = self.targetTable.selectionModel().selectedRows()

            if len(indexes) == 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No Rows Selected")
                msgBox.setText(
                    "Please select rows from the table to pass into the off-target analysis!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            for index in indexes:
                selected_rows.append(index.row())

            for row in sorted(selected_rows):
                self.OTA.append(row)
                loc = self.targetTable.item(row, 0).text()
                seq = self.targetTable.item(row,2).text()
                strand = self.targetTable.item(row,3).text()
                pam = self.targetTable.item(row,4).text()
                score = self.targetTable.item(row,5).text()
                self.rows_and_seq_list[self.comboBoxGene.currentIndex()][seq] = row
                output = str(loc) + ';' + str(seq) + ";" + str(pam) + ";" + score + ";" + str(strand)
                f.write(output + '\n')
            f.close()
            #only make off target object if first time, otherwise just
            #reshow the object
            if self.first_boot:
                self.first_boot = False
                self.off_tar_win = OffTarget.OffTarget()
                self.off_tar_win.submitButton.clicked.connect(self.refresh_data)
            center_ui(self.off_tar_win)
            ref_org = str(GlobalSettings.mainWindow.orgChoice.currentText()) ### Set default reference organism to the organism that is being targeted
            index = self.off_tar_win.OrgcomboBox.findText(ref_org) ### Find organism in combo box list
            self.off_tar_win.OrgcomboBox.setCurrentIndex(index) ### Set combo box to appropriate index
            self.off_tar_win.show()
            self.off_tar_win.activateWindow()
            f.close()
        except Exception as e:
            logger.critical("Error in Off_Target_Analysis() in results.")
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

    #refresh data is linked to the submit button on the off target analysis UI
    def refresh_data(self):
        try:
            #setup filename based on output name given in OffTarget
            filename = self.off_tar_win.output_path

            # if the user hits submit without running the program, do nothing
            if not self.off_tar_win.run_clicked:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("File Not Found")
                msgBox.setText(
                    "There was an error with the Off Target execution. No results file was found.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()
                return

            self.off_tar_win.hide()
            filename = filename[:len(filename)-1]
            filename = filename[1:]
            filename = filename.replace(r'\\', '\\')
            filename = filename.replace('"', '')
            self.files_list.append(filename)
            try:
                out_file = open(filename, "r")
            except:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Unable to Open File")
                msgBox.setText(
                    "There was an error with the Off Target execution. The results file was either not found or not able to be opened.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return

            #read the first line : either AVG or DETAILED OUTPUT
            output_type = out_file.readline()
            output_type = output_type.strip('\r\n')

            #parse based on whether avg or detailed output
            line_cnt = 0
            headers = get_table_headers(self.targetTable) # Get table headers
            num_cols = len(headers) # Get number of columns
            if(output_type == "AVG OUTPUT"):
                for line in out_file:
                    line = line.strip('\n')
                    if (line != ''):
                        values = line.split(":")
                        row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                        OT = QtWidgets.QTableWidgetItem()
                        OT.setData(QtCore.Qt.EditRole, values[1])
                        self.targetTable.setItem(row, num_cols-2, OT) # Set the OT score to the second to last column
                        self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
                        line_cnt += 1
            else:
                details_bool = False
                temp_list = []
                values = []
                for line in out_file:
                    line = line.strip('\n')
                    if(line.find(':') != -1):
                        if(details_bool == True):
                            self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list
                            details_bool = False
                            temp_list = []
                        values = line.split(":")
                        self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
                        row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                        OT = QtWidgets.QTableWidgetItem()
                        OT.setData(QtCore.Qt.EditRole, values[1])
                        self.targetTable.setItem(row, num_cols - 2, OT)
                        line_cnt += 1
                    elif line != "":
                        details_bool = True
                        temp_list.append(line)
                        details = QtWidgets.QPushButton()
                        details.setText("Details")
                        details.clicked.connect(self.show_details)
                        self.targetTable.setCellWidget(row, num_cols - 1, details)
                        line_cnt += 1
                if(details_bool == True):
                    self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list


                #make sure OT output file had lines
                if line_cnt < 1:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("File Empty")
                    msgBox.setText(
                        "There was an error with the Off Target execution. No results were found in the results file.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()

                    return

                #print(self.detail_output_list)
                self.targetTable.resizeColumnsToContents()
                out_file.close()
        except Exception as e:
            logger.critical("Error in refresh_data() in results.")
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

    def show_details(self):
        try:
            #create msg box popup for the details | used html to make it easier to style with bold
            button = self.sender()
            index = self.targetTable.indexAt(button.pos())
            msg = QtWidgets.QMessageBox()
            msg.setWindowTitle("Details")
            msg.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            key = str(self.targetTable.item(index.row(),2).text())
            temp_str = ''
            for items in self.detail_output_list[self.comboBoxGene.currentIndex()][key]:
                temp_str += items + "<br>"

            chromo_str = "<html><b>Reference gRNA:</b><br>Location, Sequence, Strand, PAM, On Score<br></html>"
            input_str = self.targetTable.item(index.row(),0).text() + ', ' + key + ', ' + \
                        self.targetTable.item(index.row(),3).text() + ', ' + self.targetTable.item(index.row(),4).text() + \
                        ', ' + self.targetTable.item(index.row(),5).text() + "<br><br>"
            detail_str = "<html><b>Off-Target Hits:</b><br>Off Score, Chromosome, Location, Sequence<br></html>"
            msg.setText(chromo_str + input_str + detail_str + temp_str)
            msg.exec()

        except Exception as e:
            logger.critical("Error in show_details() in results.")
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

    #function to clear highlights in gene viewer and un-select any rows highlighted in main table
    def clear_highlighted_guides(self):
        try:
            # Clear guides from gene viewer
            self.change_indices()

            # Clear selected rows in table
            self.checkBoxSelectAll.setChecked(False)
            self.targetTable.selectionModel().clearSelection()
        except Exception as e:
            logger.critical("Error in clear_highlighted_guides() in results.")
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

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        try:
            GlobalSettings.mainWindow.closeFunction()
            event.accept()
        except Exception as e:
            logger.critical("Error in closeEvent() in results.")
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

    ############################################
    # All Gene Viewer functions are below!
    ############################################
    # It will go based on the lengths stored in the comboGeneBox dictionary
    def load_gene_viewer(self):
        try:
            if GlobalSettings.mainWindow.annotation_files.currentText() != "None":
                sequence = ""
                # for each gene selected from the results window
                for item in self.featureDict:
                    sequence, chrom_len = self.sequence_finder(self.featureDict[item])
                    self.featureNTDict[item] = sequence
                    self.chromDict[item] = chrom_len
                self.lineEditStart.setEnabled(True)
                self.lineEditEnd.setEnabled(True)
                self.change_start_end_button.setEnabled(True)
                self.displayGeneViewer.setEnabled(True)
                self.displayGeneViewer.setChecked(0)
                self.checkGeneViewer()
            else:
                self.displayGeneViewer.setEnabled(False) # Disable Gene Viewer if no annotation file has been selected
                return
        except Exception as e:
            logger.critical("Error in load_gene_viewer() in results.")
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

    # This function takes the location data for a feature- or position-based search
    # and returns the appropriate NT sequence from a GenBank formatted annotation file using BioPython's SeqIO.parse() method
    def sequence_finder(self, location_data):
        try:
            ### Start up the function
            chrom_index = location_data[0]-1 # This is the chromosome we need to pull sequence data from. (Python indexing, so chromosome 1 is index 0)
            start = location_data[1]
            end = location_data[2]
            ### Pull the sequence information from the GenBank file if present
            if self.annotation_path != "":
                parser = SeqIO.parse(self.annotation_path,'genbank') # Initialize parser object for GenBank file
                for i,record in enumerate(parser): # Loop through chromosomes 
                    if chrom_index == i: # If this is the correct chromosome
                        chrom_seq = str(record.seq).strip()
                        ### Get appropriate sequence and padding (for visualizing gRNAs that appear at extreme ends of region)
                        if (start - 30) >= 0: # Check to make sure there is enough 5' end of gene to pull the padding from, so indexing error isn't raised
                            five_prime_tail = chrom_seq[(start-30):start]
                        else:
                            five_prime_tail = ""
                        if len(chrom_seq) >= (end + 30): # Check to make sure there is enough 3' end of gene to pull the padding from, so indexing error isn't raised
                            three_prime_tail = chrom_seq[end:end+30]
                        else:
                            three_prime_tail = ""
                        my_seq = chrom_seq[start:end] # Get the sequence from the specified location
                        ret_sequence = five_prime_tail.lower() + my_seq.upper() + three_prime_tail.lower() # Add padding to the sequence
                        return ret_sequence, len(chrom_seq)
                    else: # If this is not the right chromosome, go to the next one
                        continue
            else:
                return

        except Exception as e:
            logger.critical("Error in sequence_finder() in results.")
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

    #show filters UI
    def show_filter_options(self):
        try:
            self.filter_options.centerUI()
            self.filter_options.show()
            self.filter_options.activateWindow()
        except Exception as e:
            logger.critical("Error in show_filter_options() in results.")
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

    #show scoring window UI 
    def show_scoring_window(self):
        try:
            # Check to make sure gRNAs were highlighted
            selectedList = self.targetTable.selectedItems()
            if len(selectedList) <= 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Nothing Selected")
                msgBox.setText("No guides were highlighted. Please highlight the guides you want to score!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()
                return
            else:
                self.scoring_window.scaleUI()
                self.scoring_window.centerUI()
                self.scoring_window.show()
                self.scoring_window.activateWindow()

        except Exception as e:
            logger.critical("Error in show_scoring_window() in results.")
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

        except Exception as e:
            logger.critical("Error in show_scoring_window() in results.")
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

class Filter_Options(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        try:
            super(Filter_Options, self).__init__(parent)
            uic.loadUi(GlobalSettings.appdir + 'ui/filter_options.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Filter Options")
            self.minScoreLine.setText("0")
            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#filterBox{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            margin-top: 10px;}"""
            self.filterBox.setStyleSheet(groupbox_style)

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing Filter_Options class in results.")
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

    # scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            # 16px is used for 92 dpi / 1920x1080
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
            scaledWidth = int((width * 350) / 1920)
            scaledHeight = int((height * 250) / 1080)

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
            logger.critical("Error in scaleUI() in filter options in results.")
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

    # center UI on current screen
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
            logger.critical("Error in centerUI() in filter options in results.")
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

