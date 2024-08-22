from utils.Algorithms import get_table_headers
import models.GlobalSettings as GlobalSettings
import azimuth.model_comparison as az
from PyQt5 import QtWidgets, uic, QtCore
import platform
import traceback
import math
from Bio import SeqIO
import numpy as np
import warnings
import contextlib
import sys

""" Code for suppressing stdout of azimuth algorithm """
class DummyFile(object):
    def write(self, x): pass

@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = DummyFile()
    yield
    sys.stdout = save_stdout

#Ignore the warnings that Azimuth spits out
warnings.filterwarnings("ignore")
""" ************************************************ """

#global logger
logger = GlobalSettings.logger

# Class: scoring_window
# This class opens a window for the user to select which scoring algorithm to use
# Currently, the supported algorithms are:
#       * Azimuth 2.0, based on Rule Set 2 from Doench et al. 2016
#       * TBD...
# Inputs:
# * FASTA/FNA file that was used to create the CSPR DB that the gRNAs were generated from
#   - This is necessary because selected gRNAs must be found in the genomic sequence file
#     and algorithm-specific padding must be found for each gRNA. I.e. to get the Azimuth
#     score for a 23-nt gRNA+PAM sequence, 4 nt upstream of gRNA and 3 nt downstream of 
#     PAM must be included to get a final sequence length of 30-nt (4+20+3+3).
# Outputs:
# * Scores for the selected gRNAs based on the selected algorithm (updated in View Targets table)


class Scoring_Window(QtWidgets.QMainWindow):
    # init function. Sets all of the widgets/styling/buttons/variables.
    def __init__(self):
        try:
            # qt stuff
            super(Scoring_Window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/scoring_window.ui', self)
            self.progressBar.setValue(0) # Make sure the progress bar starts at 0

            # Button Connections
            self.browse_button.clicked.connect(self.fasta_browse)
            self.submit_button.clicked.connect(self.submit)

            # Initialize variables
            self.genome = "" # String to load in genome sequences. Key is chromosome number, value is sequence
            self.rev_genome = "" # String to load in reverse complement of genome sequences. Key is chromosome number, value is sequence

            # Styling
            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            margin-top: 10px;
                            font: bold 14pt 'Arial';}
                            """
            self.groupBox.setStyleSheet(groupbox_style)
            self.setWindowTitle("Select Scoring Algorithm")
            self.scaleUI()

            self.hide() # Start out hidden

        except Exception as e:
            logger.critical("Error initializing scoring_window class.")
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
            logger.critical("Error in scaleUI() in scoring_window.")
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
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()
            exit(-1)

    # browse for fasta function
    # allows user to browse for a fasta file where the full 30-nt sequence for each gRNA can be obtained 
    def fasta_browse(self):
        try:
            filed = QtWidgets.QFileDialog() # Initialize file dialog
            myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose a File")
            if (myFile[0] != ""): # Make sure file is not empty
                if not myFile[0].endswith(".fa") and not myFile[0].endswith(".fna") and not myFile[0].endswith(".fasta"): # Make sure file the correct type
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                    msgBox.setWindowTitle("File Selection Error")
                    msgBox.setText("You have selected an incorrect type of file. Please choose a FASTA/FNA file.")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()
                    return
                else:
                    self.file = myFile[0]
                    self.fasta_edit.setText(str(myFile[0]))

        except Exception as e:
            logger.critical("Error in fasta_browse() in scoring_window.")
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

    # This function loads in the selcted fasta/fna file's sequences so that they can be checked for gRNA sequences
    def load_fasta(self):
        try:
            if self.fasta_edit.text() == "": # If no file has been selected, throw an error and let user try again
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Error")
                msgBox.setText("No file has been selected. Please select a FASTA/FNA file.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()
            else:
                self.progressBar.setValue(10) # Update progress bar
                # Get the current chromosome to avoid having to load in the entire FASTA file to memory
                current_chrom = int(GlobalSettings.mainWindow.Results.featureDict[GlobalSettings.mainWindow.Results.comboBoxGene.currentText()][0])
                ### Save the chromosome of interest and its reverse complement to the strings self.genome and self.rev_genome, respectively.
                for i, record in enumerate(SeqIO.parse(self.fasta_edit.text(), "fasta")):
                    if i+1 == current_chrom: # Only save chromosome that contains the selected gRNAs 
                        self.genome = record.seq.__str__().upper()
                        self.rev_genome = record.seq.reverse_complement().__str__().upper()
                    else:
                        continue
                self.progressBar.setValue(30) # Update progress bar

        except Exception as e:
            logger.critical("Error loading in fasta file in scoring_window.")
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
    # This function searches for each selected gRNA in the FASTA file provided and if it is found, scores it according to the selected algorithm.
    # gRNAs that cannot be found are scored as -1's and reported to the user.
    def score_azimuth(self):
        try:
            ### Intialize the data structures for storing results
            guide_list = []
            reject_list = []
            full_seqs = []
            
            ### Find the full sequences for each gRNA from within the FASTA file and save them to full_seqs. Save index of sequences that could not be found.
            targettable = GlobalSettings.mainWindow.Results.targetTable # Make new variable for table name so that referencing isn't quite so tedious
            it = 0 # Initialize iterator

            self.progressBar.setValue(50) # Update progress bar
            scale = 20/len(targettable.selectionModel().selectedRows()) # Scale the next 20% of the progressBar to be completed over all the gRNAs

            for i in range(targettable.rowCount()):
                if targettable.item(i, 0).isSelected(): # If row was selected...
                    strand = targettable.item(i, 3).text() # Grab strand
                    seq = targettable.item(i, 2).text() + targettable.item(i,4).text() # Grab full seq (guide+PAM)
                    guide_list.append(seq) # Save the gRNA to guide_list
                    if strand == "+":
                        tmp = self.genome.find(seq) # Search forward strand of genome for gRNA
                        if tmp != -1: # If sequence was found
                            full_seqs.append(self.genome[tmp-4:tmp+26]) # Get full 30-nt sequence and append to full_seqsa
                        else: # If sequence was not found
                            reject_list.append(it) # Append index of gRNA to reject_list
                    else: # If strand == "-"
                        tmp = self.rev_genome.find(seq) # Search reverse strand of genome for gRNA
                        if tmp != -1: # If sequence was found
                            full_seqs.append(self.rev_genome[tmp-4:tmp+26]) # Get full 30-net sequence and append to pass list
                        else:
                            reject_list.append(it) # Append index of gRNA to reject_list 
                    increment = scale*(it+1) # Calculate progressBar increment for 1 gRNA search
                    self.progressBar.setValue(int(self.progressBar.value()+increment)) # Update progress bar
                    it += 1 # Iterate iterator
                else: # If row in table wasn't selected, skip it
                    continue 

            if len(full_seqs) != 0: # If some sequences were found...
                full_seqs = np.array(full_seqs) # Convert 30-nt seqs into a np.array
                with nostdout(): # Silence Azimuth's print statements
                    az_scores = az.predict(full_seqs)*100 # Predict scores and multiply by 100 (a la CHOP-CHOP)
                if len(reject_list) > 0: # If some sequences were not found
                    msg_string = "\n".join([guide_list[i] for i in reject_list]) # Separate missing sequences by newline char
                    ### Insert -1's for each guide that wasn't found
                    for i in reject_list:
                        az_scores = np.insert(az_scores, i, -1)
                    ### Generate a message showing the sequences that weren't found
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Information)
                    msgBox.setWindowTitle("Some Sequences Not Found!")
                    msgBox.setText("The following sequences (gRNA+PAM) were not found:\n%s\n\nThey have been given placeholder scores of -1." % msg_string)
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                    msgBox.exec()
            else: # If no sequences were found...
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Information)
                msgBox.setWindowTitle("No Sequences Found!")
                msgBox.setText("None of the selected guides were found. Please ensure you have the correct FASTA file and try again.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()
                return

            self.progressBar.setValue(80) # Update progress bar
            self.transfer_results(az_scores, "Azimuth 2.0") # Call transfer_results to fill the table with the scores

        except Exception as e:
            logger.critical("Error in score_azimuth() in scoring_window.")
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
    
    # This function searches for each selected gRNA in the FASTA file provided and if it is found, scores it according to the selected algorithm.
    # The table in View Targets is then updated to include the fresh scores. gRNAs that cannot be found are reported as -1's.
    def transfer_results(self, scores, algorithm):
        targettable = GlobalSettings.mainWindow.Results.targetTable # Make new variable for table name so that referencing isn't quite so tedious
        scale = 20/len(targettable.selectionModel().selectedRows()) # Scale the next 20% of progressBar to be completed over all the gRNAs
        headers = get_table_headers(targettable) # Get table headers
        num_cols = len(headers) # Get number of columns
        if algorithm not in headers: # If this algorithm hasn't been ran yet, add a column for it and place scores in it.
            targettable.insertColumn(num_cols-2) # Add column before off-target columns
            targettable.setHorizontalHeaderLabels(["Location", "Endonuclease", "Sequence", "Strand", "PAM", "Score", algorithm, "Off-Target", "Details"]) # Update column headers
            targettable.horizontalHeader().setSectionResizeMode(num_cols, QtWidgets.QHeaderView.Stretch) #Ensures last column goes to the edge of table

        col_index = get_table_headers(targettable).index(algorithm) # Get index of column to add scores to
        ### Add scores to algorithm column 
        it = 0
        for i in range(targettable.rowCount()):
            if targettable.item(i, 0).isSelected(): # If row was selected, add the corresponding score to the table
                score = QtWidgets.QTableWidgetItem() # Create new table item for the score to be placed into
                score.setData(QtCore.Qt.EditRole, round(float(scores[it]),2))
                targettable.setItem(i,col_index, score) # Add the item
                targettable.item(i,col_index).setSelected(True) # Make it selected by default to prevent export issues
                ## Update progress
                increment = scale*(it+1)
                self.progressBar.setValue(int(self.progressBar.value()+increment)) # Update progress bar
                it += 1


        targettable.resizeColumnsToContents()
        self.progressBar.setValue(100)
        self.close_window()
    
    def close_window(self):
        self.progressBar.setValue(0) # Reset progress bar
        self.hide()

    def submit(self):
        self.progressBar.setValue(0)
        self.load_fasta()
        self.score_azimuth()
