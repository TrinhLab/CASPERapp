from ast import Global
import GlobalSettings
import azimuth.model_comparison as az
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore, QtGui
import platform
import traceback
import math
from Bio import SeqIO
import numpy as np

#global logger
logger = GlobalSettings.logger

# Class: scoring_window
# This class opens a window for the user to select which scoring algorithm to use
# Currently, the supported algorithms are:
#       * Azimuth 2.0, based on Rule Set 2 from Doench et al. 2016
#       * TBD...
# Insert description of what this GUI does

class Scoring_Window(QtWidgets.QMainWindow):
    # init function. Sets all of the buttons
    def __init__(self):
        try:
            # qt stuff
            super(Scoring_Window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'scoring_window.ui', self)

            # button connections
            self.browse_button.clicked.connect(self.fasta_browse)
            self.submit_button.clicked.connect(self.submit)

            # Initialize variables
            self.genome = "" # String to load in genome sequences. Key is chromosome number, value is sequence
            self.rev_genome = "" # String to load in reverse complement of genome sequences. Key is chromosome number, value is sequence

            # GroupBox styling
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

            # variables
            # self.location = self.fileLocation_line_edit.text()
            # self.selected_table_items = []
            # self.window = ""
            # self.num_columns = []
            # self.locus_tag = False
            # self.gene_name = False

            self.setWindowTitle("Select Scoring Algorithm")
            self.scaleUI()
            self.hide()

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

    # launch function. Called in Results.
    # parameter expect: a list of the items selected from the window.
    def launch(self, select_items, window):
        try:
            if platform.system() == "Windows":
                self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "\\")
            else:
                self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "/")
            self.selected_table_items = select_items
            self.window = window
            self.centerUI()
            self.show()
            self.activateWindow()
        except Exception as e:
            logger.critical("Error in launch() in export to csv.")
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
                if self.window == "mt": ###Change headers for multitargeting table export
                    output_data = open(full_path, 'w')
                    header_string = 'Seed,Total Repeats,Avg. Repeats/Scaffold,Consensus Sequence,% Consensus,Full Sequence,Score,PAM,Strand\n'
                    output_data.write(header_string)
                elif self.window == "pa":
                    output_data = open(full_path, 'w')
                    header_string = 'Seed,% Coverage,Total Repeats,Avg. Repeats/Scaffold,Consensus Sequence,% Consensus,Full Sequence,Score,PAM,Strand\n'
                    output_data.write(header_string)
                else: ###Change headers for view results export
                    output_data = open(full_path, 'w')
                    if GlobalSettings.mainWindow.radioButton_Gene.isChecked(): # If the user chose to search via Feature
                        tmp = GlobalSettings.mainWindow.Results.comboBoxGene.currentText().split(":") # Check to see if the locus tag was found for the current gene
                        if len(tmp) > 1: # If locus tag exists for gene, include in output
                            header_string = 'Location,Endonuclease,Sequence,Full Sequence,Strand,PAM,Score,Off_Target,Locus_Tag,Gene_Name\n'
                            output_data.write(header_string)
                            self.locus_tag = True
                            self.gene_name = False
                        else: # If locus tag does not exist for gene, only include the gene name
                            header_string = 'Location,Endonuclease,Sequence,Full Sequence,Strand,PAM,Score,Off_Target,Gene_Name\n'
                            output_data.write(header_string)
                            self.gene_name = True
                            self.locus_tag = False
                    else: # If user searched by sequence or position, don't include locus tag or gene name
                        header_string = 'Location,Endonuclease,Sequence,Full Sequence,Strand,PAM,Score,Off_Target\n'
                        output_data.write(header_string)
                        self.gene_name = False
                        self.locus_tag = False

                # loop through and write the other data
                num_cols = len(header_string.split(",")) # Calculate the number of columns based on the header_string above
                if self.window == "vt" and self.locus_tag: #If the user is exporting data from VT and locus tag exists for current gene
                    tmp = GlobalSettings.mainWindow.Results.comboBoxGene.currentText().split(":") # Get the locus tag
                    locus_tag = str(tmp[0].strip())
                    gene_name = str(tmp[-1].strip())
                    # Get the gene name
                    i = 0 # Initialize iterator
                    for item in self.selected_table_items: # Loop through all the items in the View Targets table
                        if i == num_cols-3:
                            output_data.write(item.text()) #Off-target score
                            output_data.write(',')
                            output_data.write(locus_tag) #Locus tag
                            output_data.write(',')
                            output_data.write(gene_name) #Gene name
                            output_data.write('\n')
                            i = 0 # If reached the end of the row, reset the iterator
                        elif i == 2:
                            output_data.write(item.text()) #gRNA only
                            output_data.write(',')
                            output_data.write(self.leading_seq.text().strip() + item.text() + self.trailing_seq.text().strip()) #5' Leader + gRNA + 3' Trailer
                            output_data.write(',')
                            i += 2 # The plus 2 makes up for the extra column that is not in the table's items
                        else:
                            output_data.write(item.text())
                            output_data.write(',')
                            i += 1

                elif self.window == "vt" and self.gene_name: #If the user is exporting data from VT and locus tag doesn't exist for current gene
                    gene_name = str(GlobalSettings.mainWindow.Results.comboBoxGene.currentText().strip()) # Get the locus tag
                    i = 0 # Initialize iterator
                    for item in self.selected_table_items: # Loop through all the items in the View Targets table
                        if i == num_cols-2:
                            output_data.write(item.text())
                            output_data.write(',')
                            output_data.write(gene_name)
                            output_data.write('\n')
                            i = 0
                        elif i == 2:
                            output_data.write(item.text()) #gRNA only
                            output_data.write(',')
                            output_data.write(self.leading_seq.text().strip() + item.text() + self.trailing_seq.text().strip()) #5' Leader + gRNA + 3' Trailer
                            output_data.write(',')
                            i += 2 # The plus 2 makes up for the extra column that is not in the table's items
                        else:
                            output_data.write(item.text())
                            output_data.write(',')
                            i += 1
                elif self.window == "mt": #If the user is exporting data from multitargeting
                    i = 0 # Initialize iterator
                    for item in self.selected_table_items: # Loop through all the items in the View Targets table
                        if i == num_cols-1:
                            output_data.write(item.text())
                            output_data.write('\n')
                            i = 0
                        elif i == 3:
                            gRNA = item.text()
                            output_data.write(gRNA)
                            output_data.write(',')
                            i += 1
                        elif i == 4:
                            output_data.write(item.text()) #gRNA only
                            output_data.write(',')
                            output_data.write(self.leading_seq.text().strip() + gRNA + self.trailing_seq.text().strip()) #5' Leader + gRNA + 3' Trailer
                            output_data.write(',')
                            i += 2 # The plus 2 makes up for the extra column that is not in the table's items
                        else:
                            output_data.write(item.text())
                            output_data.write(',')
                            i += 1

                elif self.window == "pa": #If the user is exporting data from Population Analysis
                    i = 0 # Initialize iterator
                    for item in self.selected_table_items: # Loop through all the items in the View Targets table
                        if i == num_cols-1:
                            output_data.write(item.text())
                            output_data.write('\n')
                            i = 0
                        elif i == 4:
                            gRNA = item.text()
                            output_data.write(gRNA)
                            output_data.write(',')
                            i += 1
                        elif i == 5:
                            output_data.write(item.text()) #gRNA only
                            output_data.write(',')
                            output_data.write(self.leading_seq.text().strip() + gRNA + self.trailing_seq.text().strip()) #5' Leader + gRNA + 3' Trailer
                            output_data.write(',')
                            i += 2 # The plus 2 makes up for the extra column that is not in the table's items
                        else:
                            output_data.write(item.text())
                            output_data.write(',')
                            i += 1

                else: #If the user is exporting data from View Targets but is not using Feature search
                    i = 0 # Initialize iterator
                    for item in self.selected_table_items: # Loop through all the items in the View Targets table
                        if i == num_cols-1:
                            output_data.write(item.text())
                            output_data.write('\n')
                            i = 0
                        elif i == 2:
                            output_data.write(item.text()) #gRNA only
                            output_data.write(',')
                            output_data.write(self.leading_seq.text().strip() + item.text() + self.trailing_seq.text().strip()) #5' Leader + gRNA + 3' Trailer
                            output_data.write(',')
                            i += 2 # The plus 2 makes up for the extra column that is not in the table's items
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
            filed = QtWidgets.QFileDialog()
            myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose a File")
            if (myFile[0] != ""):
                if not myFile[0].endswith(".fa") and not myFile[0].endswith(".fna") and not myFile[0].endswith(".fasta"):
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

    # This function loads in the selcted fasta/fna file so that it can be checked for gRNA sequences
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
                current_chrom = int(GlobalSettings.mainWindow.Results.featureDict[GlobalSettings.mainWindow.Results.comboBoxGene.currentText()][0])
                # Save the chromosome of interest and its reverse complement to the dictionaries self.genome and self.rev_genome, respectively.
                for i, record in enumerate(SeqIO.parse(self.fasta_edit.text(), "fasta")):
                    if i+1 == current_chrom: # Only save chromosome the that contain the selected gRNAs 
                        self.genome = record.seq.__str__()
                        self.rev_genome = record.seq.reverse_complement().__str__()
                    else:
                        continue

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





    def score_sequences(self):
        try:
            ### Intialize the data structures for storing results
            score_dict = {}
            pass_list = []
            reject_list = []
            full_seqs = []
            azimuth_scores = []
            targettable = GlobalSettings.mainWindow.Results.targetTable # Make new variable for table name so that referencing isn't quite so tedious

            for i in range(targettable.rowCount()):
                if targettable.item(i, 0).isSelected(): # If row was selected...
                    strand = targettable.item(i, 3).text() 
                    seq = targettable.item(i, 2).text() + targettable.item(i,4).text()
                    if strand == "+":
                        tmp = self.genome.find(seq)
                        if tmp != 1:
                            pass_list.append(seq) # Save the 
                            full_seqs.append(self.genome[tmp-4:tmp+26]) # Get full 30-net sequence and append to pass list
                        else:
                            reject_list.append(seq) # This gRNA was not able to be found
                            azimuth_scores.append(-1) # Return a -1 for the score
                    else:
                        tmp = self.rev_genome.find(seq)
                        if tmp != 1:
                            pass_list.append(seq) # Save the 
                            full_seqs.append(self.rev_genome[tmp-4:tmp+26]) # Get full 30-net sequence and append to pass list
                        else:
                            reject_list.append(seq) # This gRNA was not able to be found
                            azimuth_scores.append(-1) # Return a -1 for the score
                else: # If row wasn't selected, skip it
                    continue 
            full_seqs = np.array(full_seqs) # Convert 30-nt seqs into a np.array
            az_scores = az.predict(full_seqs) # Predict scores
            print(len(az_scores))

            """
            TO DO
            1. Finish connecting scores back to Results table
            2. Add Error checking/edge case handling for when gRNAs cannot be found
            3. Add loading bar and text label progress output to scoring window
            4. Clean all the code up and make sure everything is neatly commented
            """

            # if any of the sequences return 0 matches, show the user which ones were not found
            # if len(failed_guides) > 0:
            #     msgBox = QtWidgets.QMessageBox()
            #     msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            #     msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            #     msgBox.setWindowTitle("Warning")
            #     msgBox.setText(
            #         "The following sequence(s) were not found in the Gene Viewer text:\n\t" + "\n".join(failed_guides))
            #     msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
            #     msgBox.exec()

        except Exception as e:
            logger.critical("Error in score_sequences() in scoring_window.")
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

    def submit(self):
        self.load_fasta()
        self.score_sequences()
        # self.close_window()
