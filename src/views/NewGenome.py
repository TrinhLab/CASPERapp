from ast import Global
import os
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
import models.GlobalSettings as GlobalSettings
from functools import partial
from utils.Algorithms import SeqTranslate
import webbrowser
import platform
import traceback
import math
from utils.ui import show_message, show_error, scale_ui, center_ui
from utils.web import ncbi_page, repo_page

logger = GlobalSettings.logger

def iter_except(function, exception):
    """Works like builtin 2-argument `iter()`, but stops on `exception`."""
    try:
        while True:
            yield function()
    except exception:
        return

#UI prompt for when the user has finished running jobs in new genome to allow them to choose where the want to proceed
class goToPrompt(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(goToPrompt, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/newgenomenavigationpage.ui', self)

            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            font: bold 14pt 'Arial';
                            margin-top: 10px;}"""
            self.groupBox.setStyleSheet(groupbox_style)
            scale_ui(self, custom_scale_width=575, custom_scale_height=175)
            self.setWindowTitle("New Genome")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.hide()

        except Exception as e:
            show_error("Unable to initialize goToPrompt class in New Genome.", e)

#New genome class to allow users to generate new CSPR files
class NewGenome(QtWidgets.QMainWindow):
    def __init__(self, info_path):
        try:
            super(NewGenome, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/NewGenome.ui', self)
            self.setWindowTitle('New Genome')
            self.setWindowTitle('New Genome')
            self.info_path = info_path

            #---Style Modifications---#

            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            font: bold 14pt 'Arial';
                            margin-top: 10px;}"""

            self.Step1.setStyleSheet(groupbox_style)
            self.Step2.setStyleSheet(groupbox_style.replace("Step1","Step2"))
            self.Step3.setStyleSheet(groupbox_style.replace("Step1","Step3"))

            #---Button Modifications---#

            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.resetButton.clicked.connect(self.reset)
            self.submitButton.clicked.connect(self.submit)
            self.browseForFile.clicked.connect(self.selectFasta)
            self.remove_job.clicked.connect(self.remove_from_queue)
            self.output_browser.setText("Waiting for program initiation...")
            self.contButton.clicked.connect(self.continue_to_main)

            self.comboBoxEndo.currentIndexChanged.connect(self.endo_settings)

            self.runButton.clicked.connect(self.run_jobs_wrapper)
            self.clearButton.clicked.connect(self.clear_all)

            self.JobsQueue = []  # holds Job classes.
            self.check_strings = []
            self.Endos = dict()
            self.file = ""

            self.process = QtCore.QProcess()
            self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
            self.process.finished.connect(self.upon_process_finishing)
            self.seqTrans = SeqTranslate()
            self.exit = False

            self.first = False
            #show functionalities on window
            self.fillEndo()

            self.num_chromo_next = False

            #Jobs Table
            self.job_Table.setShowGrid(False)
            self.job_Table.horizontalHeader().setSectionsClickable(True)
            self.job_Table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.job_Table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.job_Table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.job_Table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
            self.fin_index=0

            self.mwfg = self.frameGeometry()  ##Center window
            self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
            self.total_chrom_count = 0
            self.perc_increase = 0
            self.progress = 0

            #toolbar button actions
            self.visit_repo.triggered.connect(repo_page)
            self.go_ncbi.triggered.connect(ncbi_page)

            self.comboBoxEndo.currentIndexChanged.connect(self.changeEndos)

            ### NCBI tool
            self.NCBI_File_Search.clicked.connect(self.open_ncbi_tool)

            self.seed_length.setEnabled(False)
            self.five_length.setEnabled(False)
            self.three_length.setEnabled(False)
            self.repeats_box.setEnabled(False)

            ### User prompt class
            self.goToPrompt = goToPrompt()
            self.goToPrompt.goToMain.clicked.connect(self.continue_to_main)
            self.goToPrompt.goToMT.clicked.connect(self.continue_to_MT)
            self.goToPrompt.goToPop.clicked.connect(self.continue_to_pop)

            self.orgName.setFocus()

            ### Connect New endonuclease to New Genome
            self.actionUpload_New_Endonuclease.triggered.connect(self.launch_newEndonuclease)

            ### Set up validators for input fields:
            reg_ex1 = QtCore.QRegExp("[^/\\\\_]+") # No slashes or underscores
            reg_ex2 = QtCore.QRegExp("\\S+")
            input_validator1 = QtGui.QRegExpValidator(reg_ex1, self)
            input_validator2 = QtGui.QRegExpValidator(reg_ex2, self)
            self.orgName.setValidator(input_validator1)
            self.strainName.setValidator(input_validator1)
            self.orgCode.setValidator(input_validator2)

            scale_ui(self, custom_scale_width=850, custom_scale_height=750)
            self.first_show = True
        except Exception as e:
            show_error("Error initializing New Genome class.", e)

    def launch_newEndonuclease(self):
        try:
            GlobalSettings.mainWindow.getData()
            GlobalSettings.mainWindow.newEndonuclease.centerUI()
            GlobalSettings.mainWindow.newEndonuclease.show()
            GlobalSettings.mainWindow.newEndonuclease.activateWindow()
        except Exception as e:
            show_error("Error in launch_newEndonuclease() in New Genome.", e)

    #open the ncbi search tool window
    def open_ncbi_tool(self):
        try:
            #center ncbi on current screen
            if GlobalSettings.mainWindow.ncbi.first_show == True:
                GlobalSettings.mainWindow.ncbi.first_show = False
                GlobalSettings.mainWindow.ncbi.centerUI()
            if self.orgName.text() != "":
                GlobalSettings.mainWindow.ncbi.organism_line_edit.setText(self.orgName.text())
            if self.strainName.text() != "":
                GlobalSettings.mainWindow.ncbi.infra_name_line_edit.setText(self.strainName.text())
            GlobalSettings.mainWindow.ncbi.show()
            GlobalSettings.mainWindow.ncbi.activateWindow()
        except Exception as e:
            show_error("Error in open_ncbi_tool() in New Genome.", e)

    def remove_from_queue(self):
        try:
            while(True):
                indexes = self.job_Table.selectionModel().selectedRows()
                if len(indexes) == 0:
                    break
                self.job_Table.removeRow(indexes[0].row())
        except Exception as e:
            show_error("Error in remove_from_queue() in New Genome.", e)

    #prompt user with file browser to select fasta/fna files
    def selectFasta(self):
        try:
            filed = QtWidgets.QFileDialog()
            myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose a File")
            if (myFile[0] != ""):
                if not myFile[0].endswith(".fa") and not myFile[0].endswith(".fna") and not myFile[0].endswith(".fasta"):
                    show_message(
                        fontSize=12,
                        icon=QtWidgets.QMessageBox.Icon.Critical,
                        title="File Selection Error",
                        message="You have selected an incorrect type of file. Please choose a FASTA/FNA file."
                    )
                    return
                else:
                    self.file = myFile[0]
                    self.selectedFile.setText(str(myFile[0]))
        except Exception as e:
            show_error("Error in selectFasta() in New Genome.", e)

    #submit jobs to queue
    def submit(self):
        try:
            warning = ""
            if len(self.orgName.text()) == 0:
                warning = warning + "You need to include the organism's name."
            if len(self.file) == 0:
                warning = warning + "You need to select a file."
            if len(warning) != 0:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Required Information",
                    message=warning
                )
                return
            if len(self.strainName.text()) == 0:
                warning = warning + "\nIt is recommended to include the organism's subspecies/strain."
            if len(self.orgCode.text()) == 0:
                warning = warning + "\nYou must include an organism code."
            if len(warning) != 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                msgBox.setWindowTitle("Missing Information")
                msgBox.setText(warning + "\n\nDo you wish to continue without including this information?")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
                msgBox.exec()

                if msgBox.result() == QtWidgets.QMessageBox.No:
                    return

            #endo, pam, repeats, directionality, five length, seed length, three length, orgcode, output path, CASPERinfo path, fna path, orgName, notes, on target matrix
            args = self.Endos[self.comboBoxEndo.currentText()][0]
            args += " " + self.Endos[self.comboBoxEndo.currentText()][1]
            if self.mt.isChecked():
                args += " " + "TRUE"
            else:
                args += " " + "FALSE"

            if self.Endos[self.comboBoxEndo.currentText()][5] == "3":
                args += " " + "FALSE"
            else:
                args += " " + "TRUE"

            if self.repeats_box.isChecked():
                args += " " + "TRUE"
            else:
                args += " " + "FALSE"

            args += " " + self.Endos[self.comboBoxEndo.currentText()][2]
            args += " " + self.Endos[self.comboBoxEndo.currentText()][3]
            args += " " + self.Endos[self.comboBoxEndo.currentText()][4]
            args += " " + self.orgCode.text()
            if platform.system() == 'Windows':
                args += " " + '"' + GlobalSettings.CSPR_DB.replace("/","\\") + '\\"'
                args += " " + '"' + GlobalSettings.appdir.replace("/","\\") + "CASPERinfo" + '"'
                args += " " + '"' + self.file.replace("/","\\") + '"'
            else:
                args += " " + '"' + GlobalSettings.CSPR_DB.replace("\\","/") + '/"'
                args += " " + '"' + GlobalSettings.appdir.replace("\\","/") + "CASPERinfo" + '"'
                args += " " + '"' + self.file.replace("\\","/") + '"'

            args += " " + '"' + self.orgName.text() + " " + self.strainName.text() + '"'
            args += " " + '"' + "notes" + '"'
            args += " " + '"DATA:' + self.Endos[self.comboBoxEndo.currentText()][6] + '"'

            tmp = self.orgName.text()+ " " + self.strainName.text() + " " + self.Endos[self.comboBoxEndo.currentText()][0] + " " + self.orgCode.text()
            if tmp in self.check_strings:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Duplicate Entry",
                    message="You have submitted a duplicate entry. Consider changing the organism code or strain name to differentiate closely related strains."
                )
                return
            name = self.orgCode.text() + "_" + str(self.Endos[self.comboBoxEndo.currentText()][0])
            rowPosition = self.job_Table.rowCount()
            self.job_Table.insertRow(rowPosition)
            item = QtWidgets.QTableWidgetItem(name)
            item.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.job_Table.setItem(rowPosition, 0, item)
            self.check_strings.append(tmp)
            self.JobsQueue.append(args)
        except Exception as e:
            show_error("Error in submit() in New Genome.", e)

    #fill the endo dropdown
    def fillEndo(self):
        try:
            #disconnect signal
            try:
                self.comboBoxEndo.currentIndexChanged.disconnect()
            except:
                pass

            #clear out the endo box
            self.comboBoxEndo.clear()

            f = open(GlobalSettings.appdir + "CASPERinfo")
            while True:
                line = f.readline()
                if line.startswith('ENDONUCLEASES'):
                    while True:
                        line = f.readline()
                        if (line[0] == "-"):
                            break
                        line_tokened = line.split(";")
                        if len(line_tokened) == 10:
                            endo = line_tokened[0]
                            # Checking to see if there is more than one pam sequence in the list
                            if line_tokened[1].find(",") != -1:
                                p_pam = line_tokened[1].split(",")[0]
                            else:
                                p_pam = line_tokened[1]
                            five_length = line_tokened[2]
                            seed_length = line_tokened[3]
                            three_length = line_tokened[4]
                            dir = line_tokened[5]
                            on_target_data = line_tokened[8]
                            self.Endos[endo + " - PAM: " + p_pam] = (endo, p_pam, five_length, seed_length, three_length, dir, on_target_data)
                    break
            f.close()
            self.comboBoxEndo.addItems(self.Endos.keys())
            key = list(self.Endos.keys())[0]
            self.seed_length.setText(self.Endos[key][3])
            self.five_length.setText(self.Endos[key][2])
            self.three_length.setText(self.Endos[key][4])

            #reconnect signal
            self.comboBoxEndo.currentIndexChanged.connect(self.changeEndos)
        except Exception as e:
            show_error("Error in fillEndo() in New Genome.", e)

    #event handler for endo changing - update endo length data
    def changeEndos(self):
        try:
            key = str(self.comboBoxEndo.currentText())
            self.seed_length.setText(self.Endos[key][3])
            self.five_length.setText(self.Endos[key][2])
            self.three_length.setText(self.Endos[key][4])
        except Exception as e:
            show_error("Error in changeEndos() in New Genome.", e)

    #check if endo is 3' or 5'
    def endo_settings(self):
        try:
            # check the if it's 3' or 5', and check the box accordingly
            if int(self.seqTrans.endo_info[self.Endos[self.comboBoxEndo.currentText()][0]][3]) == 3:
                self.pamBox.setChecked(0)
            elif int(self.seqTrans.endo_info[self.Endos[self.comboBoxEndo.currentText()][0]][3]) == 5:
                self.pamBox.setChecked(1)
        except Exception as e:
            show_error("Error in endo_settings() in New Genome.", e)

    #wrapper for running jobs
    def run_jobs_wrapper(self):
        try:
            self.indexes = []
            self.job_Table.selectAll()
            indexes = self.job_Table.selectionModel().selectedRows()
            for index in sorted(indexes):
                if self.job_Table.item(index.row(), 0).text() != "":
                    self.indexes.append(index.row())
            self.run_job()
        except Exception as e:
            show_error("Error in run_jobs_wrapper() in New Genome.", e)

    #run job in queue
    def run_job(self):
        try:
            if len(self.indexes) > 0:
                self.progressBar.setValue(0)
                self.progress = 0
                row_index = self.indexes[0]
                name = self.job_Table.item(row_index, 0).text()
                item = QtWidgets.QTableWidgetItem(name)
                item.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.job_Table.setItem(row_index, 1, item)
                self.job_Table.setItem(row_index, 0, QtWidgets.QTableWidgetItem(""))

                def output_stdout(p):
                    line = str(p.readAll())
                    line = line[2:]
                    line = line[:len(line) - 1]
                    for lines in line.split(r"\n"):
                        lines = lines.rstrip("\n")
                        lines = lines.rstrip("\r")
                        lines = lines.rstrip(r"\n")
                        lines = lines.rstrip(r"\r")
                        lines = lines.rstrip("\r\n")
                        lines = lines.rstrip(r"\r\n")
                        if lines != "":
                            if lines.find("Number of Chromosomes/Scaffolds") != -1:
                                copy = lines
                                copy = copy.replace(" ","")
                                copy = copy[copy.find(":")+1:]
                                self.total_chrom_count = int(copy)
                                self.perc_increase = ((1 / (2 * self.total_chrom_count)) * 70)
                                self.progressBar.setValue(20)
                                self.progress = 20
                            elif lines.find("complete.") != -1:
                                self.progress += self.perc_increase
                                self.progressBar.setValue(int(self.progress))
                            elif lines.find("Processing Targets.") != -1:
                                self.progress = 70
                                self.progressBar.setValue(int(self.progress))
                            elif lines.find("Writing out uniques.") != -1:
                                self.progress = 90
                                self.progressBar.setValue(int(self.progress))
                            elif lines.find("Writing out repeats.") != -1:
                                self.progress = 95
                                self.progressBar.setValue(int(self.progress))
                            elif lines == "Finished.":
                                self.progress = 100
                                self.progressBar.setValue(int(self.progress))
                            self.output_browser.append(lines)

                job_args = self.JobsQueue[row_index]
                if platform.system() == 'Windows':
                    program = '"' + GlobalSettings.appdir + "SeqFinderFolder/Casper_Seq_Finder_Win.exe" + '" '
                elif platform.system() == 'Linux':
                    program = '"' + GlobalSettings.appdir + "SeqFinderFolder/Casper_Seq_Finder_Lin" + '" '
                else:
                    program = '"' + GlobalSettings.appdir + "SeqFinderFolder/Casper_Seq_Finder_Mac" + '" '
                program += job_args
                self.process.readyReadStandardOutput.connect(partial(output_stdout, self.process))
                self.process.start(program)
            else:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="No Jobs To Run",
                    message="No jobs are in the queue to run. Please add a job before running."
                )
        except Exception as e:
            show_error("Error in run_job() in New Genome.", e)

    #even handler for when jobs finish execution
    def upon_process_finishing(self):
        try:
            row_index = self.indexes[0]
            name = self.job_Table.item(row_index, 1).text()
            item = QtWidgets.QTableWidgetItem(name)
            item.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.job_Table.setItem(row_index, 2, item)
            self.job_Table.setItem(row_index, 1, QtWidgets.QTableWidgetItem(""))
            self.indexes.pop(0)
            if len(self.indexes) != 0:
                self.run_job()
            else:
                #prompt user if they want to analyze their new files
                center_ui(self.goToPrompt)
                self.goToPrompt.show()
                self.goToPrompt.activateWindow()
        except Exception as e:
            show_error("Error in upon_process_finishing() in New Genome.", e)

    #clear the job table
    def clear_all(self):
        try:
            self.process.kill()
            self.fin_index = 0
            self.job_Table.clearContents()
            self.job_Table.setRowCount(0)
            self.JobsQueue = []
            self.check_strings = []
            self.output_browser.clear()
            self.output_browser.setText("Waiting for program initiation...")
            self.orgName.clear()
            self.strainName.clear()
            self.orgCode.clear()
            self.selectedFile.clear()
            self.selectedFile.setPlaceholderText("Selected FASTA/FNA File")
            self.progressBar.setValue(0)
            self.first = False
        except Exception as e:
            show_error("Error in clear_all() in New Genome.", e)

    #reset the whole form
    def reset(self):
        try:
            self.orgName.clear()
            self.strainName.clear()
            self.orgCode.clear()
            self.selectedFile.clear()
            self.selectedFile.setPlaceholderText("Selected FASTA/FNA File")
            self.output_browser.clear()
            self.output_browser.setText("Waiting for program initiation...")
            self.file = ""
        except Exception as e:
            show_error("Error in reset() in New Genome.", e)
    
    #event handler for user wanting to close the window
    def closeEvent(self, event):
        try:
            # make sure that there are cspr files in the DB
            file_names = os.listdir(GlobalSettings.CSPR_DB)
            noCSPRFiles = True
            for file in file_names:
                if 'cspr' in file:
                    noCSPRFiles = False
                    break
            if noCSPRFiles == True:
                if self.exit == False:
                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                    msgBox.setWindowTitle("No CSPR file generated")
                    msgBox.setText("No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                                                            "Alternatively, you could quit the program. Would you like to quit?")
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                    msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
                    msgBox.exec()

                    if (msgBox.result() == QtWidgets.QMessageBox.No):
                        event.ignore()
                    else:
                        event.accept()
                else:
                    self.exit = False
                    event.accept()
            else:
                self.process.kill()
                self.clear_all()
                self.goToPrompt.hide()
                GlobalSettings.mainWindow.fill_annotation_dropdown()
                if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                    GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
                GlobalSettings.mainWindow.orgChoice.clear()
                GlobalSettings.mainWindow.endoChoice.clear()
                GlobalSettings.mainWindow.getData()
                GlobalSettings.MTWin.launch()
                GlobalSettings.pop_Analysis.launch()

                if GlobalSettings.mainWindow.first_show == True:
                    GlobalSettings.mainWindow.first_show = False
                    GlobalSettings.mainWindow.centerUI()
                GlobalSettings.mainWindow.show()
                event.accept()
        except Exception as e:
            show_error("Error in closeEvent() in New Genome.", e)

    #event handler for user wanting to go to Main once jobs complete
    def continue_to_main(self):
        try:
            # make sure that there are cspr files in the DB
            file_names = os.listdir(GlobalSettings.CSPR_DB)
            noCSPRFiles = True
            for file in file_names:
                if 'cspr' in file:
                    noCSPRFiles = False
                    break
            if noCSPRFiles == True:

                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                msgBox.setWindowTitle("No CSPR file generated")
                msgBox.setText(
                    "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                    "Alternatively, you could quit the program. Would you like to quit?")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
                msgBox.exec()

                if (msgBox.result() == QtWidgets.QMessageBox.Yes):
                    self.exit = True
                    self.close()
            else:
                self.process.kill()
                self.clear_all()
                self.goToPrompt.hide()
                GlobalSettings.mainWindow.fill_annotation_dropdown()
                if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                    GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
                GlobalSettings.mainWindow.orgChoice.clear()
                GlobalSettings.mainWindow.endoChoice.clear()
                GlobalSettings.mainWindow.getData()
                GlobalSettings.MTWin.launch()
                GlobalSettings.pop_Analysis.launch()

                # center main on current screen
                if GlobalSettings.mainWindow.first_show == True:
                    GlobalSettings.mainWindow.first_show = False
                    center_ui(GlobalSettings.mainWindow)
                GlobalSettings.mainWindow.show()
                self.hide()
        except Exception as e:
            show_error("Error in continue_to_main() in New Genome.", e)

    #event handler for user wanting to go to multi-targeting once jobs complete
    def continue_to_MT(self):
        try:
            # make sure that there are cspr files in the DB
            file_names = os.listdir(GlobalSettings.CSPR_DB)
            noCSPRFiles = True
            for file in file_names:
                if 'cspr' in file:
                    noCSPRFiles = False
                    break
            if noCSPRFiles == True:

                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                msgBox.setWindowTitle("No CSPR file generated")
                msgBox.setText(
                    "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                    "Alternatively, you could quit the program. Would you like to quit?")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
                msgBox.exec()



                if (msgBox.result() == QtWidgets.QMessageBox.Yes):
                    self.exit = True
                    self.close()

            else:
                self.process.kill()
                self.clear_all()
                self.goToPrompt.hide()
                GlobalSettings.mainWindow.fill_annotation_dropdown()
                if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                    GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
                GlobalSettings.mainWindow.orgChoice.clear()
                GlobalSettings.mainWindow.endoChoice.clear()
                GlobalSettings.mainWindow.getData()
                GlobalSettings.MTWin.launch()
                GlobalSettings.pop_Analysis.launch()

                # center multi-targeting on current screen
                if GlobalSettings.MTWin.first_show == True:
                    GlobalSettings.MTWin.first_show = False
                    GlobalSettings.MTWin.centerUI()

                GlobalSettings.MTWin.show()
                self.hide()
        except Exception as e:
            show_error("Error in continue_to_MT() in New Genome.", e)

    #event handler for user wanting to go to population analysis once jobs complete
    def continue_to_pop(self):
        try:
            # make sure that there are cspr files in the DB
            file_names = os.listdir(GlobalSettings.CSPR_DB)
            noCSPRFiles = True
            for file in file_names:
                if 'cspr' in file:
                    noCSPRFiles = False
                    break
            if noCSPRFiles == True:

                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
                msgBox.setWindowTitle("No CSPR file generated")
                msgBox.setText(
                    "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                    "Alternatively, you could quit the program. Would you like to quit?")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
                msgBox.exec()

                if (msgBox.result() == QtWidgets.QMessageBox.Yes):
                    self.exit = True
                    self.close()

            else:
                self.process.kill()
                self.clear_all()
                self.goToPrompt.hide()
                GlobalSettings.mainWindow.fill_annotation_dropdown()
                if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                    GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
                GlobalSettings.mainWindow.orgChoice.clear()
                GlobalSettings.mainWindow.endoChoice.clear()
                GlobalSettings.mainWindow.getData()
                GlobalSettings.MTWin.launch()
                GlobalSettings.pop_Analysis.launch()

                if GlobalSettings.pop_Analysis.first_show == True:
                    GlobalSettings.pop_Analysis.first_show = False
                    GlobalSettings.pop_Analysis.centerUI()

                GlobalSettings.pop_Analysis.show()
                self.hide()
        except Exception as e:
            show_error("Error in continue_to_pop() in New Genome.", e)