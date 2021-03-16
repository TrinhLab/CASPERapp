import sys, os
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
import GlobalSettings
from functools import partial
from Algorithms import SeqTranslate
from PyQt5.QtWidgets import *
import webbrowser
import platform

def iter_except(function, exception):
    """Works like builtin 2-argument `iter()`, but stops on `exception`."""
    try:
        while True:
            yield function()
    except exception:
        return

class NewGenome(QtWidgets.QMainWindow):


    def __init__(self, info_path):
        super(NewGenome, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'NewGenome.ui', self)
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
                        font: 15pt "Sans Serif";
                        font: bold;
                        margin-top: 10px;}"""

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1","Step2").replace("rgb(111,181,110)","rgb(77,158,89)"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1","Step3").replace("rgb(111,181,110)","rgb(53,121,93)"))

        #---Button Modifications---#

        self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.png"))
        self.whatsthisButton.clicked.connect(self.whatsthisclicked)
        self.resetButton.clicked.connect(self.reset)
        self.submitButton.clicked.connect(self.submit)
        self.browseForFile.clicked.connect(self.selectFasta)
        self.remove_job.clicked.connect(self.remove_from_queue)
        self.output_browser.setText("Waiting for program initiation...")
        self.contButton.clicked.connect(self.continue_to_main)

        self.comboBoxEndo.currentIndexChanged.connect(self.endo_settings)

        self.runButton.clicked.connect(self.run_jobs_wrapper)
        self.clearButton.clicked.connect(self.clear_job_queue)

        self.viewStatButton.setEnabled(False)

        self.JobsQueue = []  # holds Job classes.
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
        #self.show()

        self.num_chromo_next = False


        #Jobs Table
        #self.job_Table.setColumnCount(3)
        self.job_Table.setShowGrid(False)
        #self.job_Table.setHorizontalHeaderLabels(["Job Queue","Job in Progress", "Completed Jobs"])
        self.job_Table.horizontalHeader().setSectionsClickable(True)
        #self.job_Table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        #self.job_Table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
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
        self.visit_repo.triggered.connect(self.visit_repo_func)
        self.go_ncbi.triggered.connect(self.open_ncbi_web_page)

        self.comboBoxEndo.currentIndexChanged.connect(self.changeEndos)

        #some defaults pams for testing
        # self.orgName.setText("Bacillus Coagulans")
        # self.orgCode.setText("testing")

        self.seed_length.setEnabled(False)
        self.five_length.setEnabled(False)
        self.three_length.setEnabled(False)
        self.repeats_box.setEnabled(False)

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def remove_from_queue(self):
        while(True):
            indexes = self.job_Table.selectionModel().selectedRows()
            if len(indexes) == 0:
                break
            self.job_Table.removeRow(indexes[0].row())


    def selectFasta(self):

        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose a File")
        if (myFile[0] != ""):

            if not myFile[0].endswith(".fa") and not myFile[0].endswith(".fna") and not myFile[0].endswith(".gbff") and not myFile[0].endswith(".fasta"):
                QtWidgets.QMessageBox.question(self, "File Selection Error",
                                               "You have selected an incorrect type of file. "
                                               "Please choose a FASTA/FNA file.",
                                               QtWidgets.QMessageBox.Ok)
                return
            else:
                self.file = myFile[0]
                self.s_file.setText(str(myFile[0]))
        """cdir = self.lineEdit.text()
        os.chdir(mydir)
        self.gdirectory = mydir
        print(mydir)
        print(cdir)"""


    def submit(self):
        warning = ""
        if len(self.orgName.text()) == 0:
            warning = warning + "\nYou need to include the organism's name."
        if len(self.file) == 0:
            warning = warning + "\nYou need to select a file."
        if len(warning) != 0:
            QtWidgets.QMessageBox.information(self, "Required Information", warning, QtWidgets.QMessageBox.Ok)
            return

        if len(self.strainName.text()) == 0:
            warning = warning + "\nIt is recommended to include the organism's subspecies/strain."
        if len(self.orgCode.text()) == 0:
            warning = warning + "\nYou must include an organism code."
        if len(warning) != 0:
            hold = QtWidgets.QMessageBox.question(self, "Missing Information", warning +
                                                  "\n\nDo you wish to continue without including this information?"
                                                  , QtWidgets.QMessageBox.Yes |
                                                  QtWidgets.QMessageBox.No,
                                                  QtWidgets.QMessageBox.No)
            if hold == QtWidgets.QMessageBox.No:
                return



        #endo, pam, repeats, directionality, five length, seed length, three length, orgcode, output path, CASPERinfo path, fna path, orgName, notes
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

        args += " " + '"' + self.orgName.text() + '"'
        args += " " + '"' + "notes" + '"'

        name = self.orgCode.text() + "_" + str(self.Endos[self.comboBoxEndo.currentText()][0])
        rowPosition = self.job_Table.rowCount()
        self.job_Table.insertRow(rowPosition)
        item = QtWidgets.QTableWidgetItem(name)
        item.setTextAlignment(QtCore.Qt.AlignHCenter)
        self.job_Table.setItem(rowPosition, 0, item)
        self.JobsQueue.append(args)


    def fillEndo(self):
        f = open(GlobalSettings.appdir + "CASPERinfo")

        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if (line[0] == "-"):
                        break
                    line_tokened = line.split(";")
                    if len(line_tokened) == 8:
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
                        self.Endos[endo + " - PAM: " + p_pam] = (endo, p_pam, five_length, seed_length, three_length, dir)
                break
        f.close()
        self.comboBoxEndo.addItems(self.Endos.keys())
        key = list(self.Endos.keys())[0]
        self.seed_length.setText(self.Endos[key][3])
        self.five_length.setText(self.Endos[key][2])
        self.three_length.setText(self.Endos[key][4])


    def changeEndos(self):
        key = str(self.comboBoxEndo.currentText())
        self.seed_length.setText(self.Endos[key][3])
        self.five_length.setText(self.Endos[key][2])
        self.three_length.setText(self.Endos[key][4])


    def endo_settings(self):
        # check the if it's 3' or 5', and check the box accordingly
        if int(self.seqTrans.endo_info[self.Endos[self.comboBoxEndo.currentText()][0]][3]) == 3:
            self.pamBox.setChecked(0)
        elif int(self.seqTrans.endo_info[self.Endos[self.comboBoxEndo.currentText()][0]][3]) == 5:
            self.pamBox.setChecked(1)


    def findFasta(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to quit?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            sys.exit()
        else:
            pass


    def whatsthisclicked(self):
        QtWidgets.QMessageBox.information(self, "Organism Code", "The organism code is the manner in which CASPER will"
                                                                 " label its data files and references for the organism"
                                                                 " you are importing here.", QtWidgets.QMessageBox.Ok)


    def run_jobs_wrapper(self):
        self.indexes = []
        indexes = self.job_Table.selectionModel().selectedRows()
        for index in sorted(indexes):
            if self.job_Table.item(index.row(), 0).text() != "":
                self.indexes.append(index.row())
        self.run_job()


    def run_job(self):
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
                            self.progressBar.setValue(self.progress)
                        elif lines.find("Processing Targets.") != -1:
                            self.progress = 70
                            self.progressBar.setValue(self.progress)
                        elif lines.find("Writing out uniques.") != -1:
                            self.progress = 90
                            self.progressBar.setValue(self.progress)
                        elif lines.find("Writing out repeats.") != -1:
                            self.progress = 95
                            self.progressBar.setValue(self.progress)
                        elif lines == "Finished.":
                            self.progress = 100
                            self.progressBar.setValue(self.progress)
                        self.output_browser.append(lines)

            job_args = self.JobsQueue[row_index]
            if platform.system() == 'Windows':
                program = '"' + GlobalSettings.appdir + "Casper_Seq_Finder.exe" + '" '
            else:
                program = '"' + GlobalSettings.appdir + "Casper_Seq_Finder" + '" '
            program += job_args
            self.process.readyReadStandardOutput.connect(partial(output_stdout, self.process))
            self.process.start(program)
        else:
            error = QtWidgets.QMessageBox.critical(self, "No Jobs To Run",
                                                   "No jobs are in the queue to run. Please add a job before running.",
                                                   QtWidgets.QMessageBox.Ok)


    def upon_process_finishing(self):
        row_index = self.indexes[0]
        name = self.job_Table.item(row_index, 1).text()
        item = QtWidgets.QTableWidgetItem(name)
        item.setTextAlignment(QtCore.Qt.AlignHCenter)
        self.job_Table.setItem(row_index, 2, item)
        self.job_Table.setItem(row_index, 1, QtWidgets.QTableWidgetItem(""))
        self.indexes.pop(0)
        if len(self.indexes) != 0:
            self.run_job()

    def clear_job_queue(self):
        self.process.kill()
        self.fin_index = 0
        self.job_Table.clearContents()
        self.job_Table.setRowCount(0)
        self.JobsQueue = []
        self.orgName.clear()
        self.strainName.clear()
        self.orgCode.clear()
        self.output_browser.clear()
        self.s_file.clear()
        self.progressBar.setValue(0)
        self.s_file.setText("Name of File")
        self.first = False
        self.s_file.clear()

    def reset(self):
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_3.clear()
        self.first = False
        self.s_file.clear()


    def open_ncbi_web_page(self):
        webbrowser.open('https://www.ncbi.nlm.nih.gov/', new=2)


    def visit_repo_func(self):
        webbrowser.open('https://github.com/TrinhLab/CASPERapp')


    def closeEvent(self, event):
        # make sure that there are cspr files in the DB
        file_names = os.listdir(GlobalSettings.CSPR_DB)
        noCSPRFiles = True
        for file in file_names:
            if 'cspr' in file:
                noCSPRFiles = False
                break
        if noCSPRFiles == True:
            if self.exit == False:
                error = QtWidgets.QMessageBox.question(self, "No CSPR file generated",
                                                        "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                                                        "Alternatively, you could quit the program. Would you like to quit?",
                                                        QtWidgets.QMessageBox.Yes |
                                                        QtWidgets.QMessageBox.No,
                                                        QtWidgets.QMessageBox.No)
                if (error == QtWidgets.QMessageBox.No):
                    event.ignore()
                else:
                    event.accept()
            else:
                self.exit = False
                event.accept()

        else:
            self.process.kill()
            self.job_Table.clearContents()
            self.job_Table.setRowCount(0)
            self.orgName.clear()
            self.strainName.clear()
            self.orgCode.clear()
            self.output_browser.clear()
            self.s_file.setText("Name of File")
            self.progressBar.setValue(0)
            self.first = False
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            GlobalSettings.mainWindow.mwfg.moveCenter(GlobalSettings.mainWindow.cp)  ##Center window
            GlobalSettings.mainWindow.move(GlobalSettings.mainWindow.mwfg.topLeft())  ##Center window
            GlobalSettings.mainWindow.show()
            if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
            GlobalSettings.mainWindow.orgChoice.clear()
            GlobalSettings.mainWindow.endoChoice.clear()
            GlobalSettings.mainWindow.getData()
            GlobalSettings.MTWin.launch(GlobalSettings.CSPR_DB)
            GlobalSettings.pop_Analysis.launch(GlobalSettings.CSPR_DB)
            event.accept()


    def continue_to_main(self):
        # make sure that there are cspr files in the DB
        file_names = os.listdir(GlobalSettings.CSPR_DB)
        noCSPRFiles = True
        noCSPRFiles = True
        for file in file_names:
            if 'cspr' in file:
                noCSPRFiles = False
                break
        if noCSPRFiles == True:

            error = QtWidgets.QMessageBox.question(self, "No CSPR file generated",
                                                   "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                                                   "Alternatively, you could quit the program. Would you like to quit?",
                                                   QtWidgets.QMessageBox.Yes |
                                                   QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)


            if (error == QtWidgets.QMessageBox.Yes):
                self.exit = True
                self.close()

        else:
            self.process.kill()
            self.job_Table.clearContents()
            self.orgName.clear()
            self.strainName.clear()
            self.orgCode.clear()
            self.output_browser.clear()
            self.s_file.clear()
            self.progressBar.setValue(0)
            self.first = False
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            GlobalSettings.mainWindow.mwfg.moveCenter(GlobalSettings.mainWindow.cp)  ##Center window
            GlobalSettings.mainWindow.move(GlobalSettings.mainWindow.mwfg.topLeft())  ##Center window
            GlobalSettings.mainWindow.show()
            if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
            GlobalSettings.mainWindow.orgChoice.clear()
            GlobalSettings.mainWindow.endoChoice.clear()
            GlobalSettings.mainWindow.getData()
            self.hide()
