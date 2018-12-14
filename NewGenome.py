import sys, os
import subprocess as sub
from threading import Thread
from queue import Queue, Empty
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
from bioservices import KEGG
import GlobalSettings


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
        uic.loadUi('NewGenome.ui', self)
        self.setWindowTitle('New Genome')
        self.k = KEGG()
        self.info_path = info_path
        #---Button Modifications---#

        self.setWindowIcon(Qt.QIcon("cas9image.png"))
        self.whatsthisButton.clicked.connect(self.whatsthisclicked)
        self.KeggSearchButton.clicked.connect(self.updatekegglist)
        self.resetButton.clicked.connect(self.reset)
        self.submitButton.clicked.connect(self.submit)
        self.selectFastaButton.clicked.connect(self.selectFasta)
        self.JobsQueBox.setReadOnly(True)

        self.viewStatButton.setEnabled(False)

        self.JobIndexInProgress = 0
        self.JobsQueue = []  # holds Job classes.
        self.Endos = []
        self.file = ""

        self.process = QtCore.QProcess(self)




        #show functionalities on window
        self.fillEndo()
        #self.show()






    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def selectFasta(self):

        filed = QtWidgets.QFileDialog()
        documents= os.path.expanduser("~")+"/Documents"
        self.file = QtWidgets.QFileDialog.getOpenFileName(filed,"Open File",documents)
        fileName = ""
        index = 0
        self.file = self.file[0]
        while True:
            if self.file[len(self.file)-1-index]=="/" or self.file[len(self.file)-1-index]=="\\" :
                break
            fileName = self.file[len(self.file)-1-index]+fileName
            index+=1





        self.nameFile.setText(fileName)
        """cdir = self.lineEdit.text()
        os.chdir(mydir)
        self.gdirectory = mydir
        print(mydir)
        print(cdir)"""

    def submit(self):
        warning = ""
        if len(self.lineEdit_1.text())==0:
            warning = warning + "\nYou need to include the organism's name."
        if len(self.file) == 0:
            warning = warning + "\nYou need to select a file."
        if len(warning) !=0:
            QtWidgets.QMessageBox.information(self, "Required Information", warning, QtWidgets.QMessageBox.Ok)
            return
        else:
            self.JobsQueue.append(CasperJob(self.lineEdit_1.text() + self.lineEdit_2.text(),
                                            self.comboBoxEndo.currentText(), self.lineEdit_3.text(), self.nameFile))

        if len(self.lineEdit_2.text()) == 0:
            warning = warning + "\nIt is recommended to include the organism's subspecies/strain."
        if len(self.lineEdit_3.text()) == 0:
            warning = warning + "\nYou must include an organism code (KEGG code recommended)."
            return
        if len(warning)!=0:
            hold = QtWidgets.QMessageBox.question(self, "Missing information", warning+
                                        "\n\nDo you wish to continue without including this information?"
                                       , QtWidgets.QMessageBox.Yes |
                                       QtWidgets.QMessageBox.No,
                                       QtWidgets.QMessageBox.No)
            if hold == QtWidgets.QMessageBox.No:
                return

        Job = job(self.lineEdit_1.text(), self.lineEdit_2.text(), self.lineEdit_3.text(),
                  self.comboBoxEndo.currentText(), self.file)
        self.JobsQue.append(Job)
        nxtLine=""
        if len(self.JobsQueBox.toPlainText())!=0:
            nxtLine = "\n"
        self.JobsQueBox.setPlainText(self.JobsQueBox.toPlainText()+nxtLine+Job.Organism_Name)

    def fillEndo(self):
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            f = open(self.info_path + "\\CASPERinfo")
        else:
            f = open(self.info_path + "/CASPERinfo")
        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if(line[0]=="-"):
                        break
                    endo = line[:line.find("\t")]
                    while True:
                        line = line[line.find("\t")+1:]
                        index = line.find("\t")
                        if(index ==-1):
                            self.Endos.append(endo+" (Pam: " +line[:len(line)-1]+")")
                            break
                        self.Endos.append(endo+" "+"(Pam: " +line[:index]+")")

                break
        f.close()
        self.comboBoxEndo.addItems(self.Endos)

    def findFasta(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to Quit?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            sys.exit()
        else:
            pass

    def testcheckandradio(self, state):
        if state == QtCore.Qt.Checked:
            pass

    def whatsthisclicked(self):
        QtWidgets.QMessageBox.information(self, "Organism Code", "The organism code is the manner in which CASPER will"
                                                                 "label it's data files and references for the organism"
                                                                 "you are importing here. It is HIGHLY RECOMMENDED that"
                                                                 "you use the 3-4 letter code used by KEGG as this will"
                                                                 "aid in automatic accession of annotations from the"
                                                                 "database.", QtWidgets.QMessageBox.Ok)

    def updatekegglist(self):

        self.keggSuggested.clear()
        kegg_orglist = self.k.lookfor_organism(self.lineEdit_1.text())
        holder = 0
        self.keggSuggested.setColumnCount(2)

        for item in kegg_orglist:

            second_space = item[item.find(" ") + 1:].find(" ")+item.find(" ")
            code = item[item.find(" ")+1:second_space+1]
            item = item[second_space+2:]
            semi = item.find(";")
            index = 1
            while True:
                if item[semi - index] == " ":
                    break
                index = index + 1
            organism = item[:semi - index]
            self.keggSuggested.setRowCount(holder+1)
            table_code = QtWidgets.QTableWidgetItem(code)
            table_organism = QtWidgets.QTableWidgetItem(organism)
            self.keggSuggested.setItem(holder, 0, table_organism)
            self.keggSuggested.setItem(holder, 1, table_code)
            # self.keggsearchresults.insertPlainText(item)
            holder+=1
        self.keggSuggested.resizeColumnsToContents()

    def run_jobs(self):
        # Top layer for loop to go through all of the jobs in the queue:
        for job in self.JobsQueue:
            self.process.setProgram(GlobalSettings.CASPER_FOLDER_LOCATION + "/Casper_Seq_Finder_" + GlobalSettings.OPERATING_SYSTEM_ID)
            self.process.setArguments(job.get_arguments())
            if self.process.state() == QtCore.QProcess.NotRunning:
                self.JobInProgress.setText(job.name)
                self.process.start(GlobalSettings.CASPER_FOLDER_LOCATION + "/Casper_Seq_Finder_" + GlobalSettings.OPERATING_SYSTEM_ID, job.get_arguments())
            newdatafile = self.lineEdit_3.text() + self.comboBoxEndo.currentText() + ".txt"
            newdatafile = str(os.curdir) + newdatafile
            f = open(newdatafile, 'w')
            # Write the organism name and subspecies/strain if applicable to the top of the file
            f.write(self.lineEdit_1.text() + self.lineEdit_2.text() + self.comboBoxEndo.currentText())

    def reset(self):
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_3.clear()
        self.keggSuggested.clear()
        self.nameFile.setText("Name Of File")


class CasperJob:
    def __init__(self, org, endo, org_code, ref_file):
        self.name = endo + " targets in " + org
        self.organism_name = org
        self.organism_code = org_code
        self.endo_name = endo[:endo.find("PAM:")-1]
        self.endo_pam = endo[endo.find("PAM: "):]
        self.reference_file = ref_file

        # These are endonuclease specific settings that should be pulled from CASPERinfo
        self.anti = bool()
        self.sequence_length = int()
        self.seed_length = int()

    def get_arguments(self):
        ret_array = [self.endo_name, self.endo_pam, self.organism_code]
        # attach the 5' or 3' direction
        if self.anti:
            ret_array.append("TRUE")
        else:
            ret_array.append("FALSE")
        ret_array.append(GlobalSettings.CSPR_DB)
        ret_array.append(GlobalSettings.CASPER_FOLDER_LOCATION + "/CRISPRscan.txt")
        ret_array.append(self.reference_file)
        ret_array.append(self.organism_name)
        ret_array.append(self.sequence_length)
        ret_array.append(self.seed_length)


