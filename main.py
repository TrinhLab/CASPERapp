import sys
import os
import io
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from APIs import Kegg, SeqFromFasta
from bioservices import KEGG
from Results import Results


# =========================================================================================
# CLASS NAME: CMainWindow
# Inputs: Takes in the path information from the startup window and also all input parameters
# that define the search for targets e.g. endonuclease, organism genome, gene target etc.
# Outputs: The results of the target search process by generating a new Results window
# =========================================================================================


class CMainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(CMainWindow, self).__init__()
        uic.loadUi('CASPER_main.ui', self)
        self.dbpath = ""
        self.data = {}
        self.shortHand ={}
        self.orgcodes = {}  # Stores the Kegg organism code by the format {full name : organism code}
        # --- Button Modifications --- #
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.pushButton_FindTargets.clicked.connect(self.gather_settings)
        self.pushButton_ViewTargets.clicked.connect(self.view_results)
        self.pushButton_ViewTargets.setEnabled(False)

        self.pushButton_search.clicked.connect(self.search_kegg)

        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()

        # --- Menubar commands --- #
        self.actionChange_Directory.triggered.connect(self.change_directory)

        # --- Setup for Gene Entry Field --- #
        self.geneEntryField.setPlainText("Example Inputs: \n"
                                               "Gene (LocusID): YOL086C  *for Saccharomyces Cerevisiae ADH1 gene* \n"
                                               "Position: (chromosome,start,stop)(chromosome,start,stop)...\n"
                                               "Sequence: *Pure sequence. CASPER will search for targets and report off"
                                               "targets based on the genome selected if any*")


        #show functionalities on window
        self.view_my_results = Results()

        #self.show()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def gather_settings(self):
        inputstring = str(self.geneEntryField.toPlainText())
        self.progressBar.setValue(10)
        if self.radioButton_Gene.isChecked():
            ginput = inputstring.split(',')
            self.run_results("gene", ginput)
        elif self.radioButton_Position.isChecked():
            pinput = inputstring.split(';')
            self.run_results("position", pinput)
        elif self.radioButton_Sequence.isChecked():
            sinput = inputstring
            self.run_results("sequence", sinput)


    # ---- IS ONLY CALLED FROM gather_settings!!!! ---- #
    def run_results(self, inputtype, inputstring):
        kegginfo = Kegg()
        org = str(self.orgChoice.currentText())
        endo = str(self.endoChoice.currentText())
        ginfo = {}  # each entry ginfo[gene] = (chromosome number, t/f strand, start pos, end pos)
        progvalue = 15
        self.progressBar.setValue(progvalue)
        if inputtype == "gene":
            for gene in inputstring:
                g = self.shortHand[org] + ":" + gene

                hold = kegginfo.gene_locator(g)
                if hold == -1:
                    QtWidgets.QMessageBox.question(self, "Gene Database Error", "The Gene you entered could not be found in the Kegg database. Please make sure you entered everything correctly and try again.",
                                                   QtWidgets.QMessageBox.Ok)
                    progvalue = 0
                    self.progressBar.setValue(progvalue)
                    return
                ginfo[gene]=hold
                progvalue += 50/len(inputstring)
                self.progressBar.setValue(progvalue)
        if inputtype == "position":
            ginfo = inputstring[1:-1].split(",")
            self.progressBar.setValue(45)
        if inputtype == "sequence":
            self.progressBar.setValue(45)
        s = SeqFromFasta()
        filename = self.dbpath + org + ".fna"
        s.setfilename(filename)
        progvalue = 75
        self.progressBar.setValue(progvalue)
        for gene in inputstring:
            s.getsequenceandtargets(ginfo[gene], 100, 100, self.dbpath+'/'+org, endo)
            progvalue += 25/len(inputstring)
            self.progressBar.setValue(progvalue)
            self.view_my_results.loadGenesandTargets(s.getgenesequence(), ginfo[gene][2]+100, ginfo[gene][3]-100,
                                                     s.gettargets(), gene)
        self.progressBar.setValue(100)
        self.pushButton_ViewTargets.setEnabled(True)

    # This function works as a way to look up a search term in the Kegg database to potentially get the code
    # for the gene
    def search_kegg(self):
        k = KEGG()
        current_org = self.orgcodes[self.orgChoice.currentText()]
        info = k.find(current_org, self.lineEdit_search.text())
        print(info)

    # This method is for testing the execution of a button call to make sure the button is linked properly
    def testexe(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to Quit?",
                                            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            print(self.orgChoice.currentText())
            sys.exit()
        else:
            pass

    # This method checks if a check button or radio button works appropriately by printing the current organism
    def testcheckandradio(self):
         print(str(self.orgChoice.currentText()))

    # ----- CALLED IN STARTUP WINDOW ------ #
    def getData(self):
        mypath = os.getcwd()
        self.dbpath = mypath
        onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
        print(onlyfiles)
        orgsandendos = {}
        shortName = {}
        for file in onlyfiles[1:]:
            if file.find('.cspr')!=-1:
                newname = file[0:-4]
                s = newname.split('_')
                hold = open(file)
                buf = (hold.readline())
                species = buf[8:buf.find('\n')]
                endo = str(s[1])
                if species not in shortName:
                    shortName[species] = s[0]
                if species in orgsandendos:
                    orgsandendos[species].append(endo)
                else:
                    orgsandendos[species] =[endo]
                    self.orgChoice.addItem(species)


        self.data = orgsandendos
        self.shortHand= shortName
        self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])
        self.orgChoice.currentIndexChanged.connect(self.changeEndos)

        #os.chdir('/Users/brianmendoza/PycharmProjects/CASPERapp/')
        """f = open('CASPERinfo')
        while True:
            line = f.readline()
            if line.startswith('ORGA'):
                while True:
                    orginfo = f.readline()
                    orginfo = orginfo[0:-1]
                    if orginfo[0] == '-':
                        break
                    stuff = orginfo.split(":")
                    self.orgcodes[stuff[1]] = stuff[0]
                break
        f.close()

        for item in self.orgcodes:
            self.orgChoice.addItem(item)
        self.endoChoice.addItems(self.data[self.orgcodes[str(self.orgChoice.currentText())]])
        self.orgChoice.currentIndexChanged.connect(self.changeEndos)
"""
    def changeEndos(self):

        self.endoChoice.clear()
        self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])

    def change_directory(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                           self.dbpath, QtWidgets.QFileDialog.ShowDirsOnly)
        os.chdir(mydir)
        self.getData()

    @QtCore.pyqtSlot()
    def view_results(self):
        self.view_my_results.show()
        self.progressBar.setValue(0)
        self.pushButton_ViewTargets.setEnabled(False)


# ----------------------------------------------------------------------------------------------------- #

# =========================================================================================
# CLASS NAME: StartupWindow
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class StartupWindow(QtWidgets.QDialog):
    def __init__(self):
        super(StartupWindow, self).__init__()
        uic.loadUi('startupCASPER.ui', self)
        self.setWindowTitle('WELCOME TO CASPER!')
        self.setWindowModality(2)  # sets the modality of the window to Application Modal

        #---Button Modifications---#
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        pixmap = QtGui.QPixmap('mainart.jpg')
        self.labelforart.setPixmap(pixmap)
        self.pushButton_2.setDefault(True)

        self.gdirectory = os.path.expanduser("~")
        self.gdirectory = "Please select a Directory that contains .capr files"  # Temporary. Throw away at deployment
        self.lineEdit.setText(self.gdirectory)

        self.pushButton_3.clicked.connect(self.changeDir)
        self.pushButton_2.clicked.connect(self.show_window)
        self.pushButton.clicked.connect(self.errormsgmulti)

        self.show_main_window = CMainWindow()

        #show functionalities on window
        self.show()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def changeDir(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                       self.gdirectory, QtWidgets.QFileDialog.ShowDirsOnly)
        self.lineEdit.setText(mydir)
        cdir = self.lineEdit.text()
        os.chdir(mydir)
        self.gdirectory = mydir
        print(mydir)
        print(cdir)

    def errormsgmulti(self):
        QtWidgets.QMessageBox.question(self, "Under Construction...", "Sorry this functionality is still"
                                            " under construction and will be available shortly!",
                                            QtWidgets.QMessageBox.Ok)

    @QtCore.pyqtSlot()
    def show_window(self):
        if(self.gdirectory=="Please select a Directory that contains .capr files"):
            QtWidgets.QMessageBox.question(self, "Must select directory", "You must select your directory",
                                                                                      QtWidgets.QMessageBox.Ok)
        else:
            os.chdir(self.gdirectory)
            self.show_main_window.show()
            self.show_main_window.getData()
            self.close()


if __name__ == '__main__':
    app = Qt.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")
    startup = StartupWindow()
    sys.exit(app.exec_())
