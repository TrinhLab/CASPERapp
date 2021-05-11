import sys, os
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
import GlobalSettings
from PyQt5.QtGui import QIntValidator
class NewEndonuclease(QtWidgets.QDialog):

    def __init__(self):
        super(NewEndonuclease, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'newendonuclease.ui', self)
        self.setWindowTitle('New Endonuclease')
        self.error = False
        pamFlag = False;

        self.onList = []
        self.offList = []

        self.onList, self.offList = self.get_on_off_data() ### Call function to fill on- and off- data name lists

        for name in self.onList: ### Add on-target names to drop-down
            self.comboBox.addItem(str(name))

        for name in self.offList: ### Add off-target names to drop-down
            self.comboBox_2.addItem(str(name))


        self.submit_button.clicked.connect(self.submit)
        self.cancel_button.clicked.connect(self.cancel)

        self.seed_length.setValidator(QIntValidator(0,100,self.seed_length))
        self.five_length.setValidator(QIntValidator(0,100,self.five_length))
        self.three_length.setValidator(QIntValidator(0,100,self.three_length))

        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        font: 15pt "Arial";
                        font: bold;
                        margin-top: 10px;}"""

        self.groupBox.setStyleSheet(groupbox_style)
        self.groupBox_2.setStyleSheet(groupbox_style.replace("groupBox","groupBox_2"))
        self.groupBox_3.setStyleSheet(groupbox_style.replace("groupBox","groupBox_3"))

    def writeNewEndonuclease(self, newEndonucleaseStr):
        with open(GlobalSettings.appdir + 'CASPERinfo', 'r') as f, open(GlobalSettings.appdir + "new_file", 'w+') as f1:
            for line in f:
                f1.write(line)
                if 'ENDONUCLEASES' in line:
                    f1.write(newEndonucleaseStr + '\n')  # Move f1.write(line) above, to write above instead
        os.remove(GlobalSettings.appdir + "CASPERinfo")
        os.rename(GlobalSettings.appdir + "new_file",
                  GlobalSettings.appdir + "CASPERinfo")  # Rename the new file

    def submit(self):
        # This is executed when the button is pressed
        name = str(self.organism_name.text())
        abbr = str(self.abbreviation.text())
        crisprtype = str(self.crispr_type.text())
        seed_len = str(self.seed_length.text())
        five_len = str(self.five_length.text())
        three_len = str(self.three_length.text())
        pam = str(self.pam_sequence.text())
        ### Check for multiple PAMs and format if present
        if len(pam.split(','))>0:
            pam = [x.strip() for x in pam.split(',')]
            pam = ",".join(pam)
        ### Check for PAM directionality
        if self.five_pam.isChecked():
            pam_dir = str(5)
        else:
            pam_dir = str(3)
        on_scoring = str(self.comboBox.currentText())
        off_scoring = str(self.comboBox_2.currentText())
        length = len(seed_len) + len(five_len) + len(three_len)
        argument_list = [abbr, pam, five_len, seed_len, three_len, pam_dir, name, crisprtype, on_scoring, off_scoring]
        validPAM = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',',')
        self.error = False;

        ### Error checking for PAM alphabet
        for letter in pam:
            if (letter.upper() not in validPAM):
                QtWidgets.QMessageBox.question(self, "Invalid PAM", "Invalid characters in PAM Sequence.",
                                               QtWidgets.QMessageBox.Ok)
                self.exec()
                return True
        ### Error checking for filling out all fields
        for arg in argument_list:
            if ';' in arg:
                QtWidgets.QMessageBox.question(self, "Invalid Semicolon", "Invalid character used: ';'",
                                               QtWidgets.QMessageBox.Ok)
                self.exec()
                return True
            elif arg == "":
                QtWidgets.QMessageBox.question(self, "Empty Field", "Please fill in all fields.", QtWidgets.QMessageBox.Ok)
                self.exec()
                return True
            else:
                pass
            
        myString = ""
        for i, arg in enumerate(argument_list):
            if i == len(argument_list)-1: ### Last argument in list
                myString += str(arg)
            else:
                myString += str(arg) + ";"
        
        self.writeNewEndonuclease(myString)
        self.clear_all()
        self.close()

    def cancel(self):
        self.clear_all()
        self.close() 

    ### This function clears all of the line edits
    def clear_all(self):
        self.organism_name.clear()
        self.abbreviation.clear()
        self.crispr_type.clear()
        self.seed_length.clear()
        self.five_length.clear()
        self.three_length.clear()
        self.pam_sequence.clear()

    ### This function parses CASPERinfo to return the names (in lists) of all on-target and off-target scoring data
    def get_on_off_data(self):
        filename = GlobalSettings.appdir + "CASPERinfo"
        retList_on = []
        retList_off = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = str(line)
                if "ON-TARGET DATA" in line:
                    index = i
                    while "-----" not in line:
                        if "DATA:" in line:
                            retList_on.append(line.split("DATA:")[-1].strip()) ### Append name of scoring data to on-target name list
                            line = lines[index+1]
                            index += 1
                        else:
                            line = lines[index+1]
                            index += 1
                            continue
                elif "OFF-TARGET MATRICES" in line:
                    index = i
                    while "-----" not in line:
                        if "MATRIX:" in line:
                            retList_off.append(line.split("MATRIX:")[-1].strip()) ### Append name of scoring data to off-target name list
                            line = lines[index+1]
                            index += 1
                        else:
                            line = lines[index+1]
                            index += 1
                            continue
                else:
                    continue
        return retList_on, retList_off