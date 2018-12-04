
import operator
import sys
import os

from Algorithms import SeqTranslate
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic



class CoTargeting(QtWidgets.QMainWindow):

    def __init__(self, info_path):

        super(CoTargeting, self).__init__()
        uic.loadUi('CoTargeting.ui', self)
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.path= ""
        self.info_path  = info_path
        self.organisms =""
        self.shortHand = {}
        self.OrganismDrop.currentIndexChanged.connect(self.endo_Changed)
        self.AddEndoButton.clicked.connect(self.add_Endo)
        self.RemoveEndoButton.clicked.connect(self.remove_Endo)
        self.CompairButton.clicked.connect(self.compare_Endos)
        self.Endos_All = []
        self.Endos_Added = []
        self.progress = 0;
        self.comp_endo_1 = ""
        self.comp_endo_2 = ""
        self.sq = SeqTranslate()
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()
        self.AnalyzeButton.setEnabled(False)

#########UI Function Code####################

    def launch(self, organisms,path,shorthand):
        self.organisms = organisms
        self.path = path
        self.shortHand = shorthand
        self.fill_Org()
        self.AddedEndos.setColumnCount(1)
        self.show()

    def fill_Org(self):
        self.AllEndos.clear()
        self.AllEndos.setColumnCount(1)

        for item in self.organisms:

            self.OrganismDrop.addItem(item)
        """for endo in self.organisms[self.OrganismDrop.currentText()]:
            self.Endos_All.append(endo)"""



    def fill_Tables(self):
        self.AddedEndos.clear()
        self.AllEndos.clear()
        if len(self.Endos_All)>0:
            self.AllEndos.setRowCount(len(self.Endos_All))
        else:
            self.AllEndos.setRowCount(1)
        if len(self.Endos_Added)>0:
            self.AddedEndos.setRowCount(len(self.Endos_Added))
        else:
            self.AddedEndos.setRowCount(1)
        index = 0
        for endo in self.Endos_All:
            hold = QtWidgets.QTableWidgetItem(endo)
            self.AllEndos.setItem(index, 0, hold)
            index += 1
        index = 0
        for endo in self.Endos_Added:
            hold = QtWidgets.QTableWidgetItem(endo)
            self.AddedEndos.setItem(index, 0, hold)
            index += 1
        if len(self.Endos_All)==0:
            self.AddEndoButton.setEnabled(False)
        else:
            self.AddEndoButton.setEnabled(True)
        if len(self.Endos_Added)==0:
            self.RemoveEndoButton.setEnabled(False)
        else:
            self.RemoveEndoButton.setEnabled(True)

    def endo_Changed(self):
        self.Endos_All.clear()
        self.Endos_Added.clear()

        for endo in self.organisms[self.OrganismDrop.currentText()]:
            self.Endos_All.append(endo)
        self.fill_Tables()

    def add_Endo(self):
        index = self.AllEndos.currentRow()
        self.Endos_Added.append(self.Endos_All[index])
        self.Endos_All.remove(self.Endos_All[index])
        self.fill_Tables()

    def remove_Endo(self):
        index = self.AddedEndos.currentRow()
        self.Endos_All.append(self.Endos_Added[index])
        self.Endos_Added.remove(self.Endos_Added[index])
        self.fill_Tables()

##############Data Processing code####################

    def import_data(self,endo):
        file_name = self.shortHand[self.OrganismDrop.currentText()]+"_"+endo
        if self.path.find("/") != -1:
            file = open(self.path+"/"+file_name+".cspr")
        else:
            file = open(self.path + "\\" + file_name + ".cspr")
        info = file.read()
        sorted_info = {}
        info = info.split("\n")
        chromo = ""
        for line in info:

            if "GENOME" in line:
                continue
            if "REPEATS" in line:
                break
            if "CHROMOSOME" in line:
                chromo = line[12:]
                sorted_info[chromo] = {}
                continue
            sorted_info[str(chromo)][line.split(',')[0]] = line
        self.progress += round(10/len(self.Endos_Added))
        self.progressBar.setValue(self.progress)
        return sorted_info

    def compare_Endos(self):
        if len(self.Endos_Added)<2:
            QtWidgets.QMessageBox.question(self, "Insufficient endonucleases selected",
                                           "Please select at least two endonucleases to compare.",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.comp_endo_1 = self.import_data(self.Endos_Added[0])

        index= 1
        while index<len(self.Endos_Added):
            self.comp_endo_2= self.import_data(self.Endos_Added[index])
            self.comp_endo_1 = self.compare_two()
            index+=1
        self.completed_Compair()

    def compare_two(self):
        shared = {}
        for chromo in self.comp_endo_1:
            if not(chromo in self.comp_endo_2):
                continue
            for line_1 in self.comp_endo_1[chromo]:
                if line_1 in self.comp_endo_2[chromo]:
                    if chromo in shared:
                        shared[chromo][line_1] = self.comp_endo_1[chromo][line_1]
                    else:
                        shared[chromo] = {}
                        shared[chromo][line_1] = self.comp_endo_1[chromo][line_1]
                self.progress+= 60/((len(self.Endos_Added)-1)*len(self.comp_endo_1)*len(self.comp_endo_1[chromo]))
                self.progressBar.setValue(self.progress)
        return shared

    def completed_Compair(self):

        if len(self.comp_endo_1)==0:
            QtWidgets.QMessageBox.question(self, "No Matches",
                                           "No shared sequences were found in"
                                           " \nthe endonucleases you selected",
                                           QtWidgets.QMessageBox.Ok)
            self.progress = 0
            self.progressBar.setValue(self.progress)
            return
        file_name = self.shortHand[self.OrganismDrop.currentText()]+"_"
        begin = True
        for endo in self.Endos_Added:
            if begin:
                file_name = file_name+endo
                begin=False
            else:
                file_name = file_name+"&"+endo
        file = open(self.path + "\\" + file_name + ".cspr","w+")
        file.write("GENOME: "+self.OrganismDrop.currentText()+"\n")
        for chromo in self.comp_endo_1:
            file.write("CHROMOSOME #"+chromo+"\n")
            for found in self.comp_endo_1[chromo]:
                file.write(self.comp_endo_1[chromo][found]+"\n")
                self.progress = self.progress+30/(len(self.comp_endo_1)*len(self.comp_endo_1[chromo]))
                self.progressBar.setValue(self.progress)
        file.write("REPEATS\n")
        file.write("END_OF_FILE")
        file.close()
        self.progress = 100
        self.progressBar.setValue(self.progress)
        x=3
        QtWidgets.QMessageBox.question(self, "File Created",
                                       "A casper file named "+file_name+" has been"
                                        " \ncreated in your directory with all found matches",
                                       QtWidgets.QMessageBox.Ok)
        self.progress = 0
        self.progressBar.setValue(self.progress)





