"""This file runs the off target analysis for CASPER. Use the CASPEROfflist.txt file to set up the information you want
    to run with this program.
    WARNING: Running this protocol on a large number of sequences is unwise and may take significant computing power/time."""

"""Called:
O = OffTargetAlgorithm()"""


import os, sys, math, datetime
import subprocess
from PyQt5 import Qt, QtWidgets, uic, QtCore, QtGui
from functools import partial
import subprocess
import GlobalSettings
import threading
import queue
import shlex
import time
from Algorithms import SeqTranslate

#off_list, cspr_file, max_misses=4, threshold=0.05, compressed=False, output_path=None, detailed=False, average=True
class OffTarget(QtWidgets.QDialog):

    def __init__(self):

        super(OffTarget, self).__init__()
        uic.loadUi(GlobalSettings.appdir+'\OffTargetAnalysis.ui', self)
        self.show()
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()
        self.Run.clicked.connect(self.run_analysis)
        self.tolerancehorizontalSlider.valueChanged.connect(self.tol_change)
        self.tolerancehorizontalSlider.setMaximum(100)
        self.tolerancehorizontalSlider.setMinimum(0)
        self.tolerance = 0.0
        self.tolerancelineEdit.setText("0")
        #self.tolerancelineEdit.textChanged.connect(self.tol_change)
        self.pushButton.clicked.connect(self.tol_change)
        self.buttonBox.clicked.connect(self.exit)
        self.fill_data_dropdown()
        self.perc = False
        self.bool_temp = False
        self.running = False
        self.process = QtCore.QProcess()

    def fill_data_dropdown(self):
        #fill in chromosome and endo dropdowns
        onlyfiles = [f for f in os.listdir(GlobalSettings.filedir) if os.path.isfile(os.path.join(GlobalSettings.filedir , f))]
        self.orgsandendos = {}
        self.shortName = {}
        for file in onlyfiles:
            if file.find('.cspr') != -1:
                newname = file[0:-4]
                s = newname.split('_')
                hold = open(os.path.join(GlobalSettings.filedir , file))
                buf = (hold.readline())
                species = buf[8:buf.find('\n')]
                endo = str(s[1])
                if species not in self.shortName:
                    self.shortName[species] = s[0]
                if species in self.orgsandendos:
                    self.orgsandendos[species].append(endo)
                else:
                    self.orgsandendos[species] = [endo]
                    self.OrgcomboBox.addItem(species)
        self.data = self.orgsandendos
        self.shortHand = self.shortName
        temp = self.data[str(self.OrgcomboBox.currentText())]
        temp1 = []
        for i in temp:
            i = i.strip('.')
            temp1.append(i)
        self.EndocomboBox.addItems(temp1)
        self.OrgcomboBox.currentIndexChanged.connect(self.changeEndos)

        #fill in Max Mismatch dropdown
        mismatch_list = ['1','2','3','4','5','6','7','8','9','10']
        self.mismatchcomboBox.addItems(mismatch_list)

    def changeEndos(self):
        self.EndocomboBox.clear()
        temp = self.data[str(self.OrgcomboBox.currentText())]
        temp1 = []
        for i in temp:
            i = i.strip('.')
            temp1.append(i)
        self.EndocomboBox.addItems(temp1)

    def tol_change(self):
        if(self.tolerance == float(self.tolerancelineEdit.text())):
            self.tolerance = self.tolerancehorizontalSlider.value() / 100 * 0.5
            self.tolerance = round(self.tolerance, 3)
            self.tolerancelineEdit.setText(str(self.tolerance))
        else:
            self.tolerance = float(self.tolerancelineEdit.text())
            self.tolerance = round(self.tolerance, 3)
            self.tolerancehorizontalSlider.setValue(round(self.tolerance/0.5 * 100))

    def run_command(self):

        file_name_1 = self.shortName[str(self.OrgcomboBox.currentText())]
        file_name_2 = self.orgsandendos[str(self.OrgcomboBox.currentText())]
        file_name = str(file_name_1) + '_' + str(file_name_2[0]) + 'cspr'
        print(file_name)

        if (self.AVG.isChecked()):
            avg_output = r' True '
            detailed_output = r' False '
        else:
            avg_output = r' False '
            detailed_output = r' True '


        app_path = GlobalSettings.appdir
        exe_path = app_path + r'\OffTargetFolder\CasperOffTargetWindows '
        exe_path = '"' +  exe_path + '"'
        data_path = ' "' + app_path + '\\OffTargetFolder\\temp.txt' + '" ' ##
        compressed = r' True ' ##
        cspr_path = ' "' + os.getcwd() + '\\' + file_name + '" '
        output_path = ' "' + os.getcwd() + '\\' + self.FileName.text() + '_OffTargetResults.txt" '
        CASPER_info_path = r' "' + app_path + '\\CASPERinfo' + '" '
        num_of_mismathes = int(self.mismatchcomboBox.currentText())
        tolerance = self.tolerance

        cmd = exe_path + data_path + compressed + cspr_path + output_path + CASPER_info_path + str(num_of_mismathes) + ' ' + str(tolerance) + detailed_output + avg_output
        print(cmd)

        def finished():
            self.running = False
            self.progressBar.setValue(100)

        def dataReady(p):
            line = str(p.readAllStandardOutput())
            line = line[2:]
            line = line[:len(line)-1]
            for lines in filter(None,line.split(r'\r\n')):
                print(lines)
                if(lines.find("Running Off Target Algorithm for") != -1 and self.perc == False):
                    self.perc = True
                if(self.perc == True and self.bool_temp == False and lines.find("Running Off Target Algorithm for") == -1):
                    lines = lines[32:]
                    lines = lines.replace("%","")
                    if(float(lines) <= 99):
                        num = float(lines)
                        self.progressBar.setValue(num)
                    else:
                        self.bool_temp = True


        self.process.readyReadStandardOutput.connect(partial(dataReady,self.process))
        self.progressBar.setValue(0)
        QtCore.QTimer.singleShot(100, partial(self.process.start, cmd))
        self.process.finished.connect(finished)

    def run_analysis(self):
        if(self.running == False):
            self.running = True
            self.run_command()



        # self.output_file_path = "temp_off_results_file.txt"
        # self.detailed = False
        # self.average = False
        #
        # # Fill the list of variables that get passed to the subprocess:
        # self.args_array = list()  # This is the list that gets passed into the file
        # self.args_array.append("C:\\Users\\Greg...")
        # self.args_array.append(off_list)
        # if compressed:
        #     self.args_array.append("True")
        # else:
        #     self.args_array.append("False")
        # self.args_array.append(cspr_file)
        # if output_path:
        #     self.args_array.append(output_path)
        #     self.output_file_path = output_path
        # self.args_array.append("C:\\Users\\Greg...")
        # self.args_array.append(str(max_misses))
        # self.args_array.append(str(threshold))
        # if detailed:
        #     self.args_array.append("True")
        #     self.detailed = True
        # else:
        #     self.args_array.append("False")
        # if average:
        #     self.args_array.append("True")
        #     self.average = True
        # else:
        #     self.args_array.append("False")

    def exit(self):
        self.perc = False
        self.bool_temp = False
        self.running = False
        self.process.kill()
        self.hide()

    def closeEvent(self, event):
        self.process.kill()
        self.perc = False
        self.bool_temp = False
        self.running = False
        event.accept()





    def run_off_target(self):
        sub.run(self.args_array)

    def return_results(self):
        results_dict = dict()
        off_list = list()
        f = open(self.output_file_path)
        if self.average:
            sequence = f.readline().split(":")[0]
            score = f.readline()[:-1].split(":")[1]
        else:
            sequence = f.readline()[:-2]
            score = "--.--"
        if self.detailed:
            off_target = f.readline()
            while off_target != "\n":
                off_list.append(off_target[:-1].split(","))
        f.close()
        return results_dict
