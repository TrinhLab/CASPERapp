"""This file runs the off target analysis for CASPER. Use the CASPEROfflist.txt file to set up the information you want
    to run with this program.
    WARNING: Running this protocol on a large number of sequences is unwise and may take significant computing power/time."""

"""Called:
O = OffTargetAlgorithm()"""


import os, sys, math, datetime
import subprocess
from PyQt5 import Qt, QtWidgets, uic, QtCore, QtGui
import subprocess as sub
from Algorithms import SeqTranslate

#off_list, cspr_file, max_misses=4, threshold=0.05, compressed=False, output_path=None, detailed=False, average=True
class OffTarget(QtWidgets.QDialog):

    def __init__(self):

        super(OffTarget, self).__init__()
        uic.loadUi('OffTargetAnalysis.ui', self)
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()

        #cmd = "dir"
        #p = os.system(cmd)
        #print(p)












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
