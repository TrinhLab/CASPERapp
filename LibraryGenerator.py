"""This file holds the main class and any supplementary classes for forming a sgRNA library of targets."""

import os
import Algorithms
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic


class Library(QtWidgets.QDialog):

    def __init__(self):
        super(Library, self).__init__()
        uic.loadUi('library_prompt.ui', self)
        self.library_name = ""

        self.Targets = list()  # List of targets that contain information to be printed to the .csv file

    def generate_library(self):
        # get all the information from the user input boxes
        input_file = self.filename_input.text()
        target_range = (int(self.startline.text())/100, int(self.endline.text())/100)
        num_guides = int(self.numGenescomboBox.currentText())
        min_on_score_filter = int(self.minON_comboBox.currentText())
        seq_spec = self.fiveprimeseq.text()




    def print_library(self, output_file_path):
        f = open(output_file_path,'w')
        for item in self.Targets:
            f.write(item[0])
            for subsect in item:
                f.write("," + item)
