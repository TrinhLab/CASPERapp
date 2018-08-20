"""This file creates the necessary file links needed to run the off-target subprocess of C++. """

import subprocess
from Algorithms import SeqTranslate

class OffAnalysis:

    def __init__(self, off_settings_file):
        # Containers for all necessary file links
        self.csprfiles = list()
        self.concise_file = str()
        self.off_triggered = dict()  # keys are the sequences decompressed, values are the scores and other sequences

        self.retrieve_settings(off_settings_file)

    def retrieve_settings(self, path):
        for line in open(path,"r"):
            start = line.find(":")

    def append_to_concise(self):
        f = open(self.concise_file)
        fnew = open(self.concise_file[:-4] + "_withoff.csv", "w")
        for line in f:
            if line[selectindex] in self.off_triggered:
                fnew.write(line.append(self.off_triggered[line[selectindex]]))
            else:
                fnew.write(line)
        f.close()
        fnew.close()




    def run_off_target(self):
        for file in self.csprfiles:


