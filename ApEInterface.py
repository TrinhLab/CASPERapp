import os
import sys

class ApEFile:

    def __init__(self):
        self.filename = str()
        LOCUS = tuple()
        DEFINITION = '  .'
        ACCESSION = '   '
        VERSION = "     "
        SOURCE = "      ."
        ORGANISM
        COMMENT1
        COMMENT2
        FEATURES =
        ORIGIN = str()

    def Aread(self):


    def Awrite(self):
        f = open(self.filename, 'f')



class ape_feature:

    def __init__(self):
        name = str()
        pos = (0,0)  # start and stop position
        label = str()
        fwdcolor = str()
        revcolor = str()
        # graphicformat is a bit complex as it is the position of the arrow in the map
        graphicformat = "arrow_data {{0 1 2 0 0 -1} {} 0}\nwidth 5 offset 0"

