from utils.sequence_utils import SeqTranslate
import gzip

##################################################################################################################################
# Use: Use as a parser for the cspr files
# Precondition: Only to the used with .cspr files. Will not work with any other files
# This class also took some of the parsing functions from with classes (Multitargeting and Results) and stores them in here
##################################################################################################################################

class CSPRparser:
    def __init__(self, inputFileName, casper_info_path):
        self.fileName = inputFileName
        self.seqTrans = SeqTranslate(casper_info_path)

    def read_targets(self, genename, pos_tuple, endo):
        retList = []
        with open(self.fileName, 'r') as f:
            for line in f:
                if line.startswith(f">{pos_tuple[0]}"):
                    for line in f:
                        if line.startswith('>'):
                            break
                        line = line.strip().split(',')
                        if pos_tuple[1] <= abs(int(line[0])) < pos_tuple[2]:
                            strand = "-" if int(line[0]) < 0 else "+"
                            retList.append((line[0], line[1], line[2], line[3], strand, endo))
        return retList
