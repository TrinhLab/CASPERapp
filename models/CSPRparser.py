from utils.Algorithms import SeqTranslate
import gzip

##################################################################################################################################
# Use: Use as a parser for the cspr files
# Precondition: Only to the used with .cspr files. Will not work with any other files
# This class also took some of the parsing functions from with classes (Multitargeting and Results) and stores them in here
##################################################################################################################################

class CSPRparser:
    def __init__(self, inputFileName):
        self.multiSum = 0  # multitargetting sum taken from the previous version of make_graphs
        self.multiCount = 0  # multitargetting count taken from the previous version of make_graphs
        self.seqTrans = SeqTranslate()  # SeqTranslate variable. for decrompressing the data
        self.chromesomeList = list()  # list of a list for the chromesomes. As it currently stands, this variable is used in both read_chromesomes and in read_targets
        self.karystatsList = list()  # list of (ints) of the karyStats (whatever those are) to be used for the get_chrom_length function
        self.genome = ""  # genome name
        self.misc = ""  # anything from the misc line
        self.repeats = {}  # dictionary of the number of repeats. See the read_repeats function for more info
        self.seeds = {}  # dictionary of which chromesomes are repeats. See the read_repeats function for more info
        self.dec_tup_data = {}
        self.chromesomesSelectedList = list()
        # data for population analysis
        # dict:
        # key = the seed
        #       value = tuple (org name, chom #, location, sequence, pam, score, strand, endo)
        self.popData = {}

        # file path variable
        self.fileName = inputFileName

    def gen_lib_parser(self, genDict, endo):
        retDict = dict()
        for gene in genDict:
            retDict[gene] = list()
            retDict[gene] = self.read_targets('', (genDict[gene][0], genDict[gene][1], genDict[gene][2]), endo)
        return retDict

    def read_first_lines(self):
        with open(self.fileName, 'r') as f:
            i = 0
            for line in f:
                if i > 2:
                    break
                else:
                    line = str(line)
                    if i == 0:
                        colonIndex = line.find(':') + 2
                        buffer1 = line[colonIndex:]
                        self.genome = buffer1
                    elif i == 1:
                        colonIndex = line.find(':') + 2
                        k_data = line[colonIndex:]
                        k_data = k_data.split(',')
                        k_data = k_data[:-1]
                        for k in k_data:
                            self.karystatsList.append(int(k))
                    else:
                        colonIndex = line.find(':') + 2
                        buffer1 = line[colonIndex:]
                        self.misc = buffer1
                    i += 1

    def get_chromesome_names(self):
        self.chromesomesSelectedList.clear()
        with open(self.fileName,'r') as f:
            i = 0
            for line in f:
                line = str(line)
                if line == 'REPEATS':
                    break
                elif '>' in line:
                    self.chromesomesSelectedList.append(line)
                if i == 0:
                    retGen = line
                elif i == 2:
                    retMisc = line
                i += 1
        return retGen, retMisc

    def read_targets(self, genename, pos_tuple, endo):
        i = 0
        retList = []
        header = False
        with open(self.fileName, 'r') as f:
            for line in f:
                line = str(line).strip()
                if i > 2:
                    if '>' in line and '(' + str(pos_tuple[0]) + ')' in line:
                        header = True
                    elif header == True:

                        if line.find('>') != -1:
                            break
                        line = line.split(',')
                        if abs(int(line[0])) >= int(pos_tuple[1]) and abs(int(line[0])) < int(pos_tuple[2]):
                            if int(line[0]) < 0:
                                strand = "-"
                            else:
                                strand = "+"
                            retList.append((line[0], line[1], line[2], line[3], strand, endo))
                        elif int(line[0]) >= int(pos_tuple[2]):
                            break
                    elif line == 'REPEATS':
                        break
                i += 1
        return retList

    def uniq_seq_count(self):
        self.unique_targets = 0
        for chromo in self.chromesomeList:
            for data in chromo:
                if len(data) == 6:
                    self.unique_targets += 1
        return self.unique_targets