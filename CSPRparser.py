from Algorithms import SeqTranslate
import gzip

##################################################################################################################################
# CLASS NAME: CSPRparser
# Use: Use as a parser for the cspr files
# Precondition: Only to the used with .cspr files. Will not work with any other files
# This class also took some of the parsing functions from with classes (Multitargeting and Results) and stores them in here
##################################################################################################################################

class CSPRparser:

    def __init__(self, inputFileName):

        # variables used in this class
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
        print('gen lib parser')
        retDict = dict()
        for gene in genDict:
            retDict[gene] = list()
            retDict[gene] = self.read_targets('', (genDict[gene][0], genDict[gene][1], genDict[gene][2]), endo)
        return retDict


    def read_first_lines(self):
        print('read first lines')
        with gzip.open(self.fileName, 'r') as f:
            i = 0
            for line in f:
                if i > 2:
                    break
                else:
                    line = str(line)
                    line = line.strip("'b")
                    line = line[:len(line) - 5]
                    if i == 0:
                        colonIndex = line.find(':') + 2
                        buffer1 = line[colonIndex:]
                        self.genome = buffer1
                    elif i == 1:
                        colonIndex = line.find(':') + 2
                        k_data = line[colonIndex:]
                        k_data = k_data.split(',')
                        for k in k_data:
                            self.karystatsList.append(int(k))
                    else:
                        colonIndex = line.find(':') + 2
                        buffer1 = line[colonIndex:]
                        self.misc = buffer1
                    i += 1


    def read_chromesome(self, endo):
        print('read chromo')
        self.chromesomeList.clear()
        tempList = []
        fileStream = gzip.open(self.fileName, 'r')
        i = 0
        first = True
        for line in fileStream:
            if i < 3:
                i += 1
            else:
                line = str(line)
                line = line.strip("'b")
                line = line[:len(line) - 4]
                if line[0] == '>':
                    if first == True:
                        tempList.append(line)
                        first = False
                    else:
                        self.chromesomeList.append(tempList)
                        tempList = []
                        tempList.append(line)
                elif line == 'REPEATS':
                    self.chromesomeList.append(tempList)
                    break
                else:
                    tempList.append(line)


    def read_repeats(self, endoChoice):
        print('read repeats')
        fileStream = gzip.open(self.fileName, 'r')
        repeats = False
        for line in fileStream:
            line = str(line)
            line = line.strip("'b")
            line = line[:len(line) - 4]
            if line == 'REPEATS':
                repeats = True
            elif repeats == True:
                line = line.split(r"\t")
                if len(line) == 1:
                    seed = line[0]
                    self.repeats[seed] = 0
                    self.seeds[seed] = []
                    self.dec_tup_data[seed] = []
                else:
                    line = line[:-1]
                    for item in line:
                        self.repeats[seed] += 1
                        item = item.split(',')
                        self.seeds[seed].append(item)
                        temp = item[2]
                        if '+' in item[2]:
                            temp = temp.split('+')
                            sign = '+'
                        else:
                            temp = temp.split('-')
                            sign = '-'
                        full_seq = temp[0] + seed
                        pam = temp[1]
                        tup = (item[1], full_seq, pam, item[3], sign, endoChoice)
                        self.dec_tup_data[seed].append(tup)
                        self.multiSum += int(item[3])
                        self.multiCount += 1


    def get_chromesome_names(self):
        print('read chromo names')
        self.chromesomesSelectedList.clear()
        with gzip.open(self.fileName,'r') as f:
            i = 0
            for line in f:
                line = str(line)
                line = line.strip("'b")
                line = line[:len(line) - 4]
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


    def popParser(self, cspr_file, endoChoice):
        print('pop parser')
        self.popData.clear()
        referenceList = list()
        repeats = False
        with gzip.open(cspr_file, 'r') as f:
            i = 0
            for line in f:
                line = str(line)
                line = line.strip("'b")
                line = line[:len(line) - 4]
                if i == 0:
                    genomeLine = line
                    genomeLine = genomeLine.split(',')
                    retNumber = int(genomeLine[len(genomeLine) - 1])
                elif i == 2:
                    misc_line = line
                    colonIndex = misc_line.find(':') + 2
                    usefulData = misc_line[colonIndex:]
                    usefulData = usefulData.split('|')
                    usefulData.pop()
                    j = 0
                    while j < len(usefulData):
                        temp = usefulData[j].split(',')
                        referenceList.append((temp[0], temp[1]))
                        j += 1

                elif line == 'REPEATS':
                    repeats = True
                elif repeats == True:
                    line = line.split(r"\t")
                    if len(line) == 1:
                        seed = line[0]
                        self.popData[seed] = []
                    else:
                        line = line[:-1]
                        for item in line:
                            commaIndex = item.find(',')
                            chrom = item[:commaIndex]
                            orgName = referenceList[int(chrom) - 1][0]
                            item = item.split(',')
                            temp = item[2]
                            if '+' in item[2]:
                                temp = temp.split('+')
                                sign = '+'
                            else:
                                temp = temp.split('-')
                                sign = '-'
                            full_seq = temp[0] + seed
                            pam = temp[1]
                            tup = (orgName, chrom, item[1], full_seq, pam, item[3], sign, endoChoice)
                            self.popData[seed].append(tup)
                i += 1

        return retNumber, referenceList


    def MT_parser(self, endo):
        self.seeds = {}
        self.repeats = {}
        self.dec_tup_data = {}
        self.chromesomeList = []
        self.karystatsList = []
        with gzip.open(self.fileName, 'r') as f:

            #first lines
            print('first lines')
            i = 0
            for line in f:
                if i > 2:
                    break
                else:
                    line = str(line)
                    line = line.strip("'b")
                    line = line[:len(line) - 5]
                    if i == 0:
                        colonIndex = line.find(':') + 2
                        buffer1 = line[colonIndex:]
                        self.genome = buffer1
                    elif i == 1:
                        colonIndex = line.find(':') + 2
                        k_data = line[colonIndex:]
                        k_data = k_data.split(',')
                        for k in k_data:
                            self.karystatsList.append(int(k))
                    else:
                        colonIndex = line.find(':') + 2
                        buffer1 = line[colonIndex:]
                        self.misc = buffer1
                        break
                    i += 1

            #read chrommosomes
            print('chromosomes')
            tempList = []
            first = True
            for line in f:
                line = str(line)
                line = line.strip("'b")
                line = line[:len(line) - 4]
                if line[0] == '>':
                    if first == True:
                        tempList.append(line)
                        first = False
                    else:
                        self.chromesomeList.append(tempList)
                        tempList = []
                        tempList.append(line)
                elif line == 'REPEATS':
                    self.chromesomeList.append(tempList)
                    break
                else:
                    tempList.append(line)

            print('repeats')
            for line in f:
                line = str(line)
                line = line.strip("'b")
                line = line[:len(line) - 4]
                line = line.split(r"\t")
                if len(line) == 1:
                    seed = line[0]
                    self.repeats[seed] = 0
                    self.seeds[seed] = []
                    self.dec_tup_data[seed] = []
                else:
                    line = line[:-1]
                    for item in line:
                        self.repeats[seed] += 1
                        item = item.split(',')
                        self.seeds[seed].append(item)
                        temp = item[2]
                        if '+' in item[2]:
                            temp = temp.split('+')
                            sign = '+'
                        else:
                            temp = temp.split('-')
                            sign = '-'
                        full_seq = temp[0] + seed
                        pam = temp[1]
                        tup = (item[1], full_seq, pam, item[3], sign, endo)
                        self.dec_tup_data[seed].append(tup)
                        self.multiSum += int(item[3])
                        self.multiCount += 1

            print('done parsing')


    def read_targets(self, genename, pos_tuple, endo):
        i = 0
        retList = []
        print(pos_tuple)
        print(self.fileName)
        header = False
        with gzip.open(self.fileName, 'r') as f:
            for line in f:
                line = str(line)
                line = line.strip("'b")
                line = line[:len(line) - 2]
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
                            retList.append((line[0], line[1], line[2], line[3], strand))
                        elif int(line[0]) >= int(pos_tuple[2]):
                            break
                    elif line == 'REPEATS':
                        break
                i += 1
        return retList


    def uniq_seq_count(self):
        print('uniq seq count')
        self.unique_targets = 0
        for chromo in self.chromesomeList:
            for data in chromo:
                if len(data) == 6:
                    self.unique_targets += 1
        return self.unique_targets
