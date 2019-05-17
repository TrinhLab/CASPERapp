from Algorithms import SeqTranslate

class CSPRparser:
    #default ctor: currently just sets the file name and initializes all of the variables I will be using
    def __init__(self, inputFileName):
        super(CSPRparser, self).__init__()

        # variables used in this class
        self.multiSum = 0
        self.multiCount = 0
        self.seqTrans = SeqTranslate()
        self.chromesomeList = list()  # list of a list for the chromesomes
        self.karystatsList = list()  # list of (ints) of the karyStats (whatever those are) to be used for the get_chrom_length function
        self.genome = ""  # genome name
        self.misc = ""  # anything from the misc line
        self.repeats = {}  # list of all the repeats. may need to be change, still not that sure about what the repeats are
        self.seeds = {}

        #file path variable
        self.fileName = ""  # filename itself, intialized in the default ctor
        self.fileName = inputFileName

    #this function reads the first 3 lines of the file: also stores the karyStats in a list of ints
    def read_first_lines(self):
        fileStream = open(self.fileName, 'r')

        #read and parse the genome line
        self.genome = fileStream.readline()
        colonIndex = self.genome.find(':') + 2
        buffer1 = self.genome[colonIndex:]
        self.genome = buffer1

        #read and store the karystats line on its own, it is parsed down below
        buffer = fileStream.readline()

        #read and parse the misc line
        self.misc = fileStream.readline()
        colonIndex = self.misc.find(':') + 2
        buffer1 = self.misc[colonIndex:]
        self.misc = buffer1


        #now parse the karystats line
        #ignore the first bit of the string. only care about what's after the colon
        colonIndex = buffer.find(':') + 2

        #parse the line, store the numbers in the list
        for i in range(colonIndex, len(buffer)):
            bufferString1 = ""
            if buffer[i] == ',':
                bufferString1 = buffer[colonIndex:i]
                #print(bufferString1)
                colonIndex = i + 1
                self.karystatsList.append(int(bufferString1))

        #print(self.karystatsList)


#this function reads all of the chromosomes in the file
#stores the data into a list of lists. So the line starting with '>' is the first index of each sub list
    def read_chromesome(self):
        tempList = list()
        fileStream = open(self.fileName, 'r')

        #ignore the first 3 lines
        fileStream.readline()
        fileStream.readline()
        fileStream.readline()

        bufferString = fileStream.readline()
        while(True): #this loop breaks out when bufferString is REPEATS
            tempList.append(bufferString)

            if(bufferString == "REPEATS\n"):
                break
            bufferString = fileStream.readline()
            while(True): #this loop breaks out when bufferString[0] is >
                if(bufferString == "REPEATS\n"):
                    self.chromesomeList.append(tempList)
                    tempList = []
                    break

                elif(bufferString[0] == '>'): #if we get to the next chromesome, append the tempList, clear it, and break
                    self.chromesomeList.append(tempList)
                    tempList = []
                    break
                else: #else decompress the data, and append it to the list
                    #print(bufferString)
                    bufferString = self.seqTrans.decompress_csf_tuple(bufferString)
                    tempList.append(bufferString)
                    #print(bufferString)
                    bufferString = fileStream.readline()

    #this function reads just the repeats
    #it stores this data in 2 dictionaries
        #repeats dictionary is the number of dictionaries
        #seeds dictionary is each seed that is repeated
    #this function also stores the sum and count in the class itself as well
    def read_repeats(self):
        file_info = self.get_whole_file()
        split_info = file_info.split('\n')

        index = 0

        #skip first part of the file until it reaches the repeats part
        while True:
            if(split_info[index] == "REPEATS"):
                index = index + 1
                break
            else:
                index = index + 1

        #clear what is already in there
        self.repeats.clear()
        self.seeds.clear()

        self.multiSum = 0
        self.multiCount = 0

        #parse the info now and store it in the correct dictionaries
        #this essentially copies what is happening in Multitargeting.py
        while(index + 1 < len(split_info)):
            seed = self.seqTrans.decompress64(split_info[index])
            repeat =split_info[index + 1].split("\t")

            self.repeats[seed] = 0
            self.seeds[seed] = []

            for item in repeat:
                if item != "":
                    self.repeats[seed] += 1
                    sequence =item.split(',')
                    self.seeds[seed].append(sequence)
                    self.multiSum += self.seqTrans.decompress64(sequence[3])
                    self.multiCount += 1

            index = index + 2

    #this function just calls the above read functions
    def read_all(self):
        print("Reading First Lines.")
        self.read_first_lines()
        print("Reading Chromesomes.")
        self.read_chromesome()
        print("Reading Repeats.")
        self.read_repeats()

    def get_whole_file(self):
        fileStream = open(self.fileName)
        return(fileStream.read())

if __name__ == '__main__':
    parser = CSPRparser("../NewCSPRFile.cspr")
    parser.read_first_lines()

    print(sum(parser.karystatsList))

#below is testing code
    """
    print("Filename: " + parser.fileName)
    print("genome: " + parser.genome)
    print("Misc: " + parser.misc)
    print("KaryStats: \n\t" + str(parser.karystatsList))

    chromesomeFile = open("../ChromesomeList.txt", 'w')
    repeatsNumFile = open("../RepeatsNum.txt", 'w')
    seedsFile = open("../SeedsFile.txt", 'w')

    print('Writing chromesomes')
    for i in range(len(parser.chromesomeList)):
        chromesomeFile.write("List: " + str(i) + "\n")
        for j in range(len(parser.chromesomeList[i])):
            chromesomeFile.write("\t" + str(parser.chromesomeList[i][j]) + "\n")

    print('Writing Repeats')
    for repeat in parser.repeats:
        repeatsNumFile.write("Seed Number: " + str(repeat) + "\n")
        repeatsNumFile.write("\t Number of repeats: " + str(parser.repeats[repeat]) + "\n")

    print('Writing Seeds')
    for seed in parser.seeds:
        seedsFile.write("Seed Number: " + str(seed) + "\n")
        seedsFile.write("\t Repeat Sequences: " + str(parser.seeds[seed]) + "\n")

    chromesomeFile.close()
    repeatsNumFile.close()
    seedsFile.close()
    """