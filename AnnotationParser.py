###############################################################################
# This is the Annotation Parser File
# INPUTS: inputs are the annotation files to parse, either txt or gff3
# OUTPUTS: the outputs are data structures that store the parsed data
################################################################################

from Algorithms import SeqTranslate
from operator import itemgetter
import os,sys
import re


class Annotation_Parser:
    def __init__(self):
        #variables to use
        self.annotationFileName = "" #this is the variable that holds the filename itself
        self.txtLocusTag = False
        self.isGff = False
        self.isTxt = False

        #dictionary used for finding the genes in a txt annotation file
        #key: locus_tag
        #value: List of lists
        #   essentially its all based on locus tag. So the key is the locus tag, and its data is:
        #       [genomic accession, int, start, end, +\-]
        self.dict = dict()

        #parallel dictionary used for the txt annotaion file
        #key: name + symbol (space in between each word)
        #value: locus_tag (indexes dict)
        self.para_dict = dict()
    #####################END OF INIT FUNCTION

    ############################################
    # This function parses gff files and stores them in a dictionary
    # It also creates a parallel dictionary to use in searching
    # Precondition: ONLY TO BE USED WITH GFF FILES
    ############################################
    def gff_parse(self):
        self.dict.clear()
        prevFirstIndex = ""
        indexNumber = 1
        fileStream = open(self.annotationFileName)
        buffer = ""
        currentLocusTag = ""
        para_dict_key_string = ""

        while(True): #this loop breaks out once the file is empty
            buffer = fileStream.readline()

            if(buffer.startswith("#")):
                #print(buffer)
                continue
            else:
                if(len(buffer) <= 0):
                    break
                splitLine = buffer[:-1].split("\t")

                if prevFirstIndex != splitLine[0] and prevFirstIndex != "":
                    indexNumber+=1

                keyWords = splitLine[8]
                keyWordsSplit = keyWords.split(";")

                #get the original dict's data
                for i in range(len(keyWordsSplit)):
                    # only get the lines that have Name, locus_tag, or standard_name,
                    # however, Name and locus_tag are ideal
                    if keyWordsSplit[i].startswith("Name=") or keyWordsSplit[i].startswith("locus_tag="):
                        currentLocusTag = (keyWordsSplit[i][keyWordsSplit[i].find("=") + 1:])
                        values = [splitLine[0], indexNumber, splitLine[2], int(splitLine[3]) - 1, int(splitLine[4]), splitLine[6]]

                        if currentLocusTag not in self.dict:
                            self.dict[currentLocusTag] = [values]
                        elif currentLocusTag in self.dict:
                            self.dict[currentLocusTag].append(values)

                        break
                    elif keyWordsSplit[i].startswith("standard_name="):
                        currentLocusTag = (keyWordsSplit[i][keyWordsSplit[i].find("=") +1:-2])
                        values = [splitLine[0], indexNumber, splitLine[2], int(splitLine[3]) - 1, int(splitLine[4]),
                                  splitLine[6]]
                        if currentLocusTag not in self.dict:
                            self.dict[currentLocusTag] = [values]
                        elif currentLocusTag in self.dict:
                            self.dict[currentLocusTag].append(values)


                # now get the para dict's data
                for j in range(len(keyWordsSplit)):
                    if keyWordsSplit[j].startswith("Note=") or keyWordsSplit[j].startswith("product="):
                        para_dict_key_string = para_dict_key_string + " " + (
                            keyWordsSplit[j][keyWordsSplit[j].find("=") + 1:])

                # make sure that the string actually has data in it
                if(para_dict_key_string != ""):
                    if para_dict_key_string not in self.para_dict: # make a new input into the dict
                        self.para_dict[para_dict_key_string] = [currentLocusTag]
                    elif para_dict_key_string in self.para_dict:
                        if currentLocusTag not in self.para_dict[para_dict_key_string]:
                            # only append it to the dict's list if it isn't currently in there
                            self.para_dict[para_dict_key_string].append(currentLocusTag)
                para_dict_key_string = ""
                prevFirstIndex = splitLine[0]
############################END OF gff_parse

    ############################################
    # This function parses txt files and stores them in a dictionary
    # It also creates a parallel dictionary to use in searching
    # Precondition: ONLY TO BE USED WITH TXT FILES
    ############################################
    def txt_parse(self):
        self.dict.clear()
        prevGenAccession = ""
        indexNumber = 1
        fileStream = open(self.annotationFileName)
        buffer = ""
        currentLocusTag = ""
        para_dict_key_string = ""

        while(True): # this loop breaks out when buffer string is empty
            buffer = fileStream.readline()

            if(buffer.startswith("#")): #skip lines that start with #
                continue
            else:
                if(len(buffer) <= 0): # break out once we reach the end of the file
                    break

                splitLine = buffer[:-1].split("\t")

                # increment indexNumber when genomic access changes
                if prevGenAccession != splitLine[6] and prevGenAccession != "":
                    indexNumber += 1

                # if parsing on locus_tag, use the locus_tag as the key for the dict
                if self.txtLocusTag:
                    currentLocusTag = splitLine[16]
                    values = [splitLine[6], indexNumber, splitLine[0], int(splitLine[7]) - 1, int(splitLine[8]), splitLine[9]]

                    if currentLocusTag not in self.dict:
                        self.dict[currentLocusTag] = [values]
                    elif currentLocusTag in self.dict:
                        self.dict[currentLocusTag].append(values)

                # if no locus_tag, parse on product_accession, use the product_accession as the key for the dict
                elif not self.txtLocusTag:
                    currentLocusTag = splitLine[10]
                    values = [splitLine[6], indexNumber, splitLine[0], int(splitLine[7]) - 1, int(splitLine[8]), splitLine[9]]

                    if currentLocusTag not in self.dict:
                        self.dict[currentLocusTag] = [values]
                    elif currentLocusTag in self.dict:
                        self.dict[currentLocusTag].append(values)

                # set the parallel dictionary's key string
                para_dict_key_string = para_dict_key_string + " " + splitLine[13] + " " + splitLine[12] + " " + splitLine[14]

                # if the current line we're on has the data we want for the parellel dictionary, store it
                if len(para_dict_key_string) > 3:
                    if para_dict_key_string not in self.para_dict: # make a new input into the dict
                        self.para_dict[para_dict_key_string] = [currentLocusTag]
                    elif para_dict_key_string in self.para_dict:
                        if currentLocusTag not in self.para_dict[para_dict_key_string]:
                            # only append it to the dict's list if it isn't currently in there
                            self.para_dict[para_dict_key_string].append(currentLocusTag)

                para_dict_key_string = ""
                prevGenAccession = splitLine[6]
#########################END OF txt_parse

    ############################################
    # This function checks to see which file we are parsing
    # It also checks whether to parse based on locus_tag or product accession (txt files only)
    # Then it calls the respective parser functions used
    ############################################
    def find_which_file_version(self):
        if self.annotationFileName == "" :
            print("Error: No annotation file given")
            return
        if "gff" in self.annotationFileName:
            self.isGff = True
            self.gff_parse()
        elif "feature_table" in self.annotationFileName:
            # now that we know it's a txt file and not a gff, check and see if we will be parsing by locus tag or
            # product accession
            fileStream = open(self.annotationFileName)

            #skip all of the lines that start with #
            buf = fileStream.readline()
            while buf.startswith("#"):
                buf = fileStream.readline()

            # split it and see if the locus tag spot has data in it
            split = buf.split("\t")
            if split[16] != "": # if it does, we are parsing based on locus_tag
                self.txtLocusTag = True
            elif split[16] == "": # if not, we are parsing based on product accession
                self.txtLocusTag = False
            fileStream.close()
            self.isTxt = True
            self.txt_parse()
        else:
            print("Error: We cannot parse the file you have given")

####################END OF find_which_file_version

#below code is for testing purposes
"""
if __name__ == '__main__':
    myParser = Annotation_Parser()

    #testing code with diffrent file names
    #myParser.annotationFileName = "Annotation_Files/GCA_000146045.2_R64_feature_table.txt"
    #myParser.annotationFileName = "Annotation_Files/GCA_003004805.1_ASM300480v1_feature_table.txt"
    #myParser.annotationFileName = "GCA_900537225.1_YALIH222_genomic.gff"
    #myParser.annotationFileName = "Kfedtschenkoi_382_v1.1.gene.gff3"
    myParser.annotationFileName = "GCA_900537225.1_YALIH222_feature_table.txt"
    myParser.find_which_file_version()

    # this prints out the parallel dictionary
    for item in myParser.para_dict:
        print(item)
        for i in range(len(myParser.para_dict[item])):
            print("\t", myParser.para_dict[item][i])

    # this prints out the regualr dictionary
    #for item in myParser.dict:
     #   print(item)
      #  for i in range(len(myParser.dict[item])):
       #     print("\t", myParser.dict[item][i])

    # prints out the sizes of both of them
    print("Regular dict size: ", len(myParser.dict))
    print("Para dict size: ", len(myParser.para_dict))
"""