###############################################################################
# This is the Annotation Parser File
# INPUTS: inputs are the annotation files to parse. Currently, only gbff is supported.
# OUTPUTS: the outputs are data structures that store the parsed data
################################################################################

from PyQt5 import QtWidgets
import gffutils
import GlobalSettings
import os
from Bio import SeqIO
import traceback

#global logger
logger = GlobalSettings.logger

class Annotation_Parser:
    def __init__(self):
        try:
            #variables to use
            self.annotationFileName = "" #this is the variable that holds the filename itself
            self.txtLocusTag = False
            self.isGff = False
            self.isTxt = False
            self.max_chrom = 0

            #dictionary used for finding the genes in a txt annotation file
            #key: locus_tag
            #value: List of lists
            #   essentially its all based on locus tag. So the key is the locus tag, and its data is:
            #       [genomic accession, int, start, end, +\-]
            self.reg_dict = dict()

            #parallel dictionary used for the txt annotaion file
            #key: name + symbol (space in between each word)
            #value: locus_tag (indexes dict)
            self.para_dict = dict()
            
            #list of tuples containing (chromosome/scaffold # {int}, Feature matching search criteria {SeqFeature Object})
            self.results_list = list()

        except Exception as e:
            logger.critical("Error initializing Annotation_Parser class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    ### This function takes a list of lists and flattens it into a single list. Useful when dealing with a list of lists where the nested lists only have 1 entry.
    def flatten_list(self,t):
        return [item.lower() for sublist in t for item in sublist]

    ### This function finds how many chromosomes are within the selcted annotation file and returns the value 
    def get_max_chrom(self):
        parser = SeqIO.parse(self.annotationFileName, 'genbank') # Initialize parser (iterator) for each query
        for i, record in enumerate(parser):
            max_chrom = i+1
        return max_chrom

    def get_sequence_info(self, query):
        try:
            self.results_list.clear()
            parser = SeqIO.parse(self.annotationFileName, 'genbank') # Initialize parser (iterator) for each query
            for j,record in enumerate(parser): # Each record corresponds to a chromosome/scaffold in the FNA/FASTA file
                tmp = str(record.seq).find(query)
                if tmp != -1: # If match is found
                    return (j+1,tmp+1,tmp+len(query)) # Chromosome number, start index, stop index
                else:
                    tmp = str(record.seq.reverse_complement()).find(query) # Check the reverse complement now
                    if tmp != -1: # If match is found
                        return (j+1,tmp-len(query),tmp-1) # Chromosome number, start index, stop index
                    else:
                        continue
            return False
        
        except Exception as e:
            logger.critical("Error in get_sequence_info() in annotation parser.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)
    


 
    ### The workhorse function of AnnotationParser, this searches the annotation file for the user's search and returns features matching the description.
    def genbank_search(self, queries, same_search):
        index_number = 0
        try:
            if same_search: # If searching for the same thing, just return the results from last time
                return self.results_list
            else:
                self.results_list.clear()
                for i, query in enumerate(queries):
                    parser = SeqIO.parse(self.annotationFileName, 'genbank') # Initialize parser (iterator) for each query
                    for j,record in enumerate(parser): # Each record corresponds to a chromosome/scaffold in the FNA/FASTA file
                        if i == 0:
                            index_number += 1
                            for feature in record.features: # Each feature corresponds to a gene, tRNA, rep_origin, etc. in the given record (chromosome/scaffold)
                                if "translation" in feature.qualifiers:
                                    if query.lower() in " ".join(self.flatten_list(feature.qualifiers.values())[:-1]) and feature.type != "source" and feature.type != "gene": # If search matches the feature's qualifiers somewhere, save it
                                        self.results_list.append((j+1,feature))
                                    else: # If search not in the feature's qualifiers, move to the next feature
                                        continue
                                else:
                                    if query.lower() in " ".join(self.flatten_list(feature.qualifiers.values())) and feature.type != "source" and feature.type != "gene": # If search matches the feature's qualifiers somewhere, save it
                                        self.results_list.append((j+1,feature))
                                    else: # If search not in the feature's qualifiers, move to the next feature
                                        continue
                            self.max_chrom = index_number # Counts the number of chromosomes/scaffolds in the organism (only do this once, even if there are multiple queries)
                        else:
                            for feature in record.features:
                                if "translation" in feature.qualifiers:
                                    if query.lower() in " ".join(self.flatten_list(feature.qualifiers.values())[:-1]) and feature.type != "source" and feature.type != "gene": # If search matches the feature's qualifiers somewhere, save it
                                        self.results_list.append((j+1,feature))
                                    else: # If search not in the feature's qualifiers, move to the next feature
                                        continue
                                else:
                                    if query.lower() in " ".join(self.flatten_list(feature.qualifiers.values())) and feature.type != "source" and feature.type != "gene": # If search matches the feature's qualifiers somewhere, save it
                                        self.results_list.append((j+1,feature))
                                    else: # If search not in the feature's qualifiers, move to the next feature
                                        continue
                return self.results_list

        except Exception as e:
            logger.critical("Error in genbank_search() in annotation parser.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)
    



    # This function parses gff files and stores them in a dictionary
    # It also creates a parallel dictionary to use in searching
    # Precondition: ONLY TO BE USED WITH GFF FILES
    def gff_parse(self):
        try:
            self.reg_dict.clear()
            self.para_dict.clear()
            prevFirstIndex = ""
            indexNumber = 1
            fileStream = open(self.annotationFileName)
            data_base_file_name = GlobalSettings.CSPR_DB + "/" + "gff_database.db"

            # temp list will be the following each time it is put into the dictionary:
            # [Sequence ID (genomic accession or scaffold), the index number itself, the feature type (cds, gene, mrna), the start(-1), end, and the strand]
            tempList = list()
            currentLocusTag = ""
            para_dict_key_string = ""

            # initialize the data base (this is what parses it for me)
            print("Intializing the data base")
            db = gffutils.create_db(self.annotationFileName, dbfn=data_base_file_name, force=True, keep_order=True,
                                    merge_strategy='merge', sort_attribute_values=True)
            print("Finished intializing")

            # call the feature version of that data base now
            db = gffutils.FeatureDB(data_base_file_name, keep_order=True)

            # now we go through that data base and get the data we want
            for feature in db.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False,
                                           completely_within=False):
                # if the genomic accession/scaffold/chromseome changes, update the indexNumber
                if prevFirstIndex != feature.seqid and prevFirstIndex != "":
                    indexNumber += 1
                # if we find a new gene, update the locus_tag/name
                if feature.featuretype == "gene" or feature.featuretype == 'pseudogene':

                    # check and see if locus tag is in the attributes, go on the Name if locus_tag is not in there
                    if 'locus_tag' in feature.attributes:
                        currentLocusTag = feature.attributes['locus_tag'][0]
                    else:
                        currentLocusTag = feature.attributes["Name"][0]

                    # once the locus tag changes, append it to the para_dict
                    if para_dict_key_string != "":
                        if para_dict_key_string not in self.para_dict:
                            self.para_dict[para_dict_key_string] = list()
                            self.para_dict[para_dict_key_string].append(currentLocusTag)
                        else:
                            if currentLocusTag not in self.para_dict[para_dict_key_string]:
                                self.para_dict[para_dict_key_string].append(currentLocusTag)
                        para_dict_key_string = ""

                    tempList = [currentLocusTag, indexNumber, feature.featuretype, feature.start - 1, feature.end,
                                feature.strand]

                    # insert that locus tag/name into the dictionary
                    if currentLocusTag not in self.reg_dict:
                        self.reg_dict[currentLocusTag] = []
                        self.reg_dict[currentLocusTag].append(tempList)
                    elif currentLocusTag in self.reg_dict:
                        self.reg_dict[currentLocusTag].append(tempList)

                    # go through each of this child's children
                    for child in db.children(feature.id, level=None, featuretype=None, order_by=None, reverse=False,
                                             limit=None, completely_within=False):
                        tempList = [currentLocusTag, indexNumber, child.featuretype, child.start - 1, child.end, child.strand]

                        # only insert it if it hasn't been inserted before
                        if tempList not in self.reg_dict[currentLocusTag]:
                            self.reg_dict[currentLocusTag].append(tempList)

                # now go through the other ones which are not region
                elif feature.featuretype != "region" and feature.featuretype != "telomere" and feature.featuretype != "origin_of_replication":
                    tempList = [currentLocusTag, indexNumber, feature.featuretype, feature.start - 1, feature.end,
                                feature.strand]

                    # only insert if it hasn't been inserted before
                    if tempList not in self.reg_dict[currentLocusTag]:
                        self.reg_dict[currentLocusTag].append(tempList)

                        # now same as above, go through the children again
                        for child in db.children(feature.id, level=None, featuretype=None, order_by=None, reverse=False,
                                                 limit=None, completely_within=False):
                            tempList = [currentLocusTag, indexNumber, child.featuretype, child.start - 1, child.end,
                                        child.strand]

                            if tempList not in self.reg_dict[currentLocusTag]:
                                self.reg_dict[currentLocusTag].append(tempList)

                # now we need to get the para_dict up and running
                # get the stuff out of the product part
                if 'product' in feature.attributes and feature.featuretype == "CDS":
                    if para_dict_key_string == "":
                        para_dict_key_string = feature.attributes['product'][0]
                    else:
                        para_dict_key_string = para_dict_key_string + ";" + feature.attributes['product'][0]
                # get the stuff out of the Note part
                if 'Note' in feature.attributes:
                    if para_dict_key_string == "":
                        para_dict_key_string = feature.attributes['Note'][0]
                    else:
                        para_dict_key_string = para_dict_key_string + ";" + feature.attributes['Note'][0]

                prevFirstIndex = feature.seqid
            self.max_chrom = indexNumber
        except Exception as e:
            logger.critical("Error in gff_parse() in annotation parser.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # This function parses txt files and stores them in a dictionary
    # It also creates a parallel dictionary to use in searching
    # Precondition: ONLY TO BE USED WITH TXT FILES
    def txt_parse(self):
        try:
            self.reg_dict.clear()
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
                    if(len(buffer) <= 2): # break out once we reach the end of the file
                        break

                    splitLine = buffer[:-1].split("\t")

                    # increment indexNumber when genomic access changes
                    if prevGenAccession != splitLine[6] and prevGenAccession != "":
                        indexNumber += 1

                    # if parsing on locus_tag, use the locus_tag as the key for the dict
                    if self.txtLocusTag:
                        currentLocusTag = splitLine[16]
                        values = [currentLocusTag, indexNumber, splitLine[0], int(splitLine[7]) - 1, int(splitLine[8]), splitLine[9]]

                        if currentLocusTag not in self.reg_dict:
                            self.reg_dict[currentLocusTag] = [values]
                        elif currentLocusTag in self.reg_dict:
                            self.reg_dict[currentLocusTag].append(values)

                    # if no locus_tag, parse on product_accession, use the product_accession as the key for the dict
                    elif not self.txtLocusTag:
                        currentLocusTag = splitLine[10]
                        values = [currentLocusTag, indexNumber, splitLine[0], int(splitLine[7]) - 1, int(splitLine[8]), splitLine[9]]

                        if currentLocusTag not in self.reg_dict:
                            self.reg_dict[currentLocusTag] = [values]
                        elif currentLocusTag in self.reg_dict:
                            self.reg_dict[currentLocusTag].append(values)

                    if splitLine[13] != '':
                        if para_dict_key_string == '':
                            para_dict_key_string = splitLine[13] + ';'
                        else:
                            para_dict_key_string = para_dict_key_string + splitLine[13] + ';'

                    # leaving this in for now, it's related accession
                    #if splitLine[12] != '':
                     #   if para_dict_key_string == '':
                    #        para_dict_key_string = splitLine[12] + ';'
                     #   else:
                     #       para_dict_key_string = para_dict_key_string + splitLine[12] + ';'


                    if splitLine[14] != '':
                        if para_dict_key_string == '':
                            para_dict_key_string = splitLine[14] + ';'
                        else:
                            para_dict_key_string = para_dict_key_string + splitLine[14] + ';'

                    para_dict_key_string = para_dict_key_string.replace(',', '')
                    # set the parallel dictionary's key string
                    #para_dict_key_string = splitLine[13] + ";" + splitLine[12] + ";" + splitLine[14]

                    # if the current line we're on has the data we want for the parellel dictionary, store it
                    if len(para_dict_key_string) > 3:
                        if para_dict_key_string[len(para_dict_key_string) - 1] == ';':
                            para_dict_key_string = para_dict_key_string[0:len(para_dict_key_string) - 1]

                        if para_dict_key_string not in self.para_dict: # make a new input into the dict
                            self.para_dict[para_dict_key_string] = [currentLocusTag]
                        elif para_dict_key_string in self.para_dict:
                            if currentLocusTag not in self.para_dict[para_dict_key_string]:
                                # only append it to the dict's list if it isn't currently in there
                                self.para_dict[para_dict_key_string].append(currentLocusTag)

                    para_dict_key_string = ""
                    prevGenAccession = splitLine[6]
            self.max_chrom = indexNumber
        except Exception as e:
            logger.critical("Error in txt_parse() in annotation parser.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # This function checks to see which file we are parsing
    # It also checks whether to parse based on locus_tag or product accession (txt files only)
    # Then it calls the respective parser functions used
    def find_which_file_version(self):
        try:
            if self.annotationFileName == "" or GlobalSettings.mainWindow.annotation_files.currentText() == "None":
                return -1
            if "gff" in self.annotationFileName:
                ### gff file support currently deprecated
                """
                self.isGff = True
                self.gff_parse()
                """
                print("Error: Wrong annotation file format")
                return -1
                
            elif "feature_table" in self.annotationFileName:
                ### feature table file support currently deprecated
                # now that we know it's a txt file and not a gff, check and see if we will be parsing by locus tag or
                # product accession
                """
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
                """
                print("Error: Wrong annotation file format")
                return -1
            elif "gbff" or "gbk" in self.annotationFileName:
                return "gbff"
            # return -1 to throw the error window in main
            else:
                return -1
        except Exception as e:
            logger.critical("Error in find_which_file_version() in annotation parser.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(GlobalSettings.mainWindow.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

