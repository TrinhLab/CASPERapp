###############################################################################
# This is the Annotation Parser File
# INPUTS: inputs are the annotation files to parse, either txt or gff3
# OUTPUTS: the outputs are data structures that store the parsed data
################################################################################

import gffutils
import GlobalSettings
from Bio import SeqIO

class Annotation_Parser:
    def __init__(self):
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
    #####################END OF INIT FUNCTION

    ########################################
    # this function parses GBFF files
    # it stores the data in 2 dictionaries, a parallel one and a regular one
    # ONLY TO BE USED WITH GBFF FILES, NOTHING ELSE
    # testing: has been tested by Josh, and seems to work. Needs further testing by Brian
    ########################################
    def gbff_parse(self):
        # variables used
        self.reg_dict.clear()
        self.para_dict.clear()
        prevFirstIndex = ""
        indexNumber = 0
        strandChar = ""
        currentLocusTag = ""
        para_dict_key_string = ""
        values = list()

        # get all of the record or sections
        gb_record = SeqIO.parse(self.annotationFileName, 'genbank')

        # for each section
        for record in gb_record:
            indexNumber += 1
            # for each feature in that section
            for feature in record.features:
                # only change the locus_tag and update the para_dict if the feature is a gene
                if feature.type == "gene":
                    # once the locus tag changes, append it to the para_dict
                    if para_dict_key_string != "":
                        if para_dict_key_string not in self.para_dict:
                            self.para_dict[para_dict_key_string] = list()
                            self.para_dict[para_dict_key_string].append(currentLocusTag)
                        else:
                            if currentLocusTag not in self.para_dict[para_dict_key_string]:
                                self.para_dict[para_dict_key_string].append(currentLocusTag)
                        para_dict_key_string = ""

                    currentLocusTag = feature.qualifiers['locus_tag'][0]

                    # check to see if the strand is + or -
                    if feature.location.strand == -1:
                        strandChar = '-'
                    else:
                        strandChar = '+'

                    # update that one's values
                    values = [record.id, indexNumber, feature.type, int(feature.location.start) - 1,
                              int(feature.location.end), strandChar]

                    # insert
                    if currentLocusTag not in self.reg_dict:
                        self.reg_dict[currentLocusTag] = list()
                        self.reg_dict[currentLocusTag].append(values)
                    else:
                        self.reg_dict[currentLocusTag].append(values)

                # if it's not a gene, skip rep_orgin, telomere, and source
                elif feature.type != "rep_origin" and feature.type != "telomere" and feature.type != "source" and feature.type != 'assembly_gap':
                    # get the data for the normal dictionary and store it
                    if feature.location.strand == -1:
                        strandChar = '-'
                    else:
                        strandChar = '+'
                    values = [record.id, indexNumber, feature.type, int(feature.location.start) - 1,
                              int(feature.location.end), strandChar]
                    # make sure it isn't a duplicate
                    if values not in self.reg_dict[currentLocusTag]:
                        self.reg_dict[currentLocusTag].append(values)
                    # now get the para_dict's data

                    # check for the note section
                    if 'note' in feature.qualifiers:
                        if para_dict_key_string == "":
                            para_dict_key_string = feature.qualifiers['note'][0]
                        else:
                            para_dict_key_string = para_dict_key_string + ";" + feature.qualifiers['note'][0]
                    # check for the product section
                    if 'product' in feature.qualifiers:
                        if para_dict_key_string == "":
                            para_dict_key_string = feature.qualifiers['product'][0]
                        else:
                            para_dict_key_string = para_dict_key_string + ";" + feature.qualifiers['product'][0]

        # not sure if this number is correct, yet
        self.max_chrom = indexNumber

    ############################################
    # This function parses gff files and stores them in a dictionary
    # It also creates a parallel dictionary to use in searching
    # Precondition: ONLY TO BE USED WITH GFF FILES
    ############################################
    def gff_parse(self):
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
            if feature.featuretype == "gene":

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

                tempList = [feature.seqid, indexNumber, feature.featuretype, feature.start - 1, feature.end,
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
                    tempList = [child.seqid, indexNumber, child.featuretype, child.start - 1, child.end, child.strand]

                    # only insert it if it hasn't been inserted before
                    if tempList not in self.reg_dict[currentLocusTag]:
                        self.reg_dict[currentLocusTag].append(tempList)

            # now go through the other ones which are not region
            elif feature.featuretype != "region" and feature.featuretype != "telomere" and feature.featuretype != "origin_of_replication":
                tempList = [feature.seqid, indexNumber, feature.featuretype, feature.start - 1, feature.end,
                            feature.strand]

                # only insert if it hasn't been inserted before
                if tempList not in self.reg_dict[currentLocusTag]:
                    self.reg_dict[currentLocusTag].append(tempList)

                    # now same as above, go through the children again
                    for child in db.children(feature.id, level=None, featuretype=None, order_by=None, reverse=False,
                                             limit=None, completely_within=False):
                        tempList = [child.seqid, indexNumber, child.featuretype, child.start - 1, child.end,
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
############################END OF gff_parse

    ############################################
    # This function parses txt files and stores them in a dictionary
    # It also creates a parallel dictionary to use in searching
    # Precondition: ONLY TO BE USED WITH TXT FILES
    ############################################
    def txt_parse(self):
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
                    values = [splitLine[6], indexNumber, splitLine[0], int(splitLine[7]) - 1, int(splitLine[8]), splitLine[9]]

                    if currentLocusTag not in self.reg_dict:
                        self.reg_dict[currentLocusTag] = [values]
                    elif currentLocusTag in self.reg_dict:
                        self.reg_dict[currentLocusTag].append(values)

                # if no locus_tag, parse on product_accession, use the product_accession as the key for the dict
                elif not self.txtLocusTag:
                    currentLocusTag = splitLine[10]
                    values = [splitLine[6], indexNumber, splitLine[0], int(splitLine[7]) - 1, int(splitLine[8]), splitLine[9]]

                    if currentLocusTag not in self.reg_dict:
                        self.reg_dict[currentLocusTag] = [values]
                    elif currentLocusTag in self.reg_dict:
                        self.reg_dict[currentLocusTag].append(values)

                # set the parallel dictionary's key string
                para_dict_key_string = splitLine[13] + ";" + splitLine[12] + ";" + splitLine[14]

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
        self.max_chrom = indexNumber
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
        elif "gbff" in self.annotationFileName:
            self.gbff_parse()
        # return -1 to throw the error window in main
        else:
            return -1

####################END OF find_which_file_version
