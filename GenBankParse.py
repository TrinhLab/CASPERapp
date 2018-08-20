__author__ = 'brianmendoza'

from Bio import Entrez, SeqIO
import webbrowser
import re
import os


class GenBankFile:

    def __init__(self, organism):
        Entrez.email = "bmendoz1@vols.utk.edu"
        self.directory = "/Users/brianmendoza/Desktop/GenBank_files/"
        self.org = organism

    def setOrg(self, org):
        self.org = org

    def setDirectory(self, path, org):
        self.directory = path
        self.setOrg(org)

    def convertToFasta(self):
        orgfile = self.directory + self.org + ".gbff"
        output = "/Users/brianmendoza/Desktop/GenBank_files/FASTAs/" + self.org + ".fna"
        SeqIO.convert(orgfile, "genbank", output, "fasta")

    def parseAnnotation(self):
        gb_file = self.directory + self.org + ".gbff"
        records = SeqIO.parse(open(gb_file,"r"), "genbank")

        # create table for multi-targeting reference
        table = {}
        count = 0
        for record in records:
            count += 1
            chrmnumber = str(count)
            table[chrmnumber] = []
            for feature in record.features:
                if feature.type == 'CDS': # hopefully gene and CDS are the same

                    # getting the location...
                    loc = str(feature.location)
                    out = re.findall(r"[\d]+", loc)
                    start = out[0]
                    end = out[1]
                    if len(out) > 2:  # to account for "joined" domains
                        end = out[3]

                    # locus_tag and product...
                    if 'locus_tag' in feature.qualifiers:
                        ltag = feature.qualifiers['locus_tag']
                    elif 'gene' in feature.qualifiers:
                        ltag = feature.qualifiers['gene']
                    if 'product' not in feature.qualifiers:
                        prod = feature.qualifiers['note']
                    else:
                        prod = feature.qualifiers['product']
                    # adding it all up...
                    tup = (start, end, ltag, prod)
                    table[chrmnumber].append(tup)
        return table

    def getChromSequence(self, index):
        gb_file = self.directory + self.org + ".gbff"
        records = SeqIO.parse(open(gb_file,"r"), "genbank")
        count = 0
        for record in records:
            count += 1
            if count == index:
                cstr = record.seq
                return cstr

