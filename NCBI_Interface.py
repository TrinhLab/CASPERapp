"""Class Ncbi function: The instance of this class allows the user to download a genome (or genomes) from the NCBI database
 based on common biological search terms e.g. Saccharomyces cerevisiae or E. coli K-12 MG1655.
 The functionality should allow for the class to return both well-annotated complete genomes and less annotated genome
 assemblies if so desired by the user as specified by an input parameter.
 The output should be to write a .fasta file to a desired directory."""

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import time
import os
import io


class Ncbi:

    def __init__(self, email):
        Entrez.email = email

    #input file_type: e.g. 'fasta'; database: e.g. 'nucleotide'; iden: e.g. CP010053.1
    def get_file(self, database, iden):
        handle = Entrez.esummary(db=database, id=iden)
        record = Entrez.read(handle)
        print(record)
        handle.close()

    def search_database(self, searchterm, database, path):
        handle = Entrez.esearch(db=database, term=searchterm)
        record = Entrez.read(handle)
        identify = record['IdList']                                            #get the id from database
        print(identify)
        handle.close()
        allfetch = list()
        item = identify[8]
        newhandle = Entrez.efetch(db=database, id=item, rettype = 'fasta', retmode='text')  #use id to get information from database in fasta format
        mything = newhandle.read()                                                          #get informaion as string
        if not os.path.exists(path):                                                        #check to see if the directory path exist
            print("That path can not be found.")
        else:                                                                               #if it does then parse through data to get id and sequence as their own strings
            seq =[]
            seq_full = ""
            buf = io.StringIO(mything)                      #put string into buffer
            seq_ID  = buf.readline()[1:]                    #get first id
            while True:

                line = buf.readline()
                if len(line)== 0:                           #if end of file then imput last sequence and end
                    record = SeqRecord(Seq(seq_full, IUPAC.protein), id=seq_ID)
                    seq.append(record)
                    break

                if line[0] == '>':                          #if begining of another sequence then record last one as done and save new id
                    record = SeqRecord(Seq(seq_full, IUPAC.protein), id=seq_ID)
                    seq.append(record)
                    seq_ID = line[1:]
                    seq_full = ""
                else:                                       #otherwise concatinate sequence string without endline character
                    seq_full = seq_full+line[:len(line)-1]
            with open(os.path.join(path,"GeneSequence.fasta"), 'w+') as temp_file:   #create file and add write the sequences to it
                SeqIO.write(seq,temp_file,"fasta")



        #  print(item)
        #print(mything)

        #N.to_file(mything)
        allfetch.append(mything)
        newhandle.close()
        time.sleep(0.5)








N = Ncbi('bmendoz1@vols.utk.edu')

#N.get_file('nucleotide', 'CP010053.1')
N.search_database('Saccharomyces cerevisiae[orgn] AND complete genome[title]', "nucleotide","C:\\Users\\Greg Cantrall\\Documents")
