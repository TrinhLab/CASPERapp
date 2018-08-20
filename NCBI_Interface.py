"""Class Ncbi function: The instance of this class allows the user to download a genome (or genomes) from the NCBI database
 based on common biological search terms e.g. Saccharomyces cerevisiae or E. coli K-12 MG1655.
 The functionality should allow for the class to return both well-annotated complete genomes and less annotated genome
 assemblies if so desired by the user as specified by an input parameter.
 The output should be to write a .fasta file to a desired directory."""

from Bio import Entrez
import time

class Ncbi:

    def __init__(self, email):
        Entrez.email = email

    #input file_type: e.g. 'fasta'; database: e.g. 'nucleotide'; iden: e.g. CP010053.1
    def get_file(self, file_type, database, iden):
        handle = Entrez.esummary(db=database, id='1220112116')
        record = Entrez.read(handle)
        print(record)
        handle.close()

    def search_database(self, searchterm, database):
        handle = Entrez.esearch(db=database, term=searchterm)
        record = Entrez.read(handle)
        identify = record['IdList']
        print(identify)
        handle.close()
        allfetch = list()
        #for item in identify:
        item = identify[8]
        newhandle = Entrez.efetch(db=database, id=item, rettype = 'gb', retmode='text')
        mything = newhandle.read()
        print(item)
        print(mything)
        allfetch.append(mything)
        newhandle.close()
        time.sleep(0.5)


N = Ncbi('bmendoz1@vols.utk.edu')
# N.get_file('fasta', 'nucleotide', 'CP010053.1')
N.search_database('Saccharomyces cerevisiae[orgn] AND complete genome[title]', "nucleotide")
