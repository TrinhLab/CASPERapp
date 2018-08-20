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
