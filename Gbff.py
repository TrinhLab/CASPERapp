"""Parsing GBFF, GFF, and tabular genome annotations.  GFF and tabular do not have the nucleotide sequence while GBFF
has everything"""

from Bio import SeqIO, Entrez
from bs4 import BeautifulSoup
import requests
from ftplib import FTP
import time

Entrez.email = "bmendoz1@vols.utk.edu"


class GBFF_Parse:

    def __init__(self, filename, NCBI=False):
        self.filename = filename

    def convert_to_fasta(self):
        count = SeqIO.convert(self.filename, "genbank", self.filename[:-5]+"fasta", "fasta")


class Assembly:

    def __init__(self, organism, database):
        handle = Entrez.esearch(db="assembly", retmax=20, term=organism)
        record = Entrez.read(handle)
        self.database = database
        myidlist = list()
        for id in record["IdList"]:
            myidlist.append(id)
            print(id)
        handle.close()

        """for ret in myidlist:
            handle = Entrez.esummary(db="assembly", id=ret)
            gca_rec = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            print(gca_rec)
            handle.close()
            time.sleep(0.5)"""

        url = "https://www.ncbi.nlm.nih.gov/assembly/" + "GCF_000146045.2"
        source = requests.get(url)
        plain_text = source.text
        soup = BeautifulSoup(plain_text, "html.parser")

        refseq_link = str(soup.find('a', text="Download the RefSeq assembly"))
        genbank_link = str(soup.find('a', text="Download the GenBank assembly"))
        refseq_link = refseq_link[refseq_link.find("=")+2: refseq_link.find(">")-1]
        genbank_link = genbank_link[genbank_link.find("=") + 2: genbank_link.find(">") - 1]
        # Checkpoint printing
        print(refseq_link)
        print(genbank_link)
        # then go to that webpage then inspect and download with url service the gbff file
        if self.database == "GenBank":
            database_url = genbank_link
        else:
            database_url = refseq_link
        data_source = requests.get(database_url)
        soup = BeautifulSoup(data_source.text, "html.parser")
        for link in soup.find_all('a', {'class': 'icon file'}):
            print(link)

        # can save or throw away, save for those with a bunch of space or call the gff for annotation when needed


#G = GBFF_Parse("/Users/brianmendoza/Desktop/FileExamples/PantoeaYR343.gbff")
#G.convert_to_fasta()

#A = Assembly("pantoea[Organism]","RefSeq")

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd("/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/")

filename = "GCF_000146045.2_R64_genomic.fna.gz"

#local_filename = os.path.join(r"c:\myfolder", filename)
local_filename = "/Users/brianmendoza/Desktop/" + filename
lf = open(local_filename, "wb")
ftp.retrbinary("RETR " + filename, lf.write, 8*1024)
lf.close()

ftp.close()
