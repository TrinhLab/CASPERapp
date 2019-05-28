"""Parsing GBFF, GFF, and tabular genome annotations.  GFF and tabular do not have the nucleotide sequence while GBFF
has everything"""

from Bio import SeqIO, Entrez
from bs4 import BeautifulSoup
import requests
from ftplib import FTP
import time
import GlobalSettings

Entrez.email = "bmendoz1@vols.utk.edu"


class GBFF_Parse:

    def __init__(self, filename, NCBI=False):
        self.filename = filename

    def convert_to_fasta(self):
        count = SeqIO.convert(self.filename, "genbank", self.filename[:-5]+"fasta", "fasta")


class Assembly:

    def __init__(self):
        print('')
        #Moved the code from this to its own function so that I can return the database URL

    def getDataBaseURL(self, organism, database):
        handle = Entrez.esearch(db="assembly", retmax=20, term=organism)
        record = Entrez.read(handle)
        self.database = database
        myidlist = list()

        for id in record["IdList"]:
            myidlist.append(id)
            #print(id)
        handle.close()
        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(30)

        gca_rectList = list()
        for ret in myidlist:
            handle = Entrez.esummary(db="assembly", id=ret)
            gca_rec = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            #print(gca_rec)
            gca_rectList.append(gca_rec)
            handle.close()
            time.sleep(0.1)
        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(60)

        #is it ok that the URL is hard coded?
        if(gca_rectList):
            url = "https://www.ncbi.nlm.nih.gov/assembly/" + gca_rectList[0] + "/"
            source = requests.get(url)
            plain_text = source.text
            soup = BeautifulSoup(plain_text, "html.parser")

            refseq_link = str(soup.find('a', text="Download the RefSeq assembly"))
            genbank_link = str(soup.find('a', text="Download the GenBank assembly"))
            refseq_link = refseq_link[refseq_link.find("=") + 2: refseq_link.find(">") - 1]
            genbank_link = genbank_link[genbank_link.find("=") + 2: genbank_link.find(">") - 1]
            # Checkpoint printing
            #print(refseq_link)
            #print(genbank_link)
            # then go to that webpage then inspect and download with url service the gbff file
            if self.database == "GenBank":
                database_url = genbank_link
            else:
                database_url = refseq_link
            # data_source = requests.get(database_url)
            # soup = BeautifulSoup(data_source.text, "html.parser")
            for link in soup.find_all('a', {'class': 'icon file'}):
                print(link)

            #check the links and catch errors
            if(database == "RefSeq" and len(refseq_link) < 5):
                print("Error: No RefSeq file to download!")
                return ("Error: No RefSeq file to download!", list())
            elif(database == "GenBank" and len(genbank_link) < 5):
                print("Error: No GenBank file to download!")
                return ("Error: No GenBank file to download!", list())
            elif(len(genbank_link) < 5 and len(refseq_link) < 5):
                print("Error: No RefSeq or GenBank files to download")
                return ("Error: No RefSeq or GenBank files to download", list())
            else:
                return database_url, myidlist
        else:
            print('Error: no link found. Check your spelling')
            return ("No link found", list())

#Testing below
"""
#G = GBFF_Parse("/Users/brianmendoza/Desktop/FileExamples/PantoeaYR343.gbff")
#G.convert_to_fasta()

#What exacly is the pantoea[Organism]? user input?
A = Assembly("pantoea[Organism]","RefSeq")

#hard coded currently: pantoea will be from user as well as the GenBank or RefSeq
print("Getting Database Link")
database_url = A.getDataBaseURL("pantoea[Organism]", "GenBank")

#if it returns a good link
if(database_url):
    ftpLink = database_url[27:]
    #print(ftpLink)

    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd(ftpLink)

    #print(database_url.rfind('/'))
    fileNameIndex = database_url.rfind('/') + 1

    #currently hard coding the _genomic.fna.gz part in. Possibly taken in from user later on
    print("Getting the file")
    filename = database_url[fileNameIndex:]
    filename = filename + "_genomic.fna.gz"


    #store the file downloaded into the parent directory of wherever the CASPERapp is
    local_filename = "../" + filename
    #print(local_filename)
    lf = open("../" + filename, "wb")

    print("Writing the file")
    #still expecting an archive file? Not sure
    ftp.retrbinary("RETR " + filename, lf.write, 8*1024)
    lf.close()

    ftp.close()

print('Done')
"""