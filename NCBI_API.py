"""Parsing GBFF, GFF, and tabular genome annotations.  GFF and tabular do not have the nucleotide sequence while GBFF
has everything"""

from Bio import SeqIO, Entrez
from bs4 import BeautifulSoup
import requests
from ftplib import FTP
import time
import GlobalSettings
import gzip

Entrez.email = "bmendoz1@vols.utk.edu"


class GBFF_Parse:

    def __init__(self, filename, NCBI=False):
        self.filename = filename

    def convert_to_fasta(self):
        count = SeqIO.convert(self.filename, "genbank", self.filename[:-5]+"fasta", "fasta")


class Assembly:

    def __init__(self):
        self.gca_rectList = list()
        self.database_url_list = list()
        self.orgName_dict = dict()
        #Moved the code from this to its own function so that I can return the database URL

    #this function gets a list of GCA's and a list of database URLs for downloading
    def getDataBaseURL(self, organism, database):
        progValue = 15
        handle = Entrez.esearch(db="assembly", retmax=20, term=organism)
        record = Entrez.read(handle)
        self.database = database
        myidlist = list()

        for id in record["IdList"]:
            myidlist.append(id)
            progValue += 1
            GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(progValue)
            #print(id)
        handle.close()

        gca_rectList = list()
        for ret in myidlist:
            handle = Entrez.esummary(db="assembly", id=ret)
            gca_rec = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            #print(gca_rec)
            self.gca_rectList.append(gca_rec)
            handle.close()
            time.sleep(0.5)
            progValue += 1
            GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(progValue)

        for i in range(len(self.gca_rectList)):
            progValue += 1
            GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(progValue)

            url = "https://www.ncbi.nlm.nih.gov/assembly/" + self.gca_rectList[i] + "/"
            source = requests.get(url)
            plain_text = source.text
            soup = BeautifulSoup(plain_text, "html.parser")

            refseq_link = str(soup.find('a', text="Download the RefSeq assembly"))
            genbank_link = str(soup.find('a', text="Download the GenBank assembly"))
            orgName = soup.find('dd')
            refseq_link = refseq_link[refseq_link.find("=") + 2: refseq_link.find(">") - 1]
            genbank_link = genbank_link[genbank_link.find("=") + 2: genbank_link.find(">") - 1]
            # Checkpoint printing
            #print("RefSeq link: ", refseq_link)
            #print("Genbank link: ", genbank_link)
            # then go to that webpage then inspect and download with url service the gbff file
            if self.database == "GenBank":
                # check and see if GCF is in the the GCA_ID, if so, swap it with GCA
                # not sure if doing this is correct or not, but this way each description actually goes to a download link
                if "GCF" in self.gca_rectList[i]:
                    self.gca_rectList[i] = self.gca_rectList[i].replace("GCF", "GCA")
                database_url = genbank_link
            else:
                database_url = refseq_link
            # data_source = requests.get(database_url)
            # soup = BeautifulSoup(data_source.text, "html.parser")
            for link in soup.find_all('a', {'class': 'icon file'}):
                print("Link in for loop ", link)

            #check the links and catch errors
            if(database == "RefSeq" and len(refseq_link) < 5):
                print("Error: No RefSeq file to download!")
                #return ("Error: No RefSeq file to download!", list())
            elif(database == "GenBank" and len(genbank_link) < 5):
                print("Error: No GenBank file to download!")
                #return ("Error: No GenBank file to download!", list())
            elif(len(genbank_link) < 5 and len(refseq_link) < 5):
                print("Error: No RefSeq or GenBank files to download")
                #return ("Error: No RefSeq or GenBank files to download", list())
            else:
                #print(database_url)
                #return database_url
                self.database_url_list.append(database_url)
                self.orgName_dict[orgName.string] = self.gca_rectList[i]

        if len(self.database_url_list) <= 0:
            print("Error: no matches found. Please try again")
            return(self.database_url_list, self.orgName_dict)
        else:
            return (self.database_url_list, self.orgName_dict)

    # this function downloads the .gz file from the given link, and stores it in the CSPR_DB folder from Global Settings
    # NOTE: this function should only be used to download .gz files, nothing else
    # NOTE2: this function also only downloads the fna files
    def download_compressed_file(self, database_link):
        print("Downloading: ", database_link)
        ftpLink = database_link[27:]
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd(ftpLink)

        fileNameIndex = database_link.rfind('/') + 1
        filename = database_link[fileNameIndex:]
        filename = filename + "_genomic.fna.gz"

        local_filename = GlobalSettings.CSPR_DB + "/" +  filename #where the file is stored on the computer
        lf = open(local_filename, "wb")

        ftp.retrbinary("RETR " + filename, lf.write, 8 * 1024)
        lf.close()
        ftp.close()
        print("\tDone downloading")
        return local_filename

    # this function downloads the .gz file from the given link, and stores it in CSPR_DB folder from Global Settings
    # Note: this function should only be used to download .gz files and nothing else
    # Note2: this function will only be able to download feature_tables, gff, gbff files only.
    def download_compressed_annotation_file(self, database_link, fileType):
        # get to the folder on the ftp link
        print("Downloading: ", database_link)
        ftpLink = database_link[27:]
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd(ftpLink)

        # get the right file
        fileNameIndex = database_link.rfind('/') + 1
        filename = database_link[fileNameIndex:]
        filename = filename + fileType

        # open the right file in the correct location
        local_filename = GlobalSettings.CSPR_DB + "/" + filename
        lf = open(local_filename, "wb")

        # try to write the file, if it is not possible, throw an error
        # write it, and close everything up. Then return the file's location
        try:
            ftp.retrbinary("RETR " + filename, lf.write, 8 * 1024)
        except:
            lf.close()
            ftp.close()
            return
        lf.close()
        ftp.close()
        print("\tDone Downloading")
        return local_filename

    # this function decompresses .gz files into .fna files
    def decompress_file(self, file_path):
        # get the file name
        storeFileName = file_path[:-3]
        print("Decompressing: ", file_path)

        # open the file streams
        readStream = gzip.open(file_path, 'rb')
        writeStream = open(storeFileName, 'w')

        # read the compressed data and close that file
        zipFileContent = readStream.read()
        readStream.close()

        # write that data, and close that file
        tempString = str(zipFileContent)[2:-1]
        tempString = tempString.replace("\\t", "\t")
        tempList = tempString.split("\\n")
        for i in range(len(tempList)):
            writeStream.write(str(tempList[i]) + "\n")
        writeStream.close()

        print("\tDone decompressing")
        return storeFileName

    # this function decompressed an annotation file based on the type of annotation file given
    # stores the file in the CSPR_DB folder, and keeps the same name, just minus the .gz
    # may not need to check for the type of file, not sure. Still need to test
    def decompress_annotation_file(self, file_path, file_type):
        # get the compressed data
        storeFileName = file_path[:-3]
        print("Decompressing: ", file_path)
        readStream = gzip.open(file_path, 'rb')
        writeStream = open(storeFileName, 'w')
        zipFileContent = readStream.read()
        readStream.close()

        # check the file type, and decompress accordingly
        if file_type == "_genomic.gbff.gz":
            print("Decompress a gbff file here")
        elif file_type == "_genomic.gff.gz":
            tempString = str(zipFileContent)[2:-1]
            tempString = tempString.replace("\\t", "\t")
            tempList = tempString.split("\\n")
            for i in range(len(tempList)):
                writeStream.write(str(tempList[i] + "\n"))
        elif file_type == "_feature_table.txt.gz":
            tempString = str(zipFileContent)[2:-1]
            tempString = tempString.replace("\\t", "\t")
            tempList = tempString.split("\\n")
            for i in range(len(tempList)):
                writeStream.write(str(tempList[i] + "\n"))

        writeStream.close()
        print("\tDone Decompressing")
        return(storeFileName)


    # this function is similar to the getDatabaseURl function, however if the intial search yields 0 results,
    # it splits the organism and searches the rest of them as well
    def get_annotation_file(self, organism, database):
        # go through the handle and entrez for the first time
        searchOrganism = organism + "[Organism]"
        print("searching ncbi for: ", organism)
        handle = Entrez.esearch(db="assembly", retmax=20, term=searchOrganism)
        record = Entrez.read(handle)
        self.database = database
        self.database_url_list.clear()
        self.orgName_dict.clear()
        myidlist = list()

        for id in record["IdList"]:
            myidlist.append(id)
            # print(id)

        # if there's nothing in the list, go through and search the rest of them, otherwise skip this
        if len(myidlist) == 0:
            splitOrganism = organism.split(" ")
            for i in range(len(splitOrganism)):
                print("searching ncbi for: ", splitOrganism[i])
                searchOrganism = splitOrganism[i] + "[Organism]"
                handle = Entrez.esearch(db="assembly", retmax=20, term=searchOrganism)
                record = Entrez.read(handle)

                for id in record["IdList"]:
                    myidlist.append(id)

        # if the len of myIdList is still 0, then return out
        if len(myidlist) == 0:
            return (self.database_url_list, self.orgName_dict)

        gca_rectList = list()
        for ret in myidlist:
            handle = Entrez.esummary(db="assembly", id=ret)
            gca_rec = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            #print(gca_rec)
            self.gca_rectList.append(gca_rec)
            handle.close()
            time.sleep(0.5)

        for i in range(len(self.gca_rectList)):
            url = "https://www.ncbi.nlm.nih.gov/assembly/" + self.gca_rectList[i] + "/"
            source = requests.get(url)
            plain_text = source.text
            soup = BeautifulSoup(plain_text, "html.parser")

            refseq_link = str(soup.find('a', text="Download the RefSeq assembly"))
            genbank_link = str(soup.find('a', text="Download the GenBank assembly"))
            orgName = soup.find('dd')
            refseq_link = refseq_link[refseq_link.find("=") + 2: refseq_link.find(">") - 1]
            genbank_link = genbank_link[genbank_link.find("=") + 2: genbank_link.find(">") - 1]
            # Checkpoint printing
            # print("RefSeq link: ", refseq_link)
            # print("Genbank link: ", genbank_link)
            if self.database == "GenBank":
                # check and see if GCF is in the the GCA_ID, if so, swap it with GCA
                # not sure if doing this is correct or not, but this way each description actually goes to a download link
                if "GCF" in self.gca_rectList[i]:
                    self.gca_rectList[i] = self.gca_rectList[i].replace("GCF", "GCA")
                database_url = genbank_link
            else:
                database_url = refseq_link
            # data_source = requests.get(database_url)
            # soup = BeautifulSoup(data_source.text, "html.parser")
            for link in soup.find_all('a', {'class': 'icon file'}):
                print("Link in for loop ", link)

            # check the links and catch errors
            if (database == "RefSeq" and len(refseq_link) < 5):
                print("Error: No RefSeq file to download!")
                # return ("Error: No RefSeq file to download!", list())
            elif (database == "GenBank" and len(genbank_link) < 5):
                print("Error: No GenBank file to download!")
                # return ("Error: No GenBank file to download!", list())
            elif (len(genbank_link) < 5 and len(refseq_link) < 5):
                print("Error: No RefSeq or GenBank files to download")
                # return ("Error: No RefSeq or GenBank files to download", list())
            else:
                # print(database_url)
                # return database_url
                self.database_url_list.append(database_url)
                self.orgName_dict[orgName.string] = self.gca_rectList[i]

        if len(self.database_url_list) <= 0:
            print("Error: no matches found. Please try again")
            return (self.database_url_list, self.orgName_dict)
        else:
            return (self.database_url_list, self.orgName_dict)


# testing below
#myNCBI = Assembly()
#myNCBI.get_annotation_file("Bacillus Subtilis[Organism]", "GenBank")
