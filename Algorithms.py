"""Algorithms.py File.  This file contains the following classes: SeqTranslate.
    SeqTranslate Class. Used for interpreting base64 representations of the target locations as well as their sequences.
    To interpret these run the class instance at the bottom of the file with the desired base64 representation into the
    decompress_tuple function."""

import GlobalSettings
class SeqTranslate:

    def __init__(self):
        # Modification of MIME base64 coding so that +- can be used for strand direction
        self.base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/"

        self.endo_info = dict()
        self.endo_import()


    # used to convert numbers in base4 back to nucleotides
    def int2nt(self, num):
        if num == 0:
            return 'A'
        elif num == 1:
            return 'T'
        elif num == 2:
            return 'C'
        elif num == 3:
            return 'G'
        else:
            return 'N'

    def nt2int(self,nt):
        if nt == 'A':
            return 0
        elif nt == 'T':
            return 1
        elif nt == 'C':
            return 2
        elif nt == 'G':
            return 3
        else:
            return 0

    def compress(self, uncompressed, base):
        compseq = 0
        if type(uncompressed) == str:
            for i in range(len(uncompressed)):
                val = self.nt2int(uncompressed[i]) * pow(4, i)  # multiplying by power-4 converts to base10
                compseq += val
            uncompressed = compseq
        compreturn = str()
        while uncompressed >= base:
            rem = uncompressed%base
            uncompressed = int(uncompressed/base)
            compreturn = self.base_array_64[rem] + compreturn
        compreturn = self.base_array_64[uncompressed] + compreturn
        return compreturn

    def to_generic_compressed(self, seqobj):
        # Passed as a tuple from the repeats section of the .cspr file
        if type(seqobj) == list:
            gencomp = seqobj[0] + "." + seqobj[1][1:]
        else:
            split = seqobj.find("+")
            if split != -1:
                gencomp = seqobj[:split] + "." + seqobj[split+1:]
            else:
                split = seqobj.find("-")
                gencomp = seqobj[:split] + "." + seqobj[split+1:]
        return gencomp

    # Decompresses the base64 representation into base10.  If toseq is true it returns the sequence itself (nucleotides)
    def decompress64(self, base64seq, slength=0, toseq=False):
        base10seq = int()
        if isinstance(base64seq, str):
            for i in range(len(base64seq)):
                power = len(base64seq) - (i+1)
                index = self.base_array_64.find(base64seq[i])
                if index != -1:
                    base10seq += index*pow(64, power)
        else:
            base10seq = base64seq
        if toseq:
            seq = str()
            number = base10seq
            while number >= 4:
                rem = number % 4
                number = int(number/4)
                seq += self.int2nt(rem)
            seq += self.int2nt(number)
            for i in range(len(seq), slength):
                seq += 'A'
            return seq
        else:
            return base10seq

    def decompress_csf_tuple(self, locseq, bool=False, endo="spCas9"):
        # Lookup endonuclease sequence lengths for parsing
        if(bool == False):
            mytuple = locseq[:-1].split(",")
        else:
            mytuple = locseq.split(",")
            front_seq = mytuple[3]

        loc = self.decompress64(mytuple[0])
        seq = mytuple[1]
        scr = self.decompress64(mytuple[2])
        strand = seq.find("+")
        if strand != -1:
            dira = "+"
            sequence = seq[:strand]
            pam = seq[strand+1:]
        else:
            seq = seq.split("-")
            sequence = seq[0]
            pam = seq[1]
            dira = "-"
        if bool:
            seqlength = int(self.endo_info[endo][2]) - int(
                self.endo_info[endo][1])  # gets the tail sequence length for processing repeats
        else:
            seqlength = int(self.endo_info[endo][2])  # gets the total sequence length
        pamlength = len(self.endo_info[endo][0].split(",")[0])  # gets the length of the primary PAM
        #print(seqlength,pamlength)
        sequence = self.decompress64(sequence, seqlength, True)
        pam = self.decompress64(pam, pamlength, True)
        # The for loops fixes the problem of A's not being added to the end because they are removed on compression
        if(bool == True):
            sequence = sequence + front_seq
        return int(loc), str(sequence), pam, int(scr), dira, endo

    def endo_import(self):
        f = open(GlobalSettings.appdir + "CASPERinfo")
        while True:
            line = f.readline()
            if line.startswith("ENDONUCLEASES"):
                break
        while True:
            line = f.readline()
            if line.startswith("-"):
                break
            else:
                myinfo = line.split(";")
                self.endo_info[myinfo[0]] = myinfo[1:]  # first is PAM list, second is seed length, third is tot length




#S = SeqTranslate()
#print(S.decompress_csf_tuple("Dx,|S62qFEz+Qy,k", endo='asCas12'))
#print(S.decompress64("C86",False))
#print(S.compress(440159,64))

