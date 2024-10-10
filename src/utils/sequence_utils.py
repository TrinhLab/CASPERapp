"""
Algorithms.py

This file contains utility classes and functions for sequence translation and manipulation.
"""

from PyQt6 import QtWidgets
import os
from typing import List, Tuple, Union

def get_table_headers(table: QtWidgets.QTableWidget) -> List[str]:
    """
    Get the headers from a QtTableWidget.
    
    Args:
        table (QtWidgets.QTableWidget): The table widget.
    
    Returns:
        List[str]: A list of header labels.
    """
    return [table.horizontalHeaderItem(c).text() if table.horizontalHeaderItem(c) else str(c+1) 
            for c in range(table.columnCount())]

class SeqTranslate:
    """
    Class for interpreting base64 representations of target locations and their sequences.
    """

    def __init__(self, casper_info_path: str):
        self.base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/"
        self.endo_info = {}
        self.endo_import(casper_info_path)

    @staticmethod
    def int2nt(num: int) -> str:
        """Convert integer to nucleotide."""
        return 'ATCGN'[min(num, 4)]

    @staticmethod
    def nt2int(nt: str) -> int:
        """Convert nucleotide to integer."""
        return 'ATCG'.index(nt) if nt in 'ATCG' else 0

    def compress(self, uncompressed: Union[str, int], base: int) -> str:
        """Compress a sequence or number to a base64-like representation."""
        if isinstance(uncompressed, str):
            uncompressed = sum(self.nt2int(nt) * (4 ** i) for i, nt in enumerate(uncompressed))
        
        compressed = []
        while uncompressed:
            uncompressed, rem = divmod(uncompressed, base)
            compressed.append(self.base_array_64[rem])
        return ''.join(reversed(compressed))

    def to_generic_compressed(self, seqobj: Union[List[str], str]) -> str:
        """Convert a sequence object to a generic compressed format."""
        if isinstance(seqobj, list):
            return f"{seqobj[0]}.{seqobj[1][1:]}"
        split = seqobj.find('+')
        if split == -1:
            split = seqobj.find('-')
        return f"{seqobj[:split]}.{seqobj[split+1:]}"

    def decompress64(self, base64seq: Union[str, int], slength: int = 0, toseq: bool = False) -> Union[int, str]:
        """Decompress a base64 representation to base10 or sequence."""
        if isinstance(base64seq, str):
            base10seq = sum(self.base_array_64.index(char) * (64 ** power) 
                            for power, char in enumerate(reversed(base64seq)) if char in self.base_array_64)
        else:
            base10seq = base64seq

        if toseq:
            seq = ''.join(self.int2nt(base10seq % 4) for _ in range(max(slength, len(str(base10seq)))))
            return seq[:slength] if slength else seq
        return base10seq

    def decompress_csf_tuple(self, locseq: str, bool: bool = False, endo: str = "spCas9") -> Tuple[int, str, str, int, str, str]:
        """Decompress a CSF tuple."""
        mytuple = locseq[:-1].split(",") if not bool else locseq.split(",")
        
        loc = self.decompress64(mytuple[0])
        seq = mytuple[1]
        scr = self.decompress64(mytuple[2])
        
        strand = seq.find("+")
        if strand != -1:
            dira, sequence, pam = "+", seq[:strand], seq[strand+1:]
        else:
            sequence, pam = seq.split("-")
            dira = "-"
        
        seqlength = int(self.endo_info[endo][2]) - int(self.endo_info[endo][1]) if bool else int(self.endo_info[endo][2])
        pamlength = len(self.endo_info[endo][0].split(",")[0])
        
        sequence = self.decompress64(sequence, seqlength, True)
        pam = self.decompress64(pam, pamlength, True)
        
        if bool:
            sequence += mytuple[3]
        
        return int(loc), str(sequence), pam, int(scr), dira, endo

    def endo_import(self, casper_info_path: str):
        """Import endonuclease information from CASPERinfo file."""
        with open(casper_info_path, 'r') as f:
            for line in f:
                if line.startswith("ENDONUCLEASES"):
                    break
            for line in f:
                if line.startswith("-"):
                    break
                myinfo = line.strip().split(";")
                self.endo_info[myinfo[0]] = myinfo[1:]

# Example usage:
# S = SeqTranslate(os.path.join(settings.get_app_dir(), "CASPERinfo"))
# print(S.decompress_csf_tuple("Dx,|S62qFEz+Qy,k", endo='asCas12'))
# print(S.decompress64("C86", False))
# print(S.compress(440159, 64))

