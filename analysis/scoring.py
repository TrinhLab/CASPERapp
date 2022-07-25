import os
import csv
import numpy as np
from Bio import SeqIO
import azimuth.model_comparison as az

## Read in sequences to be scored and save to list:
seqs = []
spamReader = csv.reader(open('dnaJ_Doench2016.tsv', 'r'), delimiter='\t')
for i,row in enumerate(spamReader):
    if i == 0:
        pass
    else:
        seqs.append(str(row[1]).upper())


# Load in genome sequence
for record in SeqIO.parse("ecoli_BW25113.fasta", "fasta"):
    genome = record.seq.__str__()
    rev_genome = record.seq.reverse_complement().__str__()
it = 0

# Pull the full 30-nt sequence from the genome for each gRNA
full_seqs = []
for seq in seqs:
    tmp = genome.find(seq)
    if tmp != -1:
        full_seqs.append(genome[tmp-4:tmp+26])
    else:
        tmp = rev_genome.find(seq)
        full_seqs.append(rev_genome[tmp-4:tmp+26])

full_seqs = np.array(full_seqs)
tmp = az.predict(full_seqs)
print(tmp)

    # os.system(os.getcwd() + '/dist/doench_2016 --seq "' + str(seq) + '"')
