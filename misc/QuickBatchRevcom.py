# to run quick things such as getting revcoms of sequences

def revcom(sequence):
    revseq = ""
    change = {'A': 'T',
              'T': 'A',
              'G': 'C',
              'C': 'G'}
    for nt in sequence:
        rnt = change[nt]
        revseq = rnt + revseq
    return revseq


f = open("/Users/brianmendoza/Desktop/primer_sequences.txt")
for line in f:
    myseq = line[:-1].upper()
    print(revcom(myseq))
