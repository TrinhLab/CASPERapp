"""
This code was adapted by David Dooley from the code provided in Doench et al. 2016.

FULL CITATION:
Doench, J., Fusi, N., Sullender, M. et al.
Optimized sgRNA design to maximize activity and minimize
off-target effects of CRISPR-Cas9. Nat Biotechnol 34, 184–191 (2016).
https://doi.org/10.1038/nbt.3437

"""
#Calculates the Rule set 2 score for the given 30-mer
#Input: 1. 30mer sgRNA+context sequence, NNNN[sgRNA sequence]NGGNNN
#       2. Amino acid cut position, for full model prediction only
#       3. Percent peptide, for full model prediction only
#Output: Rule set 2 score

import pandas as pd
import csv, argparse, sys
import pickle
import model_comparison

def score(sequence):
    seq = sequence.upper()
    if len(seq)!=30: 
        print "Please enter a 30mer sequence."
        sys.exit(1)
    model_file_1 = '../saved_models/V3_model_nopos.pickle'
    try:
        with open(model_file, 'rb') as f:
            model = pickle.load(f)    
    except:
        raise Exception("Could not find model stored to file %s" % model_file)

    if seq[25:27] == 'GG':
        score = model_comparison.predict(seq, aa_cut, per_peptide, model=model)
        print 'Rule set 2 score: %.4f'% (score)
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'
