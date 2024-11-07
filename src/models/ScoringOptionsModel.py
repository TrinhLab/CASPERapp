import numpy as np
from Bio import SeqIO
import traceback
import warnings
import contextlib
import sys
import os

class ScoringOptionsModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.genome = ""
        self.rev_genome = ""
        self.fasta_path = ""
        
    def load_fasta(self, fasta_path, chromosome):
        """Load FASTA file for the specified chromosome"""
        try:
            self.fasta_path = fasta_path
            
            # Extract chromosome number from ID (e.g., "NZ_CP032679.1" -> "1")
            chrom_num = int(chromosome.split('.')[-1]) if '.' in chromosome else 1
            self.logger.debug(f"Looking for chromosome number: {chrom_num}")
            
            # Load only the required chromosome
            for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
                if i + 1 == chrom_num:  # 1-based indexing
                    self.genome = str(record.seq).upper()
                    self.rev_genome = str(record.seq.reverse_complement()).upper()
                    self.logger.debug(f"Loaded chromosome {chrom_num} sequence of length {len(self.genome)}")
                    break
                    
            return True
            
        except Exception as e:
            self.logger.error(f"Error loading FASTA file: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return False

    def score_sequences(self, targets, algorithm="Azimuth 2.0"):
        """Score the target sequences using the specified algorithm"""
        try:
            guide_list = []
            reject_list = []
            full_seqs = []
            
            # Process each target
            for i, target in enumerate(targets):
                strand = target['strand']
                sequence = target['sequence'] + target['pam']
                guide_list.append(sequence)
                
                # Search for sequence in genome
                if strand == "+":
                    pos = self.genome.find(sequence)
                    if pos != -1:
                        full_seqs.append(self.genome[pos-4:pos+26])
                    else:
                        reject_list.append(i)
                else:
                    pos = self.rev_genome.find(sequence)
                    if pos != -1:
                        full_seqs.append(self.rev_genome[pos-4:pos+26])
                    else:
                        reject_list.append(i)
            
            # Score sequences if any were found
            if full_seqs:
                full_seqs = np.array(full_seqs)
                
                # Score using selected algorithm
                if algorithm == "Azimuth 2.0":
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        # Add utils directory to Python path
                        utils_path = os.path.join(self.global_settings.get_src_dir_path(), 'utils')
                        if utils_path not in sys.path:
                            sys.path.append(utils_path)
                            
                        from azimuth import model_comparison as az
                        scores = az.predict(full_seqs) * 100
                else:
                    raise ValueError(f"Unknown algorithm: {algorithm}")
                
                # Insert -1 scores for rejected sequences
                for i in reject_list:
                    scores = np.insert(scores, i, -1)
                    
                return scores, reject_list, guide_list
                
            return None, reject_list, guide_list
            
        except Exception as e:
            self.logger.error(f"Error scoring sequences: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return None, [], []
