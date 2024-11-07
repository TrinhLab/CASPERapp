from models.CSPRparser import CSPRparser
from models.HomeWindowModel import HomeWindowModel
from models.AnnotationParser import AnnotationParser
import os
from Bio import SeqIO
from Bio.Seq import Seq
from functools import lru_cache
import threading
from collections import defaultdict
import re
import traceback
import time
import logging
from multiprocessing import Pool
from functools import partial
from multiprocessing import cpu_count

class ViewTargetsModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.targets = []
        self.cspr_parser = None
        self.annotation_parser = None  # Will be initialized when needed
        self.gene_sequence = ""
        self.highlighted_sequence = ""
        self.gene_info = {}
        self.available_genes = []
        self.filter_options = {}
        self.scoring_options = {}
        self.annotation_path = ""
        self.current_gene_start = 0
        self.current_gene_end = 0
        self.extended_sequence = ""
        self.chromosome = ""
        
        # Cache structures
        self._gene_data_cache = {}
        self._sequence_cache = {}
        self._parser_cache = {}
        self._chromosome_seqs = {}
        self._cached_targets = {}  # Add cache for targets

        # Connect to annotation file changes
        self.global_settings.annotation_file_changed.connect(self._on_annotation_file_changed)
        
        # Initialize annotation path
        self.annotation_path = os.path.join(
            self.global_settings.get_db_path(),
            'GBFF',
            self.global_settings.get_current_annotation_file()
        )
        self.logger.debug(f"Initialized annotation path: {self.annotation_path}")

    def cleanup(self):
        """Cleanup method to be called when the view is closed"""
        try:
            # Disconnect from annotation file changes
            if hasattr(self, '_annotation_signal'):
                self.global_settings.annotation_file_changed.disconnect(self._on_annotation_file_changed)
                self.global_settings.logger.debug("ViewTargetsModel disconnected from annotation file changes")
            
            # Clear caches
            self._gene_data_cache.clear()
            self._sequence_cache.clear()
            self._parser_cache.clear()
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in ViewTargetsModel cleanup: {str(e)}")

    def _on_annotation_file_changed(self, new_annotation_file):
        """Clear all caches when annotation file changes"""
        try:
            self.logger.debug(f"ViewTargetsModel clearing caches for new annotation file: {new_annotation_file}")
            self._gene_data_cache.clear()
            self._sequence_cache.clear()
            self._parser_cache.clear()
            
            # Update annotation path and parser
            self.annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', new_annotation_file)
            self.annotation_parser = AnnotationParser(self.global_settings)
            self.annotation_parser.set_annotation_file(self.annotation_path)
            
            # Clear other stored data
            self.gene_sequence = ""
            self.highlighted_sequence = ""
            self.gene_info = {}
            self.available_genes = []
            self._chromosome_seqs = {}
            
        except Exception as e:
            self.logger.error(f"Error in _on_annotation_file_changed: {str(e)}")

    def load_targets(self, selected_targets, organism, endonuclease):
        """Fast target loading with minimal file operations"""
        total_start = time.time()
        
        try:
            self.logger.debug(f"Starting load_targets with {len(selected_targets)} targets")
            
            # Store organism and endonuclease for potential reloading
            self.organism = organism
            self.endonuclease = endonuclease

            # Get CSPR parser from cache or create new one
            parser_start = time.time()
            cspr_key = f"{organism}_{endonuclease}"
            if cspr_key in self._parser_cache:
                self.cspr_parser = self._parser_cache[cspr_key]
                self.logger.debug("Using cached CSPR parser")
            else:
                org_files = self.get_organism_to_files()
                if organism not in org_files or endonuclease not in org_files[organism]:
                    self.logger.error(f"No CSPR file found for {organism} and {endonuclease}")
                    return

                cspr_file = org_files[organism][endonuclease][0]
                cspr_path = os.path.join(self.global_settings.get_db_path(), cspr_file)
                self.cspr_parser = CSPRparser(cspr_path, self.global_settings.get_casper_info_path())
                self._parser_cache[cspr_key] = self.cspr_parser
                self.logger.debug("Created new CSPR parser")
            parser_time = time.time() - parser_start
            self.logger.debug(f"CSPR parser initialization time: {parser_time:.2f} seconds")

            # Initialize targets and genes
            init_start = time.time()
            self.targets = []
            self.available_genes = set()
            init_time = time.time() - init_start
            self.logger.debug(f"Initialization time: {init_time:.2f} seconds")
            
            # Group targets by chromosome
            group_start = time.time()
            batch_targets = defaultdict(list)
            for target in selected_targets:
                chrom = target['chromosome']
                start, end = map(int, target['location'].split('-'))
                batch_targets[chrom].append({
                    'feature_name': target['feature_name'],
                    'feature_id': target['feature_id'],  # Include feature_id (locus_tag)
                    'start': start,
                    'end': end
                })
                # Store both feature_id and feature_name
                self.available_genes.add((target['feature_id'], target['feature_name']))
            group_time = time.time() - group_start
            self.logger.debug(f"Target grouping time: {group_time:.2f} seconds")

            # Process targets by chromosome
            process_start = time.time()
            target_count = 0
            for chrom, targets in batch_targets.items():
                batch_start = time.time()
                results = self.cspr_parser.read_targets_batch(chrom, targets, endonuclease)
                # Add feature_id to each result
                for result in results:
                    # Find matching target to get feature_id
                    for target in targets:
                        if (target['start'] <= result['position'] <= target['end'] and 
                            target['feature_name'] == result['feature_name']):
                            result['feature_id'] = target['feature_id']
                            break
                self.targets.extend(results)
                target_count += len(results)
                batch_time = time.time() - batch_start
                self.logger.debug(f"Chromosome {chrom} processing time: {batch_time:.2f} seconds")
            process_time = time.time() - process_start
            self.logger.debug(f"Total target processing time: {process_time:.2f} seconds")

            total_time = time.time() - total_start
            self.logger.debug(f"Total load_targets execution time: {total_time:.2f} seconds")
            self.logger.debug(f"Found {target_count} total CSPR targets")
            
        except Exception as e:
            self.logger.error(f"Error in load_targets: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def _get_chromosome_sequence(self, chromosome):
        """Get chromosome sequence on demand"""
        if not hasattr(self, '_chromosome_seqs'):
            self._chromosome_seqs = {}
        
        if chromosome not in self._chromosome_seqs:
            for record in SeqIO.parse(self.annotation_path, "genbank"):
                if record.id == chromosome:
                    self._chromosome_seqs[chromosome] = str(record.seq)
                    break
        
        return self._chromosome_seqs.get(chromosome)

    def _initialize_annotation_parser(self):
        """Initialize annotation parser if not already initialized"""
        if self.annotation_parser is None:
            self.annotation_parser = AnnotationParser(self.global_settings)
            if self.annotation_path:
                self.annotation_parser.set_annotation_file(self.annotation_path)

    def get_gene_data(self, locus_tag):
        """Get gene data with proper error handling"""
        try:
            if not locus_tag:
                self.logger.debug("No locus tag provided")
                return None
                
            # Check model cache first
            if locus_tag in self._gene_data_cache:
                return self._gene_data_cache[locus_tag]
            
            # Initialize annotation parser if not already done
            if not hasattr(self, 'annotation_parser') or self.annotation_parser is None:
                self.annotation_parser = AnnotationParser(self.global_settings)
                annotation_file = self.global_settings.get_current_annotation_file()
                annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
                self.annotation_parser.set_annotation_file(annotation_path)
                self.logger.debug(f"Initialized annotation parser with file: {annotation_path}")
            
            # Get gene data from parser with proper string conversion
            gene_data = None
            if isinstance(locus_tag, (str, int)):
                locus_tag_str = str(locus_tag).strip()
                self.logger.debug(f"Searching for locus tag: {locus_tag_str}")
                # Look up by locus tag directly
                gene_data = self.annotation_parser.get_gene_data(locus_tag_str.lower())
                
            if gene_data:
                self._gene_data_cache[locus_tag] = gene_data
                self.logger.debug(f"Found gene data: {gene_data.keys()}")
            else:
                self.logger.debug(f"No gene data found for locus tag: {locus_tag}")
                
            return gene_data
            
        except Exception as e:
            self.logger.error(f"Error getting gene data: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return None

    def get_targets(self):
        """Return all targets with their feature IDs"""
        return self.targets

    def get_available_genes(self):
        """Get list of available genes with format 'feature_id: feature_name'"""
        try:
            if hasattr(self, 'available_genes'):
                # Format as "feature_id: feature_name"
                return [f"{feature_id}: {feature_name}" 
                       for feature_id, feature_name in sorted(self.available_genes)]
            return []
        except Exception as e:
            self.logger.error(f"Error getting available genes: {str(e)}")
            return []

    # ... (other methods remain unchanged)

    def _process_target(self, target):
        """Process a single target - moved to separate method for parallel processing"""
        try:
            # Your existing target processing logic here
            # Make sure to handle any shared resources thread-safely
            pass
        except Exception as e:
            logging.error(f"Error processing target: {e}")
            return None

    def get_gene_sequence(self, locus_tag):
        """Get gene sequence with optimized caching and minimal I/O"""
        try:
            # Check sequence cache first
            cache_key = f"{locus_tag}_sequence"
            if cache_key in self._sequence_cache:
                self.logger.debug(f"Cache hit for gene sequence: {locus_tag}")
                return self._sequence_cache[cache_key]
                
            # Get gene data which includes location information
            print(f"Getting gene data for locus tag: {locus_tag}")
            gene_data = self.get_gene_data(locus_tag)
            if not gene_data or 'info' not in gene_data:
                self.logger.warning(f"No gene data found for locus tag: {locus_tag}")
                return None
                
            # Parse location string (format: "start:end(strand)")
            location = gene_data['info']['location']
            if ':' not in location:
                self.logger.warning(f"Invalid location format: {location}")
                return None
                
            # Extract start and end positions
            start = int(location.split(':')[0])
            end = int(location.split(':')[1].split('(')[0])
            chromosome = gene_data['info']['chromosome']
            
            # Get sequence from gene_data directly if available
            if 'sequence' in gene_data:
                sequence = gene_data['sequence']
                
                # Add padding (30 bases on each side)
                padding = 30
                seq_start = max(0, start - padding)
                seq_end = min(len(sequence), end + padding)
                
                # Get sequence with padding
                five_prime_pad = sequence[seq_start:start].lower() if seq_start < start else ""
                main_seq = sequence[start:end].upper()
                three_prime_pad = sequence[end:seq_end].lower() if end < seq_end else ""
                
                full_sequence = five_prime_pad + main_seq + three_prime_pad
                
                # Cache the result
                result = {
                    'sequence': full_sequence,
                    'chrom_length': len(sequence),
                    'start': start,
                    'end': end,
                    'padded_start': seq_start,
                    'padded_end': seq_end
                }
                self._sequence_cache[cache_key] = result
                
                self.logger.debug(f"Retrieved and cached sequence for locus tag {locus_tag} ({len(full_sequence)} bp)")
                return result
                
        except Exception as e:
            self.logger.error(f"Error getting gene sequence: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return None

    def get_scoring_options(self):
        """Get current scoring options"""
        try:
            if not hasattr(self, 'scoring_options'):
                self.scoring_options = {
                    'algorithm': 'Azimuth 2.0',
                    'fasta_file': '',
                    'min_score': 0,
                    'max_score': 100
                }
            return self.scoring_options
            
        except Exception as e:
            self.logger.error(f"Error getting scoring options: {str(e)}")
            return {}

    def set_scoring_options(self, options):
        """Set scoring options"""
        try:
            self.scoring_options = options
            self.logger.debug(f"Updated scoring options: {options}")
            
        except Exception as e:
            self.logger.error(f"Error setting scoring options: {str(e)}")