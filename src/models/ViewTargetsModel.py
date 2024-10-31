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
            self.global_settings.logger.debug(f"ViewTargetsModel clearing caches for new annotation file: {new_annotation_file}")
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
            self.global_settings.logger.error(f"Error in _on_annotation_file_changed: {str(e)}")

    def load_targets(self, selected_targets, organism, endonuclease):
        """Fast target loading with minimal file operations"""
        start_time = time.time()
        
        try:
            self.global_settings.logger.debug(f"Starting load_targets with {len(selected_targets)} targets")
            
            # Store organism and endonuclease for potential reloading
            self.organism = organism
            self.endonuclease = endonuclease
            
            # Get CSPR parser from cache or create new one
            parser_start = time.time()
            cspr_key = f"{organism}_{endonuclease}"
            if cspr_key in self._parser_cache:
                self.cspr_parser = self._parser_cache[cspr_key]
            else:
                org_files = self.get_organism_to_files()
                if organism not in org_files or endonuclease not in org_files[organism]:
                    self.global_settings.logger.error(f"No CSPR file found for {organism} and {endonuclease}")
                    return

                cspr_file = org_files[organism][endonuclease][0]
                cspr_path = os.path.join(self.global_settings.get_db_path(), cspr_file)
                self.cspr_parser = CSPRparser(cspr_path, self.global_settings.get_casper_info_path())
                self._parser_cache[cspr_key] = self.cspr_parser
            parser_time = time.time() - parser_start

            # Initialize targets and genes
            self.targets = []
            self.available_genes = set()
            
            # Set up annotation parser if needed
            if self.annotation_parser is None:
                annotation_start = time.time()
                self.annotation_parser = AnnotationParser(self.global_settings)
                annotation_files = self.get_annotation_files()
                if annotation_files:
                    self.annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_files[0])
                    self.annotation_parser.set_annotation_file(self.annotation_path)
                annotation_time = time.time() - annotation_start
            else:
                annotation_time = 0

            # Process targets in batches by chromosome
            processing_start = time.time()
            
            # Group targets by chromosome and prepare batch reading
            batch_targets = defaultdict(list)
            for target in selected_targets:
                chrom = target['chromosome']
                start, end = map(int, target['location'].split('-'))
                batch_targets[chrom].append({
                    'feature_name': target['feature_name'],
                    'start': start,
                    'end': end
                })
                self.available_genes.add(target['feature_name'])

            # Batch process targets for each chromosome
            target_count = 0
            for chrom, targets in batch_targets.items():
                self.chromosome = chrom
                
                # Sort targets by start position for more efficient reading
                targets.sort(key=lambda x: x['start'])
                
                # Read targets in a single batch per chromosome
                batch_results = self.cspr_parser.read_targets_batch(
                    chromosome=chrom,
                    targets=targets,
                    endonuclease=endonuclease
                )
                
                if batch_results:
                    self.targets.extend(batch_results)
                    target_count += len(batch_results)

            processing_time = time.time() - processing_start

            # Convert genes to sorted list
            self.available_genes = sorted(list(self.available_genes))
            
            total_time = time.time() - start_time
            self.global_settings.logger.debug(f"Total load_targets execution time: {total_time:.2f} seconds")
            self.global_settings.logger.debug(f"Found {target_count} total CSPR targets")

        except Exception as e:
            self.global_settings.logger.error(f"Error in load_targets: {str(e)}\n{traceback.format_exc()}")
            raise

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

    def get_gene_data(self, gene_name):
        """Get gene data with caching"""
        try:
            if not gene_name:
                self.global_settings.logger.error("No gene name provided")
                return None
                
            # Check model cache first
            if gene_name in self._gene_data_cache:
                return self._gene_data_cache[gene_name]
            
            # Make sure annotation parser is initialized
            if self.annotation_parser is None:
                self._initialize_annotation_parser()
                
            # Get gene data from parser
            gene_data = self.annotation_parser.get_gene_data(gene_name)
            if gene_data:
                self._gene_data_cache[gene_name] = gene_data
                
            return gene_data
            
        except Exception as e:
            self.global_settings.logger.error(f"Error getting gene data: {str(e)}")
            return None

    def get_targets(self):
        return self.targets

    def highlight_targets_in_gene_viewer(self, selected_targets):
        """Highlight selected targets in gene viewer"""
        try:
            self.global_settings.logger.debug("Starting highlight_targets_in_gene_viewer")
            sequence = self.extended_sequence
            if not sequence:
                self.global_settings.logger.error("No extended sequence available")
                return sequence

            self.global_settings.logger.debug(f"Extended sequence length: {len(sequence)}")
            
            # Sort targets by position for efficient highlighting
            highlights = []
            for target in selected_targets:
                self.global_settings.logger.debug(f"Processing target: {target}")
                sequence_to_find = target['sequence']
                strand = target['strand']
                
                # For negative strand, we need to use reverse complement
                if strand == '-':
                    sequence_to_find = str(Seq(sequence_to_find).reverse_complement())
                    self.global_settings.logger.debug(f"Reverse complemented sequence: {sequence_to_find}")
                
                # Search for the sequence in the gene viewer text
                sequence_upper = sequence.upper()
                target_upper = sequence_to_find.upper()
                
                self.global_settings.logger.debug(f"Searching for sequence: {target_upper}")
                
                # Find all occurrences
                pos = sequence_upper.find(target_upper)
                if pos != -1:
                    self.global_settings.logger.debug(f"Found sequence at position: {pos}")
                    color = 'red' if strand == '-' else 'green'
                    highlights.append((pos, len(sequence_to_find), color))
                else:
                    self.global_settings.logger.warning(f"Sequence not found: {target_upper}")

            if not highlights:
                self.global_settings.logger.error("No sequences could be highlighted")
                return sequence

            self.global_settings.logger.debug(f"Found {len(highlights)} sequences to highlight")

            # Build highlighted sequence
            result = []
            last_pos = 0
            for pos, length, color in highlights:
                result.append(sequence[last_pos:pos])
                result.append(f"<span style='background-color: {color};'>")
                result.append(sequence[pos:pos+length])
                result.append("</span>")
                last_pos = pos + length
            
            result.append(sequence[last_pos:])
            final_sequence = ''.join(result)
            
            self.global_settings.logger.debug(f"Final highlighted sequence length: {len(final_sequence)}")
            return final_sequence

        except Exception as e:
            self.global_settings.logger.error(f"Error highlighting targets: {str(e)}\n{traceback.format_exc()}")
            return sequence

    def get_available_genes(self):
        """Get list of available genes from the loaded targets"""
        try:
            # Return the available genes list that was populated during load_targets
            if hasattr(self, 'available_genes'):
                return self.available_genes
            
            # If not already populated, get unique genes from targets
            genes = set()
            for target in self.targets:
                if 'feature_name' in target:
                    genes.add(target['feature_name'])
            
            # Store for future use
            self.available_genes = sorted(list(genes))
            return self.available_genes
            
        except Exception as e:
            self.global_settings.logger.error(f"Error getting available genes: {str(e)}")
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
