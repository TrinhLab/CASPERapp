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
import logging

class ViewTargetsModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.guides = []
        self.cspr_parser = None
        self.annotation_parser = None 
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
        
        self._gene_data_cache = {}
        self._sequence_cache = {}
        self._parser_cache = {}
        self._chromosome_seqs = {}
        self._cached_guides = {}

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

    def load_guides(self, selected_targets, organism, endonuclease):
        """Load guides with proper error handling"""
        try:
            self.logger.debug(f"Starting load_guides with {len(selected_targets)} targets")
            
            self.organism = organism
            self.endonuclease = endonuclease

            # Get CSPR parser from cache or create new one
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

            # Initialize guides and genes
            self.guides = []
            self.available_genes = set()
            
            # Use a set to track unique guide positions
            seen_guides = set()
            
            # Create chromosome mapping by counting carets
            chrom_mapping = {}
            chrom_count = 0
            annotation_file = self.global_settings.get_current_annotation_file()
            annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
            
            for record in SeqIO.parse(annotation_path, "genbank"):
                chrom_count += 1
                chrom_mapping[record.id] = str(chrom_count)
            
            batch_guides = defaultdict(list)
            for target in selected_targets:
                # Get chromosome number from mapping if full_chromosome is available
                if 'full_chromosome' in target:
                    chrom = chrom_mapping.get(target['full_chromosome'], target['chromosome'])
                else:
                    chrom = target['chromosome']
                
                start, end = map(int, target['location'].split('-'))
                
                # Create a unique identifier for this position range
                position_key = f"{chrom}:{start}-{end}"
                
                # Only add if we haven't seen this position range before
                if position_key not in seen_guides:
                    seen_guides.add(position_key)
                    batch_guides[chrom].append({
                        'feature_name': target['feature_name'],
                        'feature_id': target['feature_id'],
                        'start': start,
                        'end': end
                    })
                    self.available_genes.add((target['feature_id'], target['feature_name']))

            # Process guides by chromosome
            unique_guides = {}  # Use dict to track unique guides by sequence
            for chrom, guides in batch_guides.items():
                results = self.cspr_parser.read_targets_batch(chrom, guides, endonuclease)
                
                # Add feature_id to each result and deduplicate
                for result in results:
                    # Create a unique key for each guide
                    guide_key = (result['sequence'], result['position'], result['strand'])
                    
                    if guide_key not in unique_guides:
                        # Find matching guide to get feature_id
                        for target in selected_targets:
                            if (target['start'] <= result['position'] <= target['end'] and 
                                target['feature_name'] == result['feature_name']):
                                result['feature_id'] = target['feature_id']
                                unique_guides[guide_key] = result
                                break

            # Convert unique guides back to list
            self.guides = list(unique_guides.values())
            
            self.logger.debug(f"Found {len(self.guides)} unique guides")
            
        except Exception as e:
            self.logger.error(f"Error in load_guides: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def _get_chromosome_sequence(self, chromosome):
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

    def get_guides(self):
        """Return all guides with their feature IDs"""
        return self.guides

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

    def _process_guide(self, guide):
        """Process a single guide - moved to separate method for parallel processing"""
        try:
            # Your existing guide processing logic here
            # Make sure to handle any shared resources thread-safely
            pass
        except Exception as e:
            logging.error(f"Error processing guide: {e}")
            return None

    def get_gene_sequence(self, identifier):
        """Get gene sequence with optimized caching and minimal I/O"""
        try:
            print(f"Getting gene sequence for identifier: {identifier}")
            # Check sequence cache first
            cache_key = f"{identifier}_sequence"
            if cache_key in self._sequence_cache:
                self.logger.debug(f"Cache hit for sequence: {identifier}")
                return self._sequence_cache[cache_key]
            
            # Check if this is a position-based search
            if "chrom" in identifier and "start:" in identifier:
                try:
                    # Parse position from the text (format: "chrom X, start: Y, end: Z")
                    parts = identifier.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly using _get_sequence_for_position
                    sequence = self._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        result = {
                            'sequence': sequence,
                            'start': start,
                            'end': end,
                            'chrom_length': len(sequence)
                        }
                        self._sequence_cache[cache_key] = result
                        self.logger.debug(f"Retrieved and cached position sequence ({len(sequence)} bp)")
                        return result
                        
                    self.logger.warning(f"No sequence found for position {chrom}:{start}-{end}")
                    return None
                    
                except Exception as e:
                    self.logger.error(f"Error parsing position or getting sequence: {str(e)}")
                    return None
            else:
                # Regular gene-based search
                self.logger.debug(f"Getting gene data for locus tag: {identifier}")
                gene_data = self.get_gene_data(identifier)
                if not gene_data or 'info' not in gene_data:
                    self.logger.warning(f"No gene data found for locus tag: {identifier}")
                    return None
                
                # Parse location string (format: "start:end(strand)")
                location = gene_data['info']['location']
                if ':' not in location:
                    self.logger.warning(f"Invalid location format: {location}")
                    return None
                
                # Extract start and end positions
                start = int(location.split(':')[0])
                end = int(location.split(':')[1].split('(')[0])
                
                # Get sequence from gene_data directly if available
                if 'sequence' in gene_data:
                    sequence = gene_data['sequence']
                    self.logger.debug(f"Got sequence of length: {len(sequence)}")
                    
                    # Format sequence with padding in lowercase
                    padding = 30
                    padded_start = max(0, start - padding)
                    padded_end = min(len(sequence), end + padding)
                    
                    # Split sequence into parts
                    five_prime_pad = sequence[:start - padded_start].lower() if start > padded_start else ""
                    main_sequence = sequence[start - padded_start:end - padded_start].upper()
                    three_prime_pad = sequence[end - padded_start:].lower()
                    
                    # Combine parts
                    formatted_sequence = five_prime_pad + main_sequence + three_prime_pad
                    
                    # Cache the result
                    result = {
                        'sequence': formatted_sequence,
                        'chrom_length': len(sequence),
                        'start': start,
                        'end': end,
                        'padded_start': padded_start,
                        'padded_end': padded_end
                    }
                    self._sequence_cache[cache_key] = result
                    
                    self.logger.debug(f"Retrieved and cached sequence for locus tag {identifier} ({len(formatted_sequence)} bp)")
                    return result
                    
                self.logger.warning(f"No sequence data found in gene_data for {identifier}")
                return None
                
        except Exception as e:
            self.logger.error(f"Error getting gene sequence: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return None

    def _get_sequence_for_position(self, chrom, start, end):
        """Get sequence for a given position with proper padding handling"""
        try:
            if not hasattr(self, 'annotation_parser') or self.annotation_parser is None:
                self.annotation_parser = AnnotationParser(self.global_settings)
                annotation_file = self.global_settings.get_current_annotation_file()
                annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
                self.annotation_parser.set_annotation_file(annotation_path)
                
            # Get the full chromosome ID by counting carets in annotation file
            full_chrom = None
            chrom_count = 0
            
            try:
                for record in SeqIO.parse(self.annotation_path, "genbank"):
                    chrom_count += 1
                    if chrom_count == int(chrom):  # Match based on position rather than ID number
                        full_chrom = record.id
                        self.logger.debug(f"Found chromosome {chrom} as {full_chrom}")
                        break
            except Exception as e:
                self.logger.error(f"Error finding chromosome by position: {str(e)}")
                return None

            if not full_chrom:
                self.logger.warning(f"Could not find chromosome at position {chrom}")
                return None

            feature_info = {
                'chromosome': full_chrom,
                'start': start-1,
                'end': end
            }
            
            self.logger.debug(f"Getting sequence for feature info: {feature_info}")
            
            sequence = self.annotation_parser._get_sequence_for_gene(feature_info)
            if sequence:
                padding = 30
                
                # Handle start position padding
                if start == 1:
                    # No padding at start if starting at position 1
                    five_prime_pad = ""
                    main_sequence = sequence[:-(padding if len(sequence) > padding else 0)].upper()
                else:
                    five_prime_pad = sequence[:padding].lower() if len(sequence) > padding else ""
                    main_sequence = sequence[padding:-padding].upper() if len(sequence) > 60 else sequence.upper()
                
                three_prime_pad = sequence[-padding:].lower() if len(sequence) > padding else ""
                
                return five_prime_pad + main_sequence + three_prime_pad
                
            return None
            
        except Exception as e:
            self.logger.error(f"Error getting sequence for position: {str(e)}")
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

    def get_gene_sequence_for_range(self, identifier, start, end):
        try:
            # For feature-based searches
            gene_data = self.get_gene_data(identifier)
            if not gene_data or 'info' not in gene_data:
                self.logger.warning(f"No gene data found for identifier: {identifier}")
                return None
                
            # Get chromosome from gene data
            chrom = gene_data['info']['chromosome'].split('.')[-1]  # Extract chromosome number
            
            # Use _get_sequence_for_position to get sequence with padding
            sequence = self._get_sequence_for_position(int(chrom), start, end)
            if sequence:
                result = {
                    'sequence': sequence,
                    'start': start,
                    'end': end
                }
                return result
                
            self.logger.warning(f"No sequence found for range {start}-{end} in chromosome {chrom}")
            return None
            
        except Exception as e:
            self.logger.error(f"Error getting gene sequence for range: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return None