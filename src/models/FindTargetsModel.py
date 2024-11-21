import time
from models.HomeWindowModel import HomeWindowModel
from models.CSPRparser import CSPRparser
from models.AnnotationParser import AnnotationParser
import os
from functools import lru_cache
import traceback
from Bio import SeqIO

class FindTargetsModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.results = {}
        self._parser_cache = {} 
        self.global_settings.annotation_file_changed.connect(self._on_annotation_file_changed)

    def _on_annotation_file_changed(self, new_annotation_file):
        """Clear caches when annotation file changes"""
        self.global_settings.logger.debug(f"FindTargetsModel clearing caches for new annotation file: {new_annotation_file}")
        self._parser_cache.clear()

    @lru_cache(maxsize=32)
    def _get_parser(self, file_path):
        """Cache CSPRparser instances for reuse"""
        if file_path not in self._parser_cache:
            self._parser_cache[file_path] = CSPRparser(file_path, self.global_settings.get_casper_info_path())
        return self._parser_cache[file_path]

    def find_targets(self, input_data):
        self.global_settings.logger.debug(f"Received input data: {input_data}")
        
        organism = input_data['organism']
        endo = input_data['endonuclease']
        org_files = self.get_organism_to_files()

        self._validate_input(organism, endo, org_files)
        
        file_path = os.path.join(self.global_settings.get_db_path(), org_files[organism][endo][0])
        parser = self._get_parser(file_path)
        
        search_types = {
            'feature': self.find_targets_by_feature,
            'position': self.find_targets_by_position,
            'sequence': self.find_targets_by_sequence
        }
        
        search_func = search_types.get(input_data['search_type'])
        if not search_func:
            error_msg = f"Invalid search type: {input_data['search_type']}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)
        
        self.results = search_func(parser, input_data)
        
        return self.results

    def _validate_input(self, organism, endo, org_files):
        if organism not in org_files:
            error_msg = f"Organism '{organism}' not found in the database. Available organisms: {list(org_files.keys())}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)
        
        if endo not in org_files[organism]:
            error_msg = f"Endonuclease '{endo}' not found for organism '{organism}'. Available endonucleases: {list(org_files[organism].keys())}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)

    def find_targets_by_feature(self, parser, input_data):
        try:
            annotation_file = (input_data.get('annotation_file') or 
                             self.global_settings.get_current_annotation_file())
            
            search_query = input_data['search_query'].strip()
            
            annotation_parser = AnnotationParser(self.global_settings)
            annotation_file_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
            annotation_parser.set_annotation_file(annotation_file_path)
            
            results_list = annotation_parser.genbank_search([search_query])
            
            formatted_results = []
            
            chrom_mapping = {}
            chrom_count = 0
            for record in SeqIO.parse(annotation_file_path, "genbank"):
                chrom_count += 1
                chrom_mapping[record.id] = str(chrom_count)
            
            for record_id, feature_info in results_list:
                location = feature_info['feature_location']
                start_end = location.split('(')[0]  # Get part before the strand
                start, end = map(int, start_end.split(':'))
                
                chrom_num = chrom_mapping.get(record_id, '1')
                
                target_info = {
                    'feature_type': 'CDS',
                    'chromosome': chrom_num,
                    'full_chromosome': record_id,
                    'feature_id': feature_info['feature_id'],
                    'feature_name': feature_info['feature_name'],
                    'feature_description': feature_info['feature_description'],
                    'location': f"{start}-{end}",
                    'start': start,
                    'end': end,
                    'strand': '+' if '(+)' in location else '-',
                    'endonuclease': input_data['endonuclease']
                }
                
                self.global_settings.logger.debug(f"Created target info: {target_info}")
                
                formatted_results.append(target_info)
                
            return formatted_results
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in find_targets_by_feature: {str(e)}")
            raise

    def find_targets_by_position(self, parser, input_data):
        try:
            queries = input_data['search_query'].strip().split('\n')
            all_results = []
            
            for query in queries:
                try:
                    chrom, start, end = map(int, query.strip().split(','))
                    
                    # Get full chromosome ID by counting carets
                    full_chrom = None
                    chrom_count = 0
                    
                    # Get annotation file path
                    annotation_file = self.global_settings.get_current_annotation_file()
                    annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
                    
                    # Find the full chromosome ID by position
                    for record in SeqIO.parse(annotation_path, "genbank"):
                        chrom_count += 1
                        if chrom_count == chrom:  # Match based on position
                            full_chrom = record.id
                            self.logger.debug(f"Found chromosome {chrom} as {full_chrom}")
                            break
                    
                    if not full_chrom:
                        self.logger.warning(f"Could not find chromosome at position {chrom}")
                        continue
                    
                    # Create target info with proper formatting
                    position_name = f"chrom {chrom}, start: {start}, end: {end}"
                    target_info = [{
                        'start': start,
                        'end': end,
                        'feature_id': position_name,
                        'feature_name': position_name,
                        'chromosome': str(chrom),  # Keep chromosome number for CSPR lookup
                        'full_chromosome': full_chrom  # Store full ID for sequence lookup
                    }]
                    
                    # Get targets using batch processing
                    self.logger.debug(f"Searching for targets in chromosome {chrom} from {start} to {end}")
                    targets = parser.read_targets_batch(str(chrom), target_info, input_data['endonuclease'])
                    
                    if targets:
                        self.logger.debug(f"Found {len(targets)} raw targets")
                        filtered_targets = []
                        guide_length = 23  # Length of guide RNA
                        
                        for target in targets:
                            target_pos = int(target['position'])
                            target_end = target_pos 
                            
                            # Include target if:
                            # 1. Target start position is within range
                            # 2. Target end position is within or equal to end position
                            if start <= target_pos and target_end <= end + 1:
                                filtered_targets.append(target)
                        
                        self.logger.debug(f"Filtered to {len(filtered_targets)} targets within range")
                        
                        # Get sequence for this region
                        sequence = self._get_sequence_for_position(chrom, start, end)
                        
                        # Format results
                        for target in filtered_targets:
                            result = {
                                'feature_type': 'Position',
                                'chromosome': str(chrom),
                                'feature_id': position_name,
                                'feature_name': position_name,
                                'feature_description': position_name,
                                'location': target['location'],
                                'start': start,
                                'end': end,
                                'strand': target['strand'],
                                'sequence': target['sequence'],
                                'pam': target['pam'],
                                'score': target['score'],
                                'endonuclease': target['endonuclease'],
                                'gene_sequence': sequence
                            }
                            all_results.append(result)
                            
                        self.logger.debug(f"Added {len(filtered_targets)} formatted results")
                    else:
                        self.logger.warning(f"No targets found for query: {query}")
                        
                except Exception as e:
                    self.logger.error(f"Error processing query {query}: {str(e)}")
                    self.logger.error(f"Stack trace: {traceback.format_exc()}")
                    continue
            
            self.logger.debug(f"Total results found: {len(all_results)}")
            return all_results
            
        except Exception as e:
            self.logger.error(f"Error in find_targets_by_position: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            raise

    def _get_sequence_for_position(self, chrom, start, end):
        """Get sequence for a given position with proper padding handling"""
        try:
            if not hasattr(self, 'annotation_parser') or self.annotation_parser is None:
                self.annotation_parser = AnnotationParser(self.global_settings)
                annotation_file = self.global_settings.get_current_annotation_file()
                annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
                self.annotation_parser.set_annotation_file(annotation_path)
                
            feature_info = {
                'chromosome': f"NZ_CP032679.{chrom}",  # Use full name 'chromosome'
                'start': start-1,  # Use full name 'start'
                'end': end  # Use full name 'end'
            }
            
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

    def find_targets_by_sequence(self, parser, input_data):
        """Search for targets by sequence"""
        try:
            sequence = input_data['search_query'].strip().upper()
            
            # Validate sequence length
            if len(sequence) < 100:
                self.logger.error("Sequence too short")
                return []
                
            # Get annotation file
            annotation_file = self.global_settings.get_current_annotation_file()
            if not annotation_file:
                self.logger.error("No annotation file selected")
                return []
                
            # Initialize annotation parser
            if not hasattr(self, 'annotation_parser') or self.annotation_parser is None:
                self.annotation_parser = AnnotationParser(self.global_settings)
                annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
                self.annotation_parser.set_annotation_file(annotation_path)
                
            # Find sequence in genome
            chrom_count = 0
            for record in SeqIO.parse(self.annotation_parser.annotation_file_name, "genbank"):
                chrom_count += 1  # Count chromosome position by caret
                record_seq = str(record.seq).upper()
                pos = record_seq.find(sequence)
                
                if pos != -1:
                    # Found the sequence
                    start = pos + 1  # 1-based position
                    end = start + len(sequence) - 1
                    
                    # Create position name
                    position_name = f"chrom {chrom_count}, start: {start}, end: {end}"
                    
                    # Create target info
                    target_info = [{
                        'start': start,
                        'end': end,
                        'feature_id': position_name,
                        'feature_name': position_name,
                        'chromosome': str(chrom_count),  # Use caret-based chromosome number
                        'full_chromosome': record.id  # Store full chromosome ID
                    }]
                    
                    # Get targets in this region
                    self.logger.debug(f"Found sequence in chromosome {chrom_count} from {start} to {end}")
                    targets = parser.read_targets_batch(str(chrom_count), target_info, input_data['endonuclease'])
                    
                    if targets:
                        self.logger.debug(f"Found {len(targets)} raw targets")
                        filtered_targets = []
                        guide_length = 23
                        
                        for target in targets:
                            target_pos = int(target['position'])
                            target_end = target_pos + guide_length
                            
                            # Include target if within sequence bounds
                            if start <= target_pos and target_end <= end + 1:
                                filtered_targets.append(target)
                        
                        self.logger.debug(f"Filtered to {len(filtered_targets)} targets within range")
                        
                        # Get sequence with padding
                        sequence_with_padding = self._get_sequence_for_position(chrom_count, start, end)
                        
                        # Format results
                        all_results = []
                        for target in filtered_targets:
                            result = {
                                'feature_type': 'Position',
                                'chromosome': str(chrom_count),
                                'feature_id': position_name,
                                'feature_name': position_name,
                                'feature_description': f"Sequence match at {position_name}",
                                'location': target['location'],
                                'start': start,
                                'end': end,
                                'strand': target['strand'],
                                'sequence': target['sequence'],
                                'pam': target['pam'],
                                'score': target['score'],
                                'endonuclease': target['endonuclease'],
                                'gene_sequence': sequence_with_padding
                            }
                            all_results.append(result)
                            
                        return all_results
                        
            self.logger.warning("Sequence not found in genome")
            return []
            
        except Exception as e:
            self.logger.error(f"Error in find_targets_by_sequence: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            raise

    def _format_results(self, targets):
        formatted_results = []
        for target in targets:
            formatted_results.append({
                'location': target[0],
                'sequence': target[1],
                'pam': target[2],
                'score': target[3],
                'strand': target[4],
                'endonuclease': target[5]
            })
        return formatted_results
