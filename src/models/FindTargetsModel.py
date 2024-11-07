import time
from models.HomeWindowModel import HomeWindowModel
from models.CSPRparser import CSPRparser
from models.AnnotationParser import AnnotationParser
import os
from functools import lru_cache

class FindTargetsModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.results = {}
        self._parser_cache = {}  # Cache for CSPRparser instances
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
        
        start_time = time.time()
        
        organism = input_data['organism']
        endo = input_data['endonuclease']
        org_files = self.get_organism_to_files()

        # Validate input data
        validate_start = time.time()
        self._validate_input(organism, endo, org_files)
        validate_time = time.time() - validate_start
        self.global_settings.logger.debug(f"Validation time: {validate_time:.2f} seconds")
        
        # Get file path and parser
        parser_start = time.time()
        file_path = os.path.join(self.global_settings.get_db_path(), org_files[organism][endo][0])
        parser = self._get_parser(file_path)
        parser_time = time.time() - parser_start
        self.global_settings.logger.debug(f"Parser initialization time: {parser_time:.2f} seconds")
        
        # Use dictionary for faster lookup
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
        
        # Perform the search
        search_start = time.time()
        self.results = search_func(parser, input_data)
        search_time = time.time() - search_start
        self.global_settings.logger.debug(f"Search execution time: {search_time:.2f} seconds")
        
        total_time = time.time() - start_time
        self.global_settings.logger.debug(f"Total find_targets time: {total_time:.2f} seconds")
        
        return self.results

    def _validate_input(self, organism, endo, org_files):
        """Validate input parameters"""
        if organism not in org_files:
            error_msg = f"Organism '{organism}' not found in the database. Available organisms: {list(org_files.keys())}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)
        
        if endo not in org_files[organism]:
            error_msg = f"Endonuclease '{endo}' not found for organism '{organism}'. Available endonucleases: {list(org_files[organism].keys())}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)

    def find_targets_by_feature(self, parser, input_data):
        """Search for features using the indexed annotation parser"""
        try:
            start_time = time.time()
            
            # Get annotation file from input data or global settings
            annotation_file = (input_data.get('annotation_file') or 
                             self.global_settings.get_current_annotation_file())
            
            search_query = input_data['search_query'].strip()
            
            # Create new annotation parser instance
            parser_start = time.time()
            annotation_parser = AnnotationParser(self.global_settings)
            annotation_file_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
            annotation_parser.set_annotation_file(annotation_file_path)
            parser_time = time.time() - parser_start
            self.global_settings.logger.debug(f"Annotation parser initialization time: {parser_time:.2f} seconds")
            
            # Use indexed search
            search_start = time.time()
            results_list = annotation_parser.genbank_search([search_query])
            search_time = time.time() - search_start
            self.global_settings.logger.debug(f"Genbank search time: {search_time:.2f} seconds")
            
            # Format results
            format_start = time.time()
            formatted_results = []
            for record_id, feature_info in results_list:
                # Extract start and end from feature_location
                location = feature_info['feature_location']
                start_end = location.split('(')[0]  # Get part before the strand
                start, end = map(int, start_end.split(':'))
                
                # Extract chromosome number from record_id (e.g., "NZ_CP132594.1" -> "1")
                chrom_num = record_id.split('.')[-1] if '.' in record_id else '1'
                
                # Create target info with feature_id and chromosome number
                target_info = {
                    'feature_type': 'CDS',
                    'chromosome': chrom_num,  # Use chromosome number
                    'full_chromosome': record_id,  # Store full chromosome name for reference
                    'feature_id': feature_info['feature_id'],
                    'feature_name': feature_info['feature_name'],
                    'feature_description': feature_info['feature_description'],
                    'location': f"{start}-{end}",
                    'start': start,
                    'end': end,
                    'strand': '+' if '(+)' in location else '-',
                    'endonuclease': input_data['endonuclease']
                }
                
                # Debug log the target info
                self.global_settings.logger.debug(f"Created target info: {target_info}")
                
                formatted_results.append(target_info)
                
            format_time = time.time() - format_start
            self.global_settings.logger.debug(f"Result formatting time: {format_time:.2f} seconds")
            
            # Debug log sample results
            if formatted_results:
                self.global_settings.logger.debug(f"Sample formatted result: {formatted_results[0]}")
                self.global_settings.logger.debug(f"Feature IDs present: {[r['feature_id'] for r in formatted_results[:5]]}")
            
            total_time = time.time() - start_time
            self.global_settings.logger.debug(f"Total find_targets_by_feature time: {total_time:.2f} seconds")
            
            return formatted_results
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in find_targets_by_feature: {str(e)}")
            raise

    def find_targets_by_position(self, parser, input_data):
        search_query = input_data['search_query']
        chrom, start, end = map(int, search_query.split(','))
        pos_tuple = (chrom, start - 1, end)  # Adjust for 0-based indexing
        
        targets = parser.read_targets(f"position_{chrom}_{start}_{end}", pos_tuple, input_data['endonuclease'])
        return self._format_results(targets)

    def find_targets_by_sequence(self, parser, input_data):
        search_query = input_data['search_query'].upper()
        annotation_file = input_data['annotation_file']
        
        self.annotation_parser.annotationFileName = os.path.join(self.global_settings.get_db_path(), annotation_file)
        sequence_info = self.annotation_parser.get_sequence_info(search_query)
        
        if sequence_info:
            chrom, start, end = sequence_info
            pos_tuple = (chrom, start - 1, end)
            targets = parser.read_targets(f"sequence_{start}_{end}", pos_tuple, input_data['endonuclease'])
            return self._format_results(targets)
        else:
            return []

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
