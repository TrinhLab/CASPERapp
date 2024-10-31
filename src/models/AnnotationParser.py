from PyQt6.QtWidgets import QMessageBox
from Bio import SeqIO
import os
import traceback
from functools import lru_cache

class AnnotationParser:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.annotation_file_name = ""
        self.available_genes = []
        self._feature_cache = {}  # Cache for feature data
        self._record_cache = {}   # Cache for SeqIO records
        self.gene_cache = {}  # Add cache for gene data

    def set_annotation_file(self, file_path):
        if self.annotation_file_name != file_path:
            self.annotation_file_name = file_path
            self.logger.debug(f"Set annotation file to: {file_path}")
            self._feature_cache.clear()  # Clear cache when file changes
            self._record_cache.clear()
            if hasattr(self, '_gene_index'):
                delattr(self, '_gene_index')
            
            # Pre-load records and build index
            records = self._get_records()
            self._build_gene_index(records)
            self._parse_available_genes()

    @lru_cache(maxsize=1)
    def _get_records(self):
        """Cache and return all records from the annotation file"""
        if not self._record_cache:
            try:
                self._record_cache = list(SeqIO.parse(self.annotation_file_name, "genbank"))
            except Exception as e:
                self.logger.error(f"Error reading annotation file: {str(e)}")
                return []
        return self._record_cache

    def genbank_search(self, queries):
        try:
            if not self.annotation_file_name:
                raise ValueError("Annotation file not set")
            
            self.logger.debug(f"Searching in annotation file: {self.annotation_file_name}")
            results_list = []
            
            # Convert queries to lowercase set for faster lookup
            queries = {q.lower() for q in queries}
            
            # Use cached records
            for record in self._get_records():
                for feature in record.features:
                    if feature.type in ['CDS', 'gene']:
                        # Create a hashable cache key using feature start and end positions
                        cache_key = (record.id, feature.type, 
                                   str(feature.location.start), 
                                   str(feature.location.end))
                        
                        # Use cached feature info if available
                        if cache_key not in self._feature_cache:
                            self._feature_cache[cache_key] = self._get_feature_info(feature)
                        
                        feature_info = self._feature_cache[cache_key]
                        
                        # Combine searchable text for single comparison
                        searchable_text = ' '.join([
                            feature_info['feature_name'].lower(),
                            feature_info['feature_id'].lower(),
                            feature_info['feature_description'].lower()
                        ])
                        
                        # Check if any query matches
                        if any(query in searchable_text for query in queries):
                            results_list.append((record.id, feature))
            
            self.logger.debug(f"Found {len(results_list)} results")
            return results_list
        except Exception as e:
            self.logger.error(f"Error in genbank_search: {str(e)}")
            raise

    def get_max_chrom(self):
        try:
            parser = SeqIO.parse(self.annotation_file_name, 'genbank')
            max_chrom = sum(1 for _ in parser)
            return max_chrom
        except Exception as e:
            self.logger.error(f"Error in get_max_chrom: {str(e)}")
            self._show_error("Error in get_max_chrom", str(e))
            return 0

    def get_sequence_info(self, query):
        # Implement this method if needed
        pass

    def find_which_file_version(self):
        try:
            if not self.annotation_file_name or os.path.basename(self.annotation_file_name) == "None":
                return -1
            if self.annotation_file_name.endswith(('.gbff', '.gbk')):
                return "gbff"
            else:
                return -1
        except Exception as e:
            self.logger.error(f"Error in find_which_file_version: {str(e)}")
            self._show_error("Error in find_which_file_version", str(e))
            return -1

    def _show_error(self, title, message):
        QMessageBox.critical(None, title, f"{message}\n\nFor more information, check the log file.")

    @staticmethod
    def flatten_list(t):
        return [item.lower() for sublist in t for item in sublist]

    def _get_feature_info(self, feature):
        return {
            'feature_id': self._get_feature_id(feature),
            'feature_name': self._get_feature_name(feature),
            'feature_location': self._get_feature_location(feature),
            'feature_description': self._get_feature_description(feature)
        }

    def _get_feature_id(self, feature):
        for key in ['locus_tag']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def _get_feature_name(self, feature):
        for key in ['gene']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"
    
    def _get_feature_location(self, feature):
        if feature.location:
            start = feature.location.start
            end = feature.location.end
            strand = '+' if feature.location.strand == 1 else '-'
            return f"{start}:{end}({strand})"
        return "N/A"

    def _get_feature_description(self, feature):
        for key in ['product', 'note']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def get_available_genes(self):
        return self.available_genes

    def get_gene_data(self, gene_identifier):
        """
        Get gene data using gene name or locus tag with optimized caching
        """
        try:
            self.logger.debug(f"AnnotationParser.get_gene_data called with identifier: {gene_identifier}")
            
            if not gene_identifier:
                self.logger.warning("Empty gene identifier provided")
                return None
                
            # Handle numeric gene identifiers
            if isinstance(gene_identifier, int) or str(gene_identifier).isdigit():
                if self.available_genes:
                    gene_identifier = self.available_genes[0]
                else:
                    return None

            # Check main cache first
            cache_key = f"gene_data_{gene_identifier}"
            if cache_key in self._feature_cache:
                return self._feature_cache[cache_key]

            # Get cached records
            records = self._get_records()
            if not records:
                return None

            # Use gene index if available
            if not hasattr(self, '_gene_index'):
                self._build_gene_index(records)

            # Try to get location from index
            if gene_identifier in self._gene_index:
                record_id, feature = self._gene_index[gene_identifier]
                for record in records:
                    if record.id == record_id:
                        sequence = str(feature.extract(record.seq))
                        feature_info = self._get_feature_info(feature)
                        
                        result = {
                            'sequence': sequence,
                            'info': feature_info
                        }
                        
                        self._feature_cache[cache_key] = result
                        return result

            return None
            
        except Exception as e:
            self.logger.error(f"Error in get_gene_data: {str(e)}")
            return None

    def _build_gene_index(self, records):
        """Build an index of genes for faster lookup"""
        self._gene_index = {}
        try:
            for record in records:
                for feature in record.features:
                    if feature.type == 'gene':
                        gene_name = self._get_feature_name(feature)
                        gene_id = self._get_feature_id(feature)
                        if gene_name != "N/A":
                            self._gene_index[gene_name] = (record.id, feature)
                        if gene_id != "N/A":
                            self._gene_index[gene_id] = (record.id, feature)
        except Exception as e:
            self.logger.error(f"Error building gene index: {str(e)}")

    def _parse_available_genes(self):
        self.available_genes = []
        try:
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                for feature in record.features:
                    if feature.type == 'gene':
                        self.available_genes.append(self._get_feature_name(feature))
        except Exception as e:
            self.logger.error(f"Error parsing available genes: {str(e)}")

    def get_full_gene_sequence(self):
        # Implement this method if needed
        pass
