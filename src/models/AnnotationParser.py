from PyQt6.QtWidgets import QMessageBox
from Bio import SeqIO
import os
import traceback
from functools import lru_cache
import json
import pickle
import time

class AnnotationParser:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.annotation_file_name = ""
        self.available_genes = []
        self._feature_cache = {}
        self._record_cache = {}
        self.gene_cache = {}
        self.index_file = None

    def set_annotation_file(self, file_path):
        """Set the annotation file and initialize/load index"""
        try:
            # Don't process if file_path is a directory or empty
            if not file_path or os.path.isdir(file_path):
                self.logger.debug(f"Invalid annotation file path: {file_path}")
                self._index = {'locus_tags': {}}  # Initialize empty index
                return

            if self.annotation_file_name != file_path:
                total_start = time.time()
                
                self.annotation_file_name = file_path
                self.logger.debug(f"Set annotation file to: {file_path}")
                
                # Set index file path
                self.index_file = f"{file_path}.index"
                
                # Load or create index
                index_start = time.time()
                if not self._load_index():
                    self.logger.debug("Index not found or outdated, creating new index...")
                    create_start = time.time()
                    self._create_index()
                    create_time = time.time() - create_start
                    self.logger.debug(f"Index creation time: {create_time:.2f} seconds")
                index_time = time.time() - index_start
                self.logger.debug(f"Total index handling time: {index_time:.2f} seconds")
                
        except Exception as e:
            self.logger.error(f"Error in set_annotation_file: {str(e)}")
            raise

    def _create_index(self):
        try:
            start_time = time.time()
            self.logger.debug("Creating gene index file...")
            
            # Initialize optimized index structure - no sequences stored
            index_data = {
                'locus_tags': {},  # Only store essential data
            }
            
            # Process records
            record_count = 0
            feature_count = 0
            
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                record_count += 1
                record_start = time.time()
                
                # Process features
                for feature in record.features:
                    if feature.type in ['CDS', 'gene']:
                        feature_count += 1
                        
                        # Get essential feature info
                        locus_tag = None
                        if 'locus_tag' in feature.qualifiers:
                            locus_tag = feature.qualifiers['locus_tag'][0]
                        elif 'gene' in feature.qualifiers:
                            locus_tag = feature.qualifiers['gene'][0]
                        
                        # Only process features with valid locus tags
                        if locus_tag and locus_tag.lower() != "n/a":
                            # Get location info
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            strand = '+' if feature.location.strand == 1 else '-'
                            
                            # Store feature info with full names
                            feature_entry = {
                                'feature_type': feature.type,
                                'chromosome': record.id,
                                'location': f"{start}:{end}({strand})",
                                'gene_name': feature.qualifiers.get('gene', ['N/A'])[0],
                                'description': feature.qualifiers.get('product', 
                                    feature.qualifiers.get('note', ['N/A']))[0],
                                'start': start,
                                'end': end
                            }
                            
                            # Store in index
                            index_data['locus_tags'][locus_tag] = feature_entry
                
                record_time = time.time() - record_start
                if record_count % 100 == 0:
                    self.logger.debug(f"Processed {record_count} records, {feature_count} features. Last record time: {record_time:.2f}s")
            
            # Save compressed index to file
            save_start = time.time()
            with open(self.index_file, 'wb') as f:
                pickle.dump(index_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            save_time = time.time() - save_start
            total_time = time.time() - start_time
            
            self._index = index_data
            
            self.logger.debug(f"Index creation complete. Records: {record_count}, Features: {feature_count}")
            self.logger.debug(f"Save time: {save_time:.2f}s, Total time: {total_time:.2f}s")
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating index: {str(e)}")
            return False

    def _load_index(self):
        """Load the index file if it exists and is newer than the GenBank file"""
        try:
            if not os.path.exists(self.index_file):
                return False
                
            # Check if index is older than GenBank file
            if os.path.getmtime(self.index_file) < os.path.getmtime(self.annotation_file_name):
                return False

            start_time = time.time()
            with open(self.index_file, 'rb') as f:
                self._index = pickle.load(f)
            load_time = time.time() - start_time
            print(f"Index file: {self._index}")
            self.logger.debug(f"Index file loaded successfully in {load_time:.2f} seconds")
            return True
            
        except Exception as e:
            self.logger.error(f"Error loading index: {str(e)}")
            return False

    def genbank_search(self, queries):
        """Search using the index file for better performance"""
        try:
            if not self.annotation_file_name or os.path.isdir(self.annotation_file_name):
                self.logger.warning("No valid annotation file set")
                return []
            
            self.logger.debug(f"Searching in annotation file: {self.annotation_file_name}")
            results_list = []
            
            # Convert queries to lowercase set for faster lookup
            queries = {q.lower() for q in queries}
            self.logger.debug(f"Search queries: {queries}")
            
            # Search through index
            if hasattr(self, '_index') and 'locus_tags' in self._index:
                # Search through features, filtering for CDS and gene types only
                for locus_tag, feature_entry in self._index['locus_tags'].items():
                    # Safely get feature type with default value
                    feature_type = feature_entry.get('feature_type', '')
                    
                    # Only process CDS and gene features
                    if feature_type not in ['CDS', 'gene']:
                        continue
                        
                    # Check gene name, locus tag, and description
                    searchable_text = ' '.join([
                        feature_entry.get('gene_name', '').lower(),
                        locus_tag.lower(),
                        feature_entry.get('description', '').lower()
                    ])
                    
                    # Check if any query matches
                    if any(query in searchable_text for query in queries):
                        info = {
                            'feature_id': locus_tag,
                            'feature_name': feature_entry.get('gene_name', 'N/A'),
                            'feature_location': feature_entry.get('location', 'N/A'),
                            'feature_description': feature_entry.get('description', 'N/A')
                        }
                        results_list.append((feature_entry.get('chromosome', ''), info))
                        
            return results_list
            
        except Exception as e:
            self.logger.error(f"Error in genbank_search: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")  # Add stack trace for better debugging
            raise

    def get_gene_data(self, gene_identifier):
        """Get gene data using the optimized index and fetch sequence on demand"""
        try:
            if not gene_identifier:
                return None
                
            # Ensure string conversion and proper formatting
            gene_identifier = str(gene_identifier).strip().lower()
            
            if hasattr(self, '_index'):
                # Try exact match first
                if gene_identifier in self._index['locus_tags']:
                    gene_info = self._index['locus_tags'][gene_identifier]
                    
                    # Get sequence from file
                    sequence = self._get_sequence_for_gene(gene_info)
                    if sequence is None:
                        return None
                    
                    # Use full names instead of shortened keys
                    expanded_info = {
                        'feature_type': gene_info['feature_type'],
                        'chromosome': gene_info['chromosome'],
                        'location': gene_info['location'],
                        'gene_name': gene_info['gene_name'],
                        'description': gene_info['description'],
                        'start': gene_info['start'],
                        'end': gene_info['end']
                    }
                    
                    return {
                        'sequence': sequence,
                        'info': expanded_info
                    }
                        
                # Try case-insensitive match
                for key, value in self._index['locus_tags'].items():
                    if str(key).lower() == gene_identifier:
                        sequence = self._get_sequence_for_gene(value)
                        if sequence is None:
                            return None
                        
                        expanded_info = {
                            'feature_type': value['feature_type'],
                            'chromosome': value['chromosome'],
                            'location': value['location'],
                            'gene_name': value['gene_name'],
                            'description': value['description'],
                            'start': value['start'],
                            'end': value['end']
                        }
                        
                        return {
                            'sequence': sequence,
                            'info': expanded_info
                        }
                        
                return None
                    
        except Exception as e:
            self.logger.error(f"Error in get_gene_data: {str(e)}")
            return None

    def _get_sequence_for_gene(self, gene_info):
        """Get sequence for a gene from the GenBank file"""
        try:
            self.logger.debug(f"Getting sequence for gene info: {gene_info} in _get_sequence_for_gene")
            # Parse the GenBank file and find the right record
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                if record.id == gene_info['chromosome']:  # Use full chromosome name
                    sequence = str(record.seq)
                    
                    # Get sequence with padding
                    padding = 30
                    start = max(0, gene_info['start'] - padding)
                    end = min(len(sequence), gene_info['end'] + padding)
                    padded_sequence = sequence[start:end]

                    self.logger.debug(f"Padded sequence: {padded_sequence}")
                    
                    return padded_sequence
                    
            return None
            
        except Exception as e:
            self.logger.error(f"Error getting sequence for gene: {str(e)}")
            return None

    @lru_cache(maxsize=1)
    def _get_records(self):
        """Cache and return all records from the annotation file"""
        start_time = time.time()
        if not self._record_cache:
            try:
                self.logger.debug("Loading records from file...")
                self._record_cache = list(SeqIO.parse(self.annotation_file_name, "genbank"))
                load_time = time.time() - start_time
                self.logger.debug(f"Time to load records: {load_time:.2f} seconds")
            except Exception as e:
                self.logger.error(f"Error reading annotation file: {str(e)}")
                return []
        return self._record_cache

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

    def get_full_gene_sequence(self):
        # Implement this method if needed
        pass

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

