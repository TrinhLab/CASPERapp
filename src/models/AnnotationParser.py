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
        try:
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
            
            # Initialize index structure
            index_data = {
                'locus_tags': {},  # Only store by locus_tag
                'sequences': {}  # Keep sequences for quick access
            }
            
            # Process records
            record_count = 0
            feature_count = 0
            
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                record_count += 1
                record_start = time.time()
                
                # Store sequence information first
                index_data['sequences'][record.id] = str(record.seq)
                
                # Process features
                for feature in record.features:
                    if feature.type in ['CDS', 'gene']:
                        feature_count += 1
                        feature_info = self._get_feature_info(feature)
                        locus_tag = feature_info['feature_id']
                        
                        # Only create feature entry if we have a valid locus_tag
                        if locus_tag and locus_tag.lower() != "n/a":
                            feature_entry = {
                                'record_id': record.id,
                                'feature_type': feature.type,
                                'chromosome': record.id,
                                'location': self._get_feature_location(feature),
                                'strand': '+' if feature.location.strand == 1 else '-',
                                'locus_tag': locus_tag,
                                'gene_name': feature_info['feature_name'],
                                'description': feature_info['feature_description'],
                                'qualifiers': {k: v[0] if isinstance(v, list) else v 
                                             for k, v in feature.qualifiers.items()}
                            }
                            
                            # Index only by locus_tag (lowercase for case-insensitive lookup)
                            index_data['locus_tags'][locus_tag.lower()] = feature_entry
                            
            record_time = time.time() - record_start
            if record_count % 100 == 0:
                self.logger.debug(f"Processed {record_count} records, {feature_count} features. Last record time: {record_time:.2f}s")
            
            # Save index to file
            save_start = time.time()
            with open(self.index_file, 'wb') as f:
                pickle.dump(index_data, f)
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
            self.logger.debug(f"Index file loaded successfully in {load_time:.2f} seconds")
            return True
            
        except Exception as e:
            self.logger.error(f"Error loading index: {str(e)}")
            return False

    def genbank_search(self, queries):
        """Search using the index file for better performance"""
        try:
            if not self.annotation_file_name:
                raise ValueError("Annotation file not set")
            
            self.logger.debug(f"Searching in annotation file: {self.annotation_file_name}")
            results_list = []
            
            # Convert queries to lowercase set for faster lookup
            queries = {q.lower() for q in queries}
            print(f"Search queries: {queries}")
            
            # Search through index
            if hasattr(self, '_index'):
                # Search through all features
                for feature_key, feature_entry in self._index['locus_tags'].items():
                    # Check gene name, locus tag, and description
                    searchable_text = ' '.join([
                        feature_entry['gene_name'].lower(),
                        feature_entry['locus_tag'].lower(),
                        feature_entry['description'].lower(),
                        # Also search through qualifiers
                        *[str(v).lower() for v in feature_entry['qualifiers'].values()]
                    ])
                    
                    # Check if any query matches
                    if any(query in searchable_text for query in queries):
                        info = {
                            'feature_id': feature_entry['locus_tag'],
                            'feature_name': feature_entry['gene_name'],
                            'feature_location': feature_entry['location'],
                            'feature_description': feature_entry['description']
                        }
                        results_list.append((feature_entry['record_id'], info))
                        
            return results_list
            
        except Exception as e:
            self.logger.error(f"Error in genbank_search: {str(e)}")
            raise

    def get_gene_data(self, gene_identifier):
        """Get gene data using the index for faster retrieval"""
        try:
            if not gene_identifier:
                return None
                
            # Ensure string conversion and proper formatting
            gene_identifier = str(gene_identifier).strip().lower()
            
            if hasattr(self, '_index'):
                # Try exact match first
                if gene_identifier in self._index['locus_tags']:
                    gene_info = self._index['locus_tags'][gene_identifier]
                    record_id = gene_info['record_id']
                    return {
                        'sequence': self._index['sequences'][record_id],
                        'info': gene_info
                    }
                    
                # Try case-insensitive match
                for key, value in self._index['locus_tags'].items():
                    if str(key).lower() == gene_identifier:
                        record_id = value['record_id']
                        return {
                            'sequence': self._index['sequences'][record_id],
                            'info': value
                        }
                        
                return None
                
        except Exception as e:
            self.logger.error(f"Error in get_gene_data: {str(e)}")
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

