import Bio
from PyQt6.QtWidgets import QMessageBox
from Bio import SeqIO
import os
import traceback
import pickle

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
            # Don't process if file_path is empty
            if not file_path:
                self.logger.warning("Empty annotation file path provided")
                self._index = {'locus_tags': {}}  # Initialize empty index
                return

            # Normalize path and remove any trailing slashes
            file_path = os.path.normpath(file_path)
            
            # Verify file exists and is a file (not a directory)
            if not os.path.isfile(file_path):
                self.logger.error(f"Invalid annotation file path: {file_path}")
                self._index = {'locus_tags': {}}
                return

            if self.annotation_file_name != file_path:
                self.annotation_file_name = file_path
                
                # Set index file path
                self.index_file = f"{file_path}.index"
                
                # Load or create index
                if not self._load_index():
                    self.logger.debug("Index not found or outdated, creating new index...")
                    self._create_index()
                
        except Exception as e:
            self.logger.error(f"Error in set_annotation_file: {str(e)}")
            raise

    def _create_index(self):
        try:
            self.logger.debug("Creating gene index file...")
            
            # Initialize optimized index structure - no sequences stored
            index_data = {
                'locus_tags': {},  # Only store essential data
            }
            
            # Process records
            record_count = 0
            feature_count = 0
            
            # Priority order for feature types (higher index = higher priority)
            feature_priority = {
                'CDS': 0,
                'gene': 1,
                'mRNA': 2,
                'tRNA': 2,
                'rRNA': 2,
                'ncRNA': 2
            }
            
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                record_count += 1
                
                # Process features
                for feature in record.features:
                    if feature.type in feature_priority:
                        feature_count += 1
                        
                        # Get essential feature info
                        locus_tag = None
                        if 'locus_tag' in feature.qualifiers:
                            locus_tag = feature.qualifiers['locus_tag'][0]
                        elif 'gene' in feature.qualifiers:
                            locus_tag = feature.qualifiers['gene'][0]
                    
                        # Only process features with valid locus tags
                        if locus_tag and locus_tag.lower() != "n/a":
                            # Get description, use product as fallback
                            description = feature.qualifiers.get('description', ['N/A'])[0]
                            if description == 'N/A' or not description:
                                description = feature.qualifiers.get('product', ['N/A'])[0]
                                if locus_tag == "BN896_RS00070":
                                    print(f"Feature description: {description}")

                            # Handle joined locations
                            if isinstance(feature.location, Bio.SeqFeature.CompoundLocation):
                                if locus_tag == "BN896_RS00070":
                                    print(f"Feature location: {feature.location}")
                                    # print(f"Feature product:" )
                                # Get all parts of the joined location
                                parts = feature.location.parts
                                # Find min start and max end across all parts
                                start = min(int(part.start) for part in parts)
                                end = max(int(part.end) for part in parts)
                                
                                # Format parts with strand info
                                formatted_parts = [
                                    f"{int(part.start)}..{int(part.end)}({'+' if part.strand == 1 else '-'})" 
                                    for part in parts
                                ]
                                
                                # If on minus strand, reverse the order of parts
                                if feature.location.strand == -1:
                                    formatted_parts.reverse()
                                    
                                # Join parts into full location string
                                full_location = ','.join(formatted_parts)
                                
                                # Get overall strand for location field
                                strand = '+' if feature.location.strand == 1 else '-'
                            else:
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                strand = '+' if feature.location.strand == 1 else '-'
                                full_location = f"{start}..{end}({strand})"
                            
                            # Get gene name, use description if gene name is N/A
                            gene_name = feature.qualifiers.get('gene', ['N/A'])[0]
                            if gene_name == 'N/A':
                                gene_name = description  # Use description as name if no gene name
                            
                            # Create new feature entry
                            feature_entry = {
                                'feature_type': feature.type,
                                'chromosome': record.id,
                                'location': f"{start}:{end}({strand})",
                                'full_location': full_location,
                                'gene_name': gene_name,
                                'description': description,
                                'start': start,
                                'end': end
                            }

                            # Update index based on modified priority logic
                            if locus_tag in index_data['locus_tags']:
                                existing_entry = index_data['locus_tags'][locus_tag]
                                existing_priority = feature_priority[existing_entry['feature_type']]
                                current_priority = feature_priority[feature.type]
                                
                                # Always create a merged entry
                                merged_entry = existing_entry.copy()
                                
                                # Update feature type only if priority is higher
                                if current_priority >= existing_priority:
                                    merged_entry['feature_type'] = feature.type
                                
                                # Always update description if new one is not N/A
                                if feature_entry['description'] != 'N/A':
                                    merged_entry['description'] = feature_entry['description']
                                    # If gene name is N/A, use the new description
                                    if merged_entry['gene_name'] == 'N/A':
                                        merged_entry['gene_name'] = feature_entry['description']
                                
                                # Update other fields if they're not 'N/A'
                                if feature_entry['gene_name'] != 'N/A':
                                    merged_entry['gene_name'] = feature_entry['gene_name']
                                
                                # Always update location information if priority is higher
                                if current_priority >= existing_priority:
                                    merged_entry.update({
                                        'location': feature_entry['location'],
                                        'full_location': feature_entry['full_location'],
                                        'start': feature_entry['start'],
                                        'end': feature_entry['end']
                                    })
                                
                                index_data['locus_tags'][locus_tag] = merged_entry
                            else:
                                # New entry
                                index_data['locus_tags'][locus_tag] = feature_entry
            
            # Save compressed index to file
            with open(self.index_file, 'wb') as f:
                pickle.dump(index_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            self._index = index_data

            # if locus tag is CAALFM_C304810CA
            if 'CAALFM_C304810CA' in index_data['locus_tags']:
                print(f"Locus tag CAALFM_C304810CA found: {index_data['locus_tags']['CAALFM_C304810CA']}")

            self.logger.debug(f"Index creation complete. Records: {record_count}, Features: {feature_count}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating index: {str(e)}")
            return False

    def _load_index(self):
        try:
            self.logger.debug(f"Loading index from: {self.index_file}")
            if not os.path.exists(self.index_file):
                return False

            # Check if index is older than GenBank file
            if os.path.getmtime(self.index_file) < os.path.getmtime(self.annotation_file_name):
                return False

            with open(self.index_file, 'rb') as f:
                self._index = pickle.load(f)
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
            
            queries = {q.lower() for q in queries}
            
            # Search through index
            if hasattr(self, '_index') and 'locus_tags' in self._index:
                for locus_tag, feature_entry in self._index['locus_tags'].items():
                    # Create searchable text including feature type
                    searchable_text = ' '.join([
                        feature_entry.get('gene_name', '').lower(),
                        locus_tag.lower(),
                        feature_entry.get('description', '').lower(),
                        feature_entry.get('feature_type', '').lower()  # Add feature type to searchable text
                    ])
                    
                    # Check if any query matches
                    if any(query in searchable_text for query in queries):
                        info = {
                            'feature_id': locus_tag,
                            'feature_name': feature_entry.get('gene_name', 'N/A'),
                            'feature_full_location': feature_entry.get('full_location', 'N/A'),
                            'feature_location': feature_entry.get('location', 'N/A'),
                            'feature_description': feature_entry.get('description', 'N/A'),
                            'feature_type': feature_entry.get('feature_type', 'CDS')
                        }
                        results_list.append((feature_entry.get('chromosome', ''), info))
                        
            return results_list
            
        except Exception as e:
            self.logger.error(f"Error in genbank_search: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
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
                        'full_location': gene_info.get('full_location', ''),  # Add full location
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
                            'full_location': value.get('full_location', ''),  # Add full location
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
            # Parse the GenBank file and find the right record
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                if record.id == gene_info['chromosome']:  # Use full chromosome name
                    sequence = str(record.seq)
                    
                    # Get sequence with padding
                    padding = 30
                    start = max(0, gene_info['start'] - padding)
                    end = min(len(sequence), gene_info['end'] + padding)
                    padded_sequence = sequence[start:end]
                    return padded_sequence
            return None
            
        except Exception as e:
            self.logger.error(f"Error getting sequence for gene: {str(e)}")
            return None

    def _get_sequence_for_position(self, chrom, start, end):
        """Get sequence for a specific position range from the GenBank file
        
        Args:
            chrom (str): Chromosome identifier
            start (int): Start position (0-based)
            end (int): End position
            
        Returns:
            str: The sequence for the specified range with padding, or None if not found
        """
        try:
            self.logger.debug(f"Getting sequence for position {chrom}:{start}-{end}")
            # Parse the GenBank file and find the right record
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                if record.id == chrom:  # Use full chromosome name
                    sequence = str(record.seq)
                    
                    # Get sequence with padding
                    padding = 30
                    padded_start = max(0, start - padding)
                    padded_end = min(len(sequence), end + padding)
                    padded_sequence = sequence[padded_start:padded_end]
                    return padded_sequence
            return None
            
        except Exception as e:
            self.logger.error(f"Error getting sequence for position: {str(e)}")
            return None