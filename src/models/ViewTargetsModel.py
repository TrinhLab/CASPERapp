from models.CSPRparser import CSPRparser
from models.BaseModel import BaseModel
from Bio import SeqIO
from collections import defaultdict
import traceback
import os

class ViewTargetsModel(BaseModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        
        # Initialize model state
        self.guides = []
        self.cspr_parser = None
        self.gene_sequence = ""
        self.highlighted_sequence = ""
        self.gene_info = {}
        self.available_genes = []
        self.filter_options = {}
        self.scoring_options = {}
        self.current_gene_start = 0
        self.current_gene_end = 0
        self.extended_sequence = ""
        self.chromosome = ""
        
        # Initialize caches
        self._gene_data_cache = {}
        self._sequence_cache = {}
        self._parser_cache = {}
        self._chromosome_seqs = {}
        self._cached_guides = {}

    def _clear_caches(self):
        """Clear all model-specific caches"""
        self._gene_data_cache.clear()
        self._sequence_cache.clear()
        self._parser_cache.clear()
        self._chromosome_seqs.clear()
        self._cached_guides.clear()
        
        # Clear other stored data
        self.gene_sequence = ""
        self.highlighted_sequence = ""
        self.gene_info = {}
        self.available_genes = []

    def _ensure_annotation_parser(self) -> bool:
        """Ensure annotation parser is initialized
        
        Returns:
            bool: True if parser is ready, False otherwise
        """
        if self.annotation_parser is None:
            try:
                self._initialize_annotation_parser()
                return True
            except Exception as e:
                self.logger.error(f"Failed to initialize annotation parser: {str(e)}")
                return False
        return True

    def load_guides(self, selected_targets, organism, endonuclease):
        """Load guides with proper error handling"""
        try:
            self.organism = organism
            self.endonuclease = endonuclease

            # Get CSPR parser from cache or create new one
            cspr_key = f"{organism}_{endonuclease}"
            if cspr_key in self._parser_cache:
                self.cspr_parser = self._parser_cache[cspr_key]
                self.logger.debug("Using cached CSPR parser")
            else:
                org_files = self.global_settings.get_organism_files()
                if organism not in org_files or endonuclease not in org_files[organism]:
                    self.logger.error(f"No CSPR file found for {organism} and {endonuclease}")
                    return

                cspr_file = org_files[organism][endonuclease][0]
                cspr_path = os.path.join(self.global_settings.get_db_path(), cspr_file)
                self.cspr_parser = CSPRparser(cspr_path, self.global_settings.get_casper_info_path())
                self._parser_cache[cspr_key] = self.cspr_parser

            # Initialize guides and genes
            self.guides = []
            self.available_genes = set()  # Clear existing genes
            
            # Use a set to track unique guide positions
            seen_guides = set()
            
            batch_guides = defaultdict(list)
            for target in selected_targets:
                # Use full_chromosome directly if available, otherwise use chromosome
                chrom = target.get('full_chromosome', target['chromosome'])
                
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
                    # Add only the feature_id from the original target
                    self.available_genes.add(target['feature_id'])

            # Process guides by chromosome
            unique_guides = {}  # Use dict to track unique guides by sequence
            for chrom, guides in batch_guides.items():
                # Use full chromosome ID for CSPR lookup
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

    def get_gene_data(self, locus_tag):
        """Get gene data with proper error handling"""
        try:
            if not locus_tag:
                self.logger.debug("No locus tag provided")
                return None
                
            # Check model cache first
            if locus_tag in self._gene_data_cache:
                return self._gene_data_cache[locus_tag]
            
            # Ensure parser is initialized
            if not self._ensure_annotation_parser():
                return None
            
            # Get gene data from parser with proper string conversion
            gene_data = None
            if isinstance(locus_tag, (str, int)):
                locus_tag_str = str(locus_tag).strip()
                # Look up by locus tag directly
                gene_data = self.annotation_parser.get_gene_data(locus_tag_str.lower())
                
            if gene_data:
                self._gene_data_cache[locus_tag] = gene_data
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

    def get_gene_sequence(self, identifier):
        """Get gene sequence with optimized caching and minimal I/O"""
        try:
            self.logger.debug(f"Getting gene sequence for identifier: {identifier}")
            self.logger.debug(f"View exons only is: {getattr(self, '_view_exons_only', False)}")
            
            # Regular gene-based search
            gene_data = self.get_gene_data(identifier)
            if not gene_data or 'info' not in gene_data:
                self.logger.warning(f"No gene data found for locus tag: {identifier}")
                return None
            
            # Check if we're in exons-only mode
            if getattr(self, '_view_exons_only', False):
                print(f"gene_data: {gene_data}")
                full_location = gene_data['info'].get('full_location', '')
                print(f"Full location: {full_location}")
                if full_location:
                    if ',' in full_location:  # Multiple exons
                        self.logger.debug(f"Processing exons from full location: {full_location}")
                        exon_sequences = []
                        full_sequence = gene_data['sequence']  # Use sequence from gene_data
                        
                        # Calculate padding offset
                        padding = 30
                        gene_start = gene_data['info']['start']
                        padded_start = max(0, gene_start - padding)
                        padding_offset = gene_start - padded_start
                        
                        print(f"gene_start: {gene_start}, padded_start: {padded_start}, padding_offset: {padding_offset}")
                        
                        # Process each exon location
                        for exon in full_location.split(','):
                            # Extract coordinates and strand
                            coords = exon.split('(')[0]  # Get part before strand
                            strand = exon.split('(')[1][0]  # Get + or - from (+ or (-
                            start, end = map(int, coords.split('..'))

                            print(f"coords: {coords}, strand: {strand}, start: {start}, end: {end}")
                            
                            # Adjust coordinates relative to gene start and account for padding
                            relative_start = start - gene_start + padding_offset
                            relative_end = end - gene_start + padding_offset
                            print(f"relative_start: {relative_start}, relative_end: {relative_end}")
                            
                            # Get exon sequence from the padded sequence
                            exon_seq = full_sequence[relative_start:relative_end]
                            
                            exon_sequences.append(exon_seq)
                        
                        # Join exon sequences
                        sequence = ''.join(exon_sequences)
                        self.logger.debug(f"Created concatenated exon sequence of length: {len(sequence)}")
                        
                        return {
                            'sequence': sequence,
                            'info': gene_data['info'],
                            'start': gene_data['info']['start'],
                            'end': gene_data['info']['end']
                        }
                    else:  # Single location/exon
                        # Return sequence without padding for single exon
                        sequence = gene_data['sequence']
                        gene_start = gene_data['info']['start']
                        gene_end = gene_data['info']['end']
                        
                        # Calculate padding offset
                        padding = 30
                        padded_start = max(0, gene_start - padding)
                        padding_offset = gene_start - padded_start
                        
                        # Get sequence without padding
                        relative_start = padding_offset
                        relative_end = len(sequence) - padding_offset
                        sequence = sequence[relative_start:relative_end]
                        
                        return {
                            'sequence': sequence,
                            'info': gene_data['info'],
                            'start': gene_data['info']['start'],
                            'end': gene_data['info']['end']
                        }
            
            # If not in exons-only mode or no exons to process, return normal sequence
            if 'sequence' in gene_data:
                sequence = gene_data['sequence']
                
                # Format sequence with padding in lowercase (only if not in exons-only mode)
                if not hasattr(self, '_view_exons_only') or not self._view_exons_only:
                    padding = 30
                    start = gene_data['info']['start']
                    end = gene_data['info']['end']
                    padded_start = max(0, start - padding)
                    padded_end = min(len(sequence), end + padding)
                    
                    # Split sequence into parts
                    five_prime_pad = sequence[:start - padded_start].lower() if start > padded_start else ""
                    main_sequence = sequence[start - padded_start:end - padded_start].upper()
                    three_prime_pad = sequence[end - padded_start:].lower()
                    
                    # Combine parts
                    formatted_sequence = five_prime_pad + main_sequence + three_prime_pad
                else:
                    formatted_sequence = sequence
                
                result = {
                    'sequence': formatted_sequence,
                    'info': gene_data['info'],
                    'start': gene_data['info']['start'],
                    'end': gene_data['info']['end'],
                    'full_location': gene_data['info'].get('full_location', '')
                }
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
            if not self._ensure_annotation_parser():
                return None
                
            feature_info = {
                'chromosome': chrom,  # Use raw chromosome ID directly
                'start': start,  # Keep as is since annotation parser handles 0-based conversion
                'end': end
            }
            
            self.logger.debug(f"Getting sequence for feature info: {feature_info}")
            
            # Use annotation parser's method directly
            sequence = self.annotation_parser._get_sequence_for_position(chrom, start, end)
            if sequence:
                padding = 30
                
                # Handle start position padding
                if start == 0:  # Already 0-based for sequence operations
                    # No padding at start if starting at position 0
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
            
            chrom = gene_data['info']['chromosome']  
            
            self.logger.debug(f"Getting sequence for chromosome: {chrom}, start: {start}, end: {end}")
            
            # Use _get_sequence_for_position to get sequence with padding
            sequence = self._get_sequence_for_position(chrom, start, end)
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

    def set_view_exons_only(self, enabled):
        """Set whether to view exons only"""
        try:
            self.logger.debug(f"Setting view exons only to: {enabled}")
            self._view_exons_only = enabled
            # Clear cache when changing view mode
            self._sequence_cache.clear()
            self.logger.debug("Cleared sequence cache")
        except Exception as e:
            self.logger.error(f"Error setting view exons only: {str(e)}")

    def get_features_for_gene(self, locus_tag):
        """Get features for a specific gene"""
        try:
            if not self.annotation_parser:
                self._initialize_annotation_parser()
            
            features = []
            gene_data = self.get_gene_data(locus_tag)
            
            if gene_data and 'info' in gene_data:
                info = gene_data['info']
                
                # Add the main gene feature
                features.append({
                    'type': info['feature_type'],
                    'start': info['start'],
                    'end': info['end'],
                    'name': info['gene_name'],
                    'id': locus_tag,
                    'strand': '+' if '(+)' in info['location'] else '-'
                })
                
                # Parse additional features from full location if available
                if 'full_location' in info and ',' in info['full_location']:
                    for i, part in enumerate(info['full_location'].split(',')):
                        coords = part.split('(')[0]
                        strand = part.split('(')[1][0]
                        start, end = map(int, coords.split('..'))
                        
                        features.append({
                            'type': 'exon',
                            'start': start,
                            'end': end,
                            'name': f'Exon {i+1}',
                            'id': f'{locus_tag}_exon_{i+1}',
                            'strand': strand
                        })
                        
            return features
            
        except Exception as e:
            self.logger.error(f"Error getting features for gene: {str(e)}")
            return []

    def get_features_for_region(self, chromosome, start, end):
        """Get features within a specific region"""
        try:
            if not self.annotation_parser:
                self._initialize_annotation_parser()
            
            features = []
            
            # Search through index for features in this region
            if hasattr(self, '_index') and 'locus_tags' in self._index:
                for locus_tag, feature_info in self._index['locus_tags'].items():
                    if (feature_info['chromosome'] == chromosome and
                        feature_info['start'] <= end and
                        feature_info['end'] >= start):
                        
                        features.append({
                            'type': feature_info['feature_type'],
                            'start': feature_info['start'],
                            'end': feature_info['end'],
                            'name': feature_info['gene_name'],
                            'id': locus_tag,
                            'strand': '+' if '(+)' in feature_info['location'] else '-'
                        })
                        
            return features
            
        except Exception as e:
            self.logger.error(f"Error getting features for region: {str(e)}")
            return []