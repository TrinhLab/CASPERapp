from models.CSPRparser import CSPRparser
from models.OffTargetModel import OffTargetModel
import os
import re
import traceback
from PyQt6.QtCore import QObject, pyqtSignal

class GenerateLibraryModel(QObject):
    progress_updated = pyqtSignal(int)  # Signal to emit progress updates
    
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.logger
        self.parser = None
        self.targets_data = {}
        self._deleted_targets = {}
        self.off_target_model = OffTargetModel(global_settings)
        
    def initialize_parser(self, cspr_file):
        """Initialize CSPR parser"""
        self.parser = CSPRparser(cspr_file, self.global_settings.get_casper_info_path())
        
    def get_organism_to_files(self):
        """Get mapping of organisms to their files from global settings"""
        return self.global_settings.get_organism_files()

    def generate_library(self, selected_targets, settings):
        """Generate library with given settings"""
        try:
            self.logger.debug(f"Generating library with settings: {settings}")
            
            # Process targets based on settings
            processed_targets = self._process_targets(
                selected_targets,
                settings['min_score'],
                settings['five_prime_seq'],
                settings['target_range_start'],
                settings['target_range_end']
            )
            
            if settings.get('find_off_targets'):
                # Write targets to temp file for off-target analysis
                self._write_targets_to_temp(processed_targets)
                
                # Get organism and endonuclease from home window
                if hasattr(self.global_settings, '_current_home_window'):
                    organism = self.global_settings._current_home_window.view.combo_box_organism.currentText()
                    endonuclease = self.global_settings._current_home_window.view.combo_box_endonuclease.currentText()
                else:
                    raise ValueError("Could not access home window to get organism and endonuclease")
                
                if not organism or not endonuclease:
                    raise ValueError("Could not determine organism or endonuclease from home window")
                
                self.logger.debug(f"Using organism: {organism} and endonuclease: {endonuclease} for off-target analysis")
                
                # Setup off-target parameters
                off_target_params = {
                    'organism': organism,
                    'endonuclease': endonuclease,
                    'max_mismatches': 4,  # Default value from old implementation
                    'tolerance': 0.05,  # Default value from old implementation
                    'average_output': True,
                    'save_output': False,
                    'output_filename': '',
                    'targets': selected_targets,
                    'annotation_file': self.global_settings.get_current_annotation_file()
                }
                
                # Connect to off-target model signals
                self.off_target_model.progress_updated.connect(self._handle_off_target_progress)
                self.off_target_model.results_ready.connect(lambda results: self._handle_off_target_results(results, processed_targets, settings))
                
                # Start off-target analysis
                self.off_target_model.start_analysis(off_target_params)
                return True
            else:
                # Generate output for each target
                output_data = self._generate_output(
                    processed_targets,
                    settings['guides_per_gene'],
                    settings['space_between_guides']
                )
                
                self.logger.debug(f"Output data: {output_data}")
                
                # Write output to file
                self._write_output(output_data, settings)
                return True
            
        except Exception as e:
            self.logger.error(f"Error generating library: {str(e)}")
            self.logger.error(traceback.format_exc())
            raise
            
    def _write_targets_to_temp(self, processed_targets):
        """Write targets to temp file for off-target analysis"""
        try:
            temp_path = os.path.join(self.global_settings.get_db_path(), 'temp.txt')
            
            with open(temp_path, 'w') as f:
                for gene in processed_targets:
                    for target in processed_targets[gene]:
                        # Format: position;sequence;pam;score;strand
                        entry = f"{target['position']};{target['sequence']};{target['pam']};{target['score']};{target['strand']}\n"
                        f.write(entry)
                        
            self.logger.debug(f"Wrote targets to temp file: {temp_path}")
            
        except Exception as e:
            self.logger.error(f"Error writing targets to temp file: {str(e)}")
            raise
            
    def _handle_off_target_progress(self, value, status):
        """Handle progress updates from off-target analysis"""
        self.progress_updated.emit(value)
        
    def _handle_off_target_results(self, results, processed_targets, settings):
        """Handle results from off-target analysis"""
        try:
            scores_dict, _ = results
            
            # Update targets with off-target scores
            for gene in processed_targets:
                for target in processed_targets[gene]:
                    if target['sequence'] in scores_dict:
                        target['off_target_score'] = scores_dict[target['sequence']]
            
            # Generate output with updated targets
            output_data = self._generate_output(
                processed_targets,
                settings['guides_per_gene'],
                settings['space_between_guides']
            )
            
            # Write output to file
            self._write_output(output_data, settings)
            
            # Clean up temp file
            temp_path = os.path.join(self.global_settings.get_db_path(), 'temp.txt')
            if os.path.exists(temp_path):
                os.remove(temp_path)
                
        except Exception as e:
            self.logger.error(f"Error handling off-target results: {str(e)}")
            raise
        
    def _process_targets(self, targets, min_score, five_prime_seq, start_range, end_range):
        """Process and filter targets based on criteria"""
        processed = {}
        self._deleted_targets = {}  # Store deleted targets for modify parameters option
        
        self.logger.debug(f"Processing {len(targets)} targets with min_score={min_score}")
        
        for target in targets:
            # Use feature_name (gene name) as key instead of feature_id
            gene_name = target.get('feature_name', target.get('feature_id'))
            
            if gene_name not in processed:
                processed[gene_name] = []
                self._deleted_targets[gene_name] = []
                
            # Create a copy of the target data
            target_data = target.copy()
            
            # Check if target passes filters
            if self._passes_filters(target_data, min_score, five_prime_seq, start_range, end_range):
                processed[gene_name].append(target_data)
            else:
                self._deleted_targets[gene_name].append(target_data)
                
        # Log first 5 targets for each gene for debugging
        for gene in processed:
            if processed[gene]:
                self.logger.debug(f"First 5 targets for gene {gene}: {[t['score'] for t in processed[gene][:5]]}")
                self.logger.debug(f"First 5 target positions for gene {gene}: {[t['position'] for t in processed[gene][:5]]}")
                
        return processed
        
    def _passes_filters(self, target, min_score, five_prime_seq, start_range, end_range):
        """Check if target passes all filters"""
        try:
            # Score filter - convert score to float and compare
            target_score = float(target.get('score', 0))
            if target_score < min_score:
                self.logger.debug(f"Target failed score filter: {target_score} < {min_score}, target: {target}")
                return False
            
            self.logger.debug(f"Target passed score filter: {target_score} >= {min_score}, target: {target}")
            
            # Poly-T filter (5-10 consecutive T's)
            if re.search("T{5,10}", target['sequence']):
                self.logger.debug(f"Target failed poly-T filter: {target['sequence']}")
                return False
            
            self.logger.debug(f"Target passed poly-T filter: {target['sequence']}")
            
            # 5' sequence filter
            if five_prime_seq and not target['sequence'].startswith(five_prime_seq.upper()):
                self.logger.debug(f"Target failed 5' sequence filter")
                return False
            
            # Range filter
            if start_range != 0 or end_range != 100:
                position_ratio = self._calculate_position_ratio(target)
                if not (start_range/100 <= position_ratio <= end_range/100):
                    self.logger.debug(f"Target failed range filter: {position_ratio}")
                    return False
                
            return True
            
        except Exception as e:
            self.logger.error(f"Error in _passes_filters: {str(e)}")
            self.logger.error(f"Target data: {target}")
            return False
        
    def _calculate_position_ratio(self, target):
        """Calculate relative position ratio in gene"""
        try:
            # Get gene length from feature info
            gene_length = abs(int(target.get('end', 0)) - int(target.get('start', 0)))
            if gene_length == 0:
                return 0
            
            # Calculate position relative to gene start
            target_pos = abs(int(target['position']))
            gene_start = int(target.get('start', 0))
            
            # If gene is on reverse strand, flip the position calculation
            if target.get('strand', '+') == '-':
                ratio = (target.get('end', 0) - target_pos) / gene_length
            else:
                ratio = (target_pos - gene_start) / gene_length
            
            self.logger.debug(f"Position ratio: {ratio} (pos={target_pos}, start={gene_start}, len={gene_length}, strand={target.get('strand', '+')}")
            return ratio
            
        except Exception as e:
            self.logger.error(f"Error calculating position ratio: {str(e)}")
            return 0
        
    def _generate_output(self, processed_targets, guides_per_gene, space_between):
        """Generate output with proper spacing between guides"""
        output = {}
        
        for gene_id, targets in processed_targets.items():
            output[gene_id] = []
            
            # First sort by score (descending)
            targets.sort(key=lambda x: float(x['score']), reverse=True)
            self.logger.debug(f"First 5 targets positions for gene {gene_id} by score: {[(t['position'], t['score']) for t in targets[:5]]}")
            
            # Then sort by position
            targets.sort(key=lambda x: abs(int(x['position'])))
            self.logger.debug(f"First 5 targets positions for gene {gene_id}: {[t['position'] for t in targets[:5]]}")
            
            i = 0  # Counter for selected guides
            vec_index = 0  # Index for current target being considered
            prev_target = None
            
            while i < guides_per_gene:
                if len(targets) == 0 or vec_index >= len(targets):
                    break
                
                current = targets[vec_index]
                
                # For first target, just add it
                if prev_target is None:
                    output[gene_id].append(current)
                    prev_target = current
                    i += 1
                else:
                    # Check spacing from previous target
                    distance = abs(int(current['position']) - int(prev_target['position']))
                    
                    if distance >= space_between:
                        # Look ahead for better scoring targets within this space
                        best_target = current
                        look_ahead_index = vec_index + 1
                        
                        while look_ahead_index < len(targets):
                            next_target = targets[look_ahead_index]
                            next_distance = abs(int(next_target['position']) - int(prev_target['position']))
                            
                            # If we've gone too far, break
                            if next_distance >= space_between:
                                break
                                
                            # If this target has better score
                            if float(next_target['score']) > float(best_target['score']):
                                best_target = next_target
                                
                            look_ahead_index += 1
                            
                        output[gene_id].append(best_target)
                        prev_target = best_target
                        i += 1
                        
                        # Move vec_index past the selected target's position
                        while vec_index < len(targets) and abs(int(targets[vec_index]['position'])) <= abs(int(best_target['position'])):
                            vec_index += 1
                        continue
                
                vec_index += 1
            
            # Sort final output by position
            output[gene_id].sort(key=lambda x: abs(int(x['position'])))
            
            # If gene is on negative strand, reverse the order
            if output[gene_id] and output[gene_id][0].get('strand', '+') == '-':
                output[gene_id].reverse()
                
            self.logger.debug(f"Selected targets positions for gene {gene_id}: {[t['position'] for t in output[gene_id]]}")
        
        return output
        
    def _write_output(self, output_data, settings):
        """Write library to output file"""
        output_file = settings['output_file']
        if not output_file.endswith('.csv'):
            output_file += '.csv'
            
        with open(output_file, 'w') as f:
            # Write header
            headers = ['Gene Name', 'Sequence', 'On-Target Score']
            if settings.get('find_off_targets'):
                headers.append('Off-Target Score')
            headers.extend(['Location', 'PAM', 'Strand'])
            f.write(','.join(headers) + '\n')
            
            # Write data
            for gene_id, targets in output_data.items():
                if not targets:  # Skip genes with no targets
                    continue
                    
                for i, target in enumerate(targets, 1):
                    # Use feature_name (gene name) instead of feature_id
                    gene_name = target.get('feature_name', gene_id)
                    
                    # Format gene name with index
                    tag_id = f"{gene_name}-{i}"
                    tag_id = tag_id.replace(',', '')  # Remove any commas
                    
                    # If target was modified by parameters, add asterisks
                    if settings.get('modify_params') and target.get('modified'):
                        tag_id = "**" + tag_id
                    
                    row = [
                        tag_id,
                        target['sequence'],
                        str(target['score'])
                    ]
                    
                    if settings.get('find_off_targets'):
                        row.append(str(target.get('off_target_score', '')))
                        
                    row.extend([
                        str(abs(int(target['position']))),
                        target['pam'],
                        target['strand'][0]  # Just take first character of strand
                    ])
                    
                    f.write(','.join(row) + '\n')
