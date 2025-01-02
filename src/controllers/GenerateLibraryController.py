from views.GenerateLibraryView import GenerateLibraryView
from models.GenerateLibraryModel import GenerateLibraryModel
from PyQt6.QtCore import QObject
import os

class GenerateLibraryController(QObject):
    def __init__(self, global_settings, selected_targets=None):
        # Initialize QObject first
        super(GenerateLibraryController, self).__init__()
        
        self.global_settings = global_settings
        self.logger = global_settings.logger
        self.model = GenerateLibraryModel(global_settings)
        self.view = GenerateLibraryView(global_settings)
        
        # Get selected targets from global settings if not provided
        if selected_targets is None and hasattr(self.global_settings, '_current_selected_targets'):
            selected_targets = self.global_settings._current_selected_targets
        
        self.selected_targets = selected_targets or []
        
        self.view.ledFileName.setText('eck_12_spCas9_lib')
        
        # Get CSPR file path
        if self.selected_targets and len(self.selected_targets) > 0:
            # Get organism name from the first target's chromosome
            chrom = self.selected_targets[0].get('full_chromosome', '')
            if chrom:
                # Extract organism name from chromosome ID
                org_name = chrom.split('.')[0]
                default_filename = f"{org_name}_lib.csv"
                self.view.ledFileName.setText(default_filename)
                
                # Get CSPR file path and initialize parser
                org_files = self.model.get_organism_to_files()
                endonuclease = self.selected_targets[0].get('endonuclease', '').lower()
                if org_name in org_files and endonuclease in org_files[org_name]:
                    cspr_file = os.path.join(
                        self.global_settings.get_db_path(),
                        org_files[org_name][endonuclease][0]
                    )
                    self.model.initialize_parser(cspr_file)
                    
                    # Get guide data for each target
                    processed_targets = []
                    for target in self.selected_targets:
                        target_info = [{
                            'feature_id': target['feature_id'],
                            'feature_name': target['feature_name'],
                            'start': target['start'],
                            'end': target['end'],
                            'chromosome': target['chromosome']
                        }]
                        guides = self.model.parser.read_targets_batch(
                            target['chromosome'],
                            target_info,
                            target['endonuclease']
                        )
                        if guides:
                            for guide in guides:
                                guide.update({
                                    'feature_id': target['feature_id'],
                                    'feature_name': target['feature_name'],
                                    'feature_type': target['feature_type'],
                                    'feature_description': target['feature_description'],
                                    'start': target['start'],
                                    'end': target['end']
                                })
                            processed_targets.extend(guides)
                    
                    self.selected_targets = processed_targets
                else:
                    self.selected_targets = selected_targets
        
        self._connect_signals()
        
    def _connect_signals(self):
        """Connect view signals to controller methods"""
        try:
            self.view.submit_clicked.connect(self._handle_submit)
        except Exception as e:
            self.logger.error(f"Error connecting signals: {str(e)}")
        
    def show(self):
        """Show the generate library window"""
        try:
            # Store reference to prevent garbage collection
            self.global_settings.main_window._current_generate_library_controller = self
            
            # Show the view
            self.view.show()
            
        except Exception as e:
            self.logger.error(f"Error showing generate library window: {str(e)}")
        
    def _handle_submit(self, settings):
        """Handle submit button click"""
        try:
            self.logger.debug(f"Handling submit with settings: {settings}")
            
            # Validate settings
            self._validate_settings(settings)
            
            # Get guide data for targets if not already processed
            if not hasattr(self, 'processed_targets'):
                self.processed_targets = []
                
                # Get CSPR file path
                if self.selected_targets and len(self.selected_targets) > 0:
                    first_target = self.selected_targets[0]
                    self.logger.debug(f"First target: {first_target}")
                    
                    # Get organism name from the main window's home window view
                    home_window = self.global_settings._current_home_window
                    org_name = home_window.view.combo_box_organism.currentText()
                    endonuclease = first_target['endonuclease'].lower()
                    
                    self.logger.debug(f"Looking for CSPR file for {org_name} and {endonuclease}")
                    
                    # Get CSPR file path
                    org_files = self.model.get_organism_to_files()
                    if org_name in org_files:
                        # Debug available endonucleases
                        self.logger.debug(f"Available endonucleases for {org_name}: {list(org_files[org_name].keys())}")
                        
                        # Try both lowercase and original case
                        cspr_file = None
                        if endonuclease in org_files[org_name]:
                            cspr_file = org_files[org_name][endonuclease][0]
                        elif first_target['endonuclease'] in org_files[org_name]:
                            cspr_file = org_files[org_name][first_target['endonuclease']][0]
                        
                        if cspr_file:
                            cspr_path = os.path.join(self.global_settings.get_db_path(), cspr_file)
                            self.logger.debug(f"Using CSPR file: {cspr_path}")
                            
                            self.model.initialize_parser(cspr_path)
                            
                            # Process each target to get guide data
                            for target in self.selected_targets:
                                # Get start and end from location if not present
                                if 'start' not in target or 'end' not in target:
                                    if 'location' in target:
                                        start, end = map(int, target['location'].split('-'))
                                        target['start'] = start
                                        target['end'] = end
                                
                                target_info = [{
                                    'feature_id': target['feature_id'],
                                    'feature_name': target['feature_name'],
                                    'start': int(target['start']),
                                    'end': int(target['end']),
                                    'chromosome': target['chromosome']
                                }]
                                
                                self.logger.debug(f"Searching guides for target: {target_info}")
                                
                                # Get guides for this target
                                guides = self.model.parser.read_targets_batch(
                                    target['chromosome'],
                                    target_info,
                                    target['endonuclease']
                                )
                                
                                if guides:
                                    self.logger.debug(f"Found {len(guides)} guides for target {target['feature_id']}")
                                    # Add target info to each guide
                                    for guide in guides:
                                        guide.update({
                                            'feature_id': target['feature_id'],
                                            'feature_name': target['feature_name'],
                                            'feature_type': target.get('feature_type', 'CDS'),
                                            'feature_description': target.get('feature_description', ''),
                                            'start': int(target['start']),
                                            'end': int(target['end'])
                                        })
                                    self.processed_targets.extend(guides)
                                else:
                                    self.logger.warning(f"No guides found for target {target['feature_id']}")
                                    
                            self.logger.debug(f"Processed {len(self.processed_targets)} total guides from CSPR file")
                        else:
                            raise ValueError(f"Could not find CSPR file for endonuclease {endonuclease}")
                    else:
                        raise ValueError(f"Could not find organism {org_name} in database")
                
                # Connect to model's progress signal if off-target analysis is enabled
                if settings.get('find_off_targets'):
                    self.model.progress_updated.connect(self._handle_progress)
                
                # Generate library using processed targets
                success = self.model.generate_library(
                    self.processed_targets if hasattr(self, 'processed_targets') else self.selected_targets,
                    settings
                )
                
                if success and not settings.get('find_off_targets'):
                    self.view.show_success("Library generated successfully!")
                    self.view.close()
                    
        except ValueError as e:
            self.logger.error(f"Validation error: {str(e)}")
            self.view.show_error("Invalid Input", str(e))
        except Exception as e:
            self.logger.error(f"Error generating library: {str(e)}")
            self.view.show_error(
                "Error",
                f"An error occurred while generating the library: {str(e)}"
            )
            
    def _validate_settings(self, settings):
        """Validate library generation settings"""
        try:
            if not settings['output_file']:
                raise ValueError("Please specify an output file")
                
            if settings['target_range_start'] >= settings['target_range_end']:
                raise ValueError("Start range must be less than end range")
                
            if settings['target_range_start'] < 0 or settings['target_range_end'] > 100:
                raise ValueError("Target range must be between 0 and 100")
                
            if settings['space_between_guides'] < 0:
                raise ValueError("Space between guides must be positive")
                
            if settings.get('find_off_targets'):
                max_score = settings.get('max_off_target_score')
                if max_score is None or not 0 <= max_score <= 0.5:
                    raise ValueError("Maximum off-target score must be between 0 and 0.5 inclusive")
                    
        except Exception as e:
            self.logger.error(f"Settings validation error: {str(e)}")
            raise
            
    def _handle_progress(self, value):
        """Handle progress updates from model"""
        try:
            self.view.progBar.setValue(value)
        except Exception as e:
            self.logger.error(f"Error updating progress: {str(e)}")
