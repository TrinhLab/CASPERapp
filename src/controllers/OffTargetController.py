from models.OffTargetModel import OffTargetModel
from utils.ui import show_error
from views.OffTargetView import OffTargetView
from PyQt6.QtCore import QObject, pyqtSlot, pyqtSignal
import os

class OffTargetController(QObject):
    # Update signal to emit tuple of (scores, details)
    off_target_results_ready = pyqtSignal(tuple)  # Emits (scores_dict, details_dict)

    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()

        try:
            self.model = OffTargetModel(global_settings)
            self.view = OffTargetView(global_settings)
            
            self._init_ui()
            
            self._setup_connections()
            
            self.global_settings.theme_changed.connect(self.on_theme_changed)
        except Exception as e:
            show_error(self.global_settings, "Error initializing OffTargetController", str(e))

    def _init_ui(self):
        try:
            organisms = self.model.get_organisms()
            self.view.set_combo_box_organism(organisms)
            
            if organisms:
                endonucleases = self.model.get_endonucleases(organisms[0])
                self.view.set_combo_box_endonuclease(endonucleases)
            
            self.view.set_combo_box_max_mismatches()
            
        except Exception as e:
            self.logger.error(f"Error setting up initial data: {str(e)}")
            self.view.show_error("Setup Error", str(e))

    def _setup_connections(self):
        try:
            self.view.push_button_submit.clicked.connect(self._on_submit_clicked)
            
            # Connect to model's signals
            self.model.results_ready.connect(lambda results: self._on_results_ready(results))
            self.model.progress_updated.connect(self._on_progress_updated)
            
        except Exception as e:
            self.logger.error(f"Error connecting signals: {str(e)}")

    @pyqtSlot(str)
    def _on_organism_changed(self, organism):
        """Handle organism selection change"""
        try:
            endonucleases = self.model.get_endonucleases(organism)
            self.view.set_endonucleases(endonucleases)
        except Exception as e:
            self.logger.error(f"Error updating endonucleases: {str(e)}")

    @pyqtSlot()
    def _on_submit_clicked(self):
        """Handle submit button click"""
        try:
            # Get analysis parameters
            parameters = self.view.get_analysis_parameters()
            
            # Add stored targets
            if hasattr(self, '_targets'):
                parameters['targets'] = self._targets
            else:
                raise ValueError("No targets available for analysis")
            
            # Validate parameters
            if not self._validate_parameters(parameters):
                return
            
            self.model._write_targets_to_temp(parameters['targets'])
            
            # Start analysis
            self.model.start_analysis(parameters)
            
            # Update UI
            self.view.update_progress_bar(1, "Starting analysis...")
            self.view.push_button_submit.setEnabled(False)
                
        except Exception as e:
            self.logger.error(f"Error in submit action: {str(e)}")
            self.view.show_error("Submit Error", str(e))

    def _validate_parameters(self, parameters):
        """Validate analysis parameters"""
        try:
            if parameters['save_output']:
                if not parameters['output_filename']:
                    self.view.show_warning(
                        "Input Error",
                        "Please enter a valid output file name."
                    )
                    return False
                    
                output_path = os.path.join(
                    self.global_settings.get_db_path(),
                    parameters['output_filename']
                )
                if os.path.exists(output_path):
                    self.view.show_warning(
                        "File Error",
                        "Output file already exists. Please choose a new name."
                    )
                    return False
                    
            return True
            
        except Exception as e:
            self.logger.error(f"Error validating parameters: {str(e)}")
            return False

    @pyqtSlot()
    def _on_cancel_clicked(self):
        """Handle cancel button click"""
        try:
            self.model.stop_analysis()
            self._cleanup_temp_files()
            self.view.close()
        except Exception as e:
            self.logger.error(f"Error canceling analysis: {str(e)}")

    def _cleanup_temp_files(self):
        """Clean up temporary files"""
        try:
            # Get path to temp file
            temp_path = os.path.join(
                self.global_settings.get_off_target_dir_path(),
                'temp.txt'
            )
            
            # Clean up temp file
            if os.path.exists(temp_path):
                os.remove(temp_path)
                self.logger.debug(f"Removed temp file: {temp_path}")
                
            # Clean up local output file
            local_output = os.path.join(
                self.global_settings.get_off_target_dir_path(),
                'local_output.txt'
            )
            if os.path.exists(local_output):
                os.remove(local_output)
                self.logger.debug(f"Removed local output file: {local_output}")
                
        except Exception as e:
            self.logger.error(f"Error cleaning up temp files: {str(e)}")

    def _on_results_ready(self, results):
        """Handle results from model"""
        try:
            scores, details = results  # Unpack the tuple of results
            if scores:  # Make sure we have valid results
                # Store results
                self._off_target_results = scores
                self._off_target_details = details
                
                # Emit results signal immediately with both scores and details
                self.off_target_results_ready.emit((scores, details))
                
                # Close the window
                self.view.close()
                
                self.logger.debug(f"Received and emitted {len(scores)} off-target results with {len(details)} detailed results")
                
            else:
                self.logger.warning("Received empty results")
                self.view.show_warning(
                    "No Results",
                    "No off-target analysis results were generated."
                )
            
        except Exception as e:
            self.logger.error(f"Error handling results: {str(e)}")
            self.view.show_error("Results Error", str(e))

    def show(self):
        """Show the view and bring to front"""
        # Reset view state before showing
        self.view.prog_bar.setValue(0)
        self.view.push_button_submit.setEnabled(True)
        
        # Show and bring window to front
        self.view.show()
        self.view.raise_()  # Bring window to front
        self.view.activateWindow()  # Give window focus
        self.view.apply_theme()

    def initialize_analysis(self, parameters):
        """Initialize analysis with parameters from ViewTargets"""
        try:
            # Reset view state
            self.view.prog_bar.setValue(0)
            self.view.push_button_submit.setEnabled(True)
            
            # Reset model state
            if hasattr(self.model, '_current_parameters'):
                delattr(self.model, '_current_parameters')
            if hasattr(self, '_off_target_results'):
                delattr(self, '_off_target_results')
            
            # Set organism in view
            organism_index = self.view.combo_box_organism.findText(parameters['organism'])
            if organism_index >= 0:
                self.view.combo_box_organism.setCurrentIndex(organism_index)
                
            # Set endonuclease in view
            if 'endonuclease' in parameters:
                endo_index = self.view.combo_box_endonuclease.findText(parameters['endonuclease'])
                if endo_index >= 0:
                    self.view.combo_box_endonuclease.setCurrentIndex(endo_index)
            
            # Validate and set annotation file
            if 'annotation_file' not in parameters:
                raise ValueError("No annotation file provided in parameters")
                
            annotation_file = parameters['annotation_file']
            if not annotation_file:
                raise ValueError("Empty annotation file path provided")
                
            # Set annotation file in global settings
            self.global_settings.set_current_annotation_file(annotation_file)
            self.logger.debug(f"Set annotation file to: {annotation_file}")
            
            # Store targets for analysis
            if 'guides' in parameters:
                self._targets = parameters['guides']
            else:
                raise ValueError("No guides provided for analysis")
            
            # Reset radio buttons and fields
            self.view.radio_button_average_output_no.setChecked(True)
            self.view.radio_button_save_output_no.setChecked(True)
            self.view.line_edit_output_file.clear()
            self.view.line_edit_output_file.setEnabled(False)
            
            # Reset combo box for mismatches
            self.view.set_combo_box_max_mismatches()
            
            # Enable submit button
            self.view.push_button_submit.setEnabled(True)
            
        except Exception as e:
            self.logger.error(f"Error initializing analysis: {str(e)}")
            self.view.show_error("Initialization Error", str(e))

    @pyqtSlot(str)
    def on_theme_changed(self, theme):
        """Handle theme change"""
        self.view.apply_theme()

    def _on_progress_updated(self, value, status):
        """Handle progress updates"""
        try:
            self.view.update_progress_bar(value, status)
        except Exception as e:
            self.logger.error(f"Error updating progress: {str(e)}")
