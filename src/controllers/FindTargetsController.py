from models.FindTargetsModel import FindTargetsModel
from views.FindTargetsView import FindTargetsView
from PyQt6.QtWidgets import QMessageBox
from PyQt6.QtCore import QTimer
import time

class FindTargetsController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = FindTargetsModel(self.global_settings)
        self.view = FindTargetsView(self.global_settings)
        self.organism = None
        self.endonuclease = None
        self._input_data = None
        self._current_annotation_file = None
        
        # Connect to annotation file changes
        self.global_settings.annotation_file_changed.connect(self._on_annotation_file_changed)
        self._connect_signals()

    def _on_annotation_file_changed(self, new_annotation_file):
        """Handle annotation file changes by clearing and updating results"""
        try:
            self.global_settings.logger.debug(f"FindTargetsController received new annotation file: {new_annotation_file}")
            self._current_annotation_file = new_annotation_file
            
            # Clear the current results
            if self.view and hasattr(self.view, 'results_table'):
                self.view.clear_results()
                
                # If we have previous input data, rerun the search with the new annotation file
                if self._input_data:
                    self._input_data['annotation_file'] = new_annotation_file
                    self._process_input_data(self._input_data)
                
        except Exception as e:
            self.global_settings.logger.error(f"Error handling annotation file change: {str(e)}")

    def _connect_signals(self):
        """Connect view signals"""
        if self.view:
            self.view.push_button_view_targets.clicked.connect(self.view_targets)

    def find_targets(self, input_data):
        """Process input data and update existing view or create new one"""
        try:
            start_time = time.time()
            
            # Get current annotation file
            current_annotation = self.global_settings.get_current_annotation_file()
            input_data['annotation_file'] = current_annotation
            self._current_annotation_file = current_annotation
            self._input_data = input_data.copy()  # Store a copy of the input data
            
            # Process data and update view
            self._process_input_data(input_data)
            
            # If there's no existing tab, create one
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Find Targets")
            if not existing_tab:
                main_window.open_new_tab("Find Targets", self)
            
            total_time = time.time() - start_time
            self.global_settings.logger.debug(f"Total time to process find targets: {total_time:.2f} seconds")
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in find_targets: {str(e)}")
            raise

    def _process_input_data(self, input_data):
        """Process input data and update view"""
        try:
            start_time = time.time()
            
            self.global_settings.logger.debug(f"FindTargetsController processing input data: {input_data}")
            self.organism = input_data['organism']
            self.endonuclease = input_data['endonuclease']
            
            # Get new results
            search_start = time.time()
            results = self.model.find_targets(input_data)
            search_time = time.time() - search_start
            self.global_settings.logger.debug(f"Time to search: {search_time:.2f} seconds")
            self.global_settings.logger.debug(f"Found {len(results) if results else 0} targets")
            
            # Update view with new results
            view_start = time.time()
            if results:
                self.view.display_results(results)
            view_time = time.time() - view_start
            self.global_settings.logger.debug(f"Time to update view: {view_time:.2f} seconds")
            
            total_time = time.time() - start_time
            self.global_settings.logger.debug(f"Total time to process data: {total_time:.2f} seconds")
            
        except Exception as e:
            self.global_settings.logger.error(f"Error processing input data: {str(e)}")
            if self.view:
                QMessageBox.critical(self.view, "Error", f"An error occurred while processing data: {str(e)}")

    def view_targets(self):
        try:
            if not self.view:
                return
                
            selected_targets = self.view.get_selected_targets()
            print(f"Selected targets: {selected_targets}")
            print(f"Organism: {self.organism}")
            print(f"Endonuclease: {self.endonuclease}")
            if not selected_targets:
                QMessageBox.warning(self.view, "No Selection", "Please select targets to view.")
                return
            
            view_targets_controller = self.global_settings.get_view_targets_window()
            view_targets_controller.load_targets(selected_targets, self.organism, self.endonuclease)
            self.global_settings.main_window.open_new_tab("View Targets", view_targets_controller)
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in view_targets: {str(e)}")
            if self.view:
                QMessageBox.critical(self.view, "Error", f"An error occurred while viewing targets: {str(e)}")
