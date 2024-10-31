from models.FindTargetsModel import FindTargetsModel
from views.FindTargetsView import FindTargetsView
from PyQt6.QtWidgets import QMessageBox
from PyQt6.QtCore import QTimer

class FindTargetsController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = None
        self.view = None
        self.organism = None
        self.endonuclease = None
        self._input_data = None
        self._current_annotation_file = None
        
        # Connect to annotation file changes
        self.global_settings.annotation_file_changed.connect(self._on_annotation_file_changed)

    def _on_annotation_file_changed(self, new_annotation_file):
        """Handle annotation file changes by reprocessing data"""
        try:
            self.global_settings.logger.debug(f"FindTargetsController received new annotation file: {new_annotation_file}")
            self._current_annotation_file = new_annotation_file
            
            # Clear existing view and model
            self.view = None
            self.model = None
            
        except Exception as e:
            self.global_settings.logger.error(f"Error handling annotation file change: {str(e)}")

    def _connect_signals(self):
        """Connect view signals"""
        if self.view:
            self.view.push_button_view_targets.clicked.connect(self.view_targets)

    def find_targets(self, input_data):
        """Initialize view and process input data"""
        try:
            # Get current annotation file
            current_annotation = self.global_settings.get_current_annotation_file()
            input_data['annotation_file'] = current_annotation
            self._current_annotation_file = current_annotation
            
            # Always create new instances
            self.model = FindTargetsModel(self.global_settings)
            self.view = FindTargetsView(self.global_settings)
            self._connect_signals()
            
            self._input_data = input_data
            
            # Find existing Find Targets tab
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Find Targets")
            
            if existing_tab:
                # Remove the existing tab
                tab_index = main_window.view.tab_widget.indexOf(existing_tab)
                main_window.view.tab_widget.removeTab(tab_index)
            
            # Process data and create new tab
            self._process_input_data(input_data)
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in find_targets: {str(e)}")
            raise

    def _process_input_data(self, input_data):
        """Process input data and update view"""
        try:
            if not self.view:
                return
                
            self.global_settings.logger.debug(f"FindTargetsController processing input data: {input_data}")
            self.organism = input_data['organism']
            self.endonuclease = input_data['endonuclease']
            
            # Get new results
            results = self.model.find_targets(input_data)
            self.global_settings.logger.debug(f"Found {len(results) if results else 0} targets")
            
            # Update view with new results
            if results:
                self.view.display_results(results)
            
            # Add new tab with updated view
            main_window = self.global_settings.main_window
            main_window.open_new_tab("Find Targets", self)
            
        except Exception as e:
            self.global_settings.logger.error(f"Error processing input data: {str(e)}")
            if self.view:
                QMessageBox.critical(self.view, "Error", f"An error occurred while processing data: {str(e)}")

    def view_targets(self):
        """Handle view targets button click"""
        try:
            if not self.view:
                return
                
            selected_targets = self.view.get_selected_targets()
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
