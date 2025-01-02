from models.FindTargetsModel import FindTargetsModel
from utils.ui import show_error
from views.FindTargetsView import FindTargetsView
from PyQt6.QtWidgets import QMessageBox
from views.LoadingDialog import LoadingDialog
from PyQt6.QtWidgets import QApplication
import os

class FindTargetsController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = FindTargetsModel(self.global_settings)
        self.view = FindTargetsView(self.global_settings)
        self.organism = None
        self.endonuclease = None
        self._input_data = None
        self._current_annotation_file = None
        self.logger = self.global_settings.logger
        
        # Connect to annotation file changes
        self.global_settings.annotation_file_changed.connect(self._on_annotation_file_changed)
        self._connect_signals()

    def _on_annotation_file_changed(self, new_annotation_file):
        """Handle annotation file changes by clearing and updating results"""
        try:
            self.global_settings.logger.debug(f"FindTargetsController received new annotation file: {new_annotation_file}")
            self._current_annotation_file = new_annotation_file
            
            # Only process if we have a valid annotation file and input data
            if new_annotation_file and self._input_data:
                # Verify annotation file exists
                annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', new_annotation_file)
                if not os.path.isfile(annotation_path):
                    # Try without GBFF subdirectory
                    annotation_path = os.path.join(self.global_settings.get_db_path(), new_annotation_file)
                    if not os.path.isfile(annotation_path):
                        self.logger.warning(f"Annotation file not found at {annotation_path}")
                        return
                
                # Clear the current results
                if self.view and hasattr(self.view, 'results_table'):
                    self.view.clear_results()
                    
                    # Update input data with new annotation file
                    self._input_data['annotation_file'] = new_annotation_file
                    self._process_input_data(self._input_data)
                
        except Exception as e:
            self.global_settings.logger.error(f"Error handling annotation file change: {str(e)}")
            # Don't raise the error since this is an event handler

    def _connect_signals(self):
        """Connect view signals"""
        if self.view:
            self.view.push_button_view_targets.clicked.connect(self.view_targets)

    def find_targets(self, input_data):
        """Process input data and update existing view or create new one"""
        try:
            current_annotation = self.global_settings.get_current_annotation_file()
            input_data['annotation_file'] = current_annotation
            self._current_annotation_file = current_annotation
            self._input_data = input_data.copy()
            
            # Process data and update view
            self._process_input_data(input_data)
            
            # If there's no existing tab, create one
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Find Targets")
            if not existing_tab:
                main_window.open_new_tab("Find Targets", self)
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in find_targets: {str(e)}")
            raise

    def _process_input_data(self, input_data):
        """Process input data and update view"""
        try:
            self.global_settings.logger.debug(f"FindTargetsController processing input data: {input_data}")
            self.organism = input_data['organism']
            self.endonuclease = input_data['endonuclease']
            
            # Get new results
            results = self.model.find_targets(input_data)
            self.global_settings.logger.debug(f"Found {len(results) if results else 0} targets")
            
            if results:
                self.view.display_results(results)
            
        except Exception as e:
            self.global_settings.logger.error(f"Error processing input data: {str(e)}")
            if self.view:
                QMessageBox.critical(self.view, "Error", f"An error occurred while processing data: {str(e)}")

    def view_targets(self):
        try:
            if not self.view:
                return
            
            selected_targets = self.view.get_selected_targets()
            if not selected_targets:
                QMessageBox.warning(self.view, "No Selection", "Please select targets to view.")
                return

            # Create loading dialog with the main window as parent
            main_window = self.global_settings.main_window
            loading_dialog = LoadingDialog(main_window.view)
            loading_dialog.show()
            loading_dialog.set_progress(0)
            QApplication.processEvents()

            try:
                # Get the current Find Targets tab number
                current_tab_index = main_window.view.tab_widget.currentIndex()
                current_tab_title = main_window.view.tab_widget.tabText(current_tab_index)
                
                # Extract number from Find Targets tab (if any)
                view_targets_title = "View Targets"
                if current_tab_title != "Find Targets":
                    try:
                        number = current_tab_title.split()[-1]
                        view_targets_title = f"View Targets {number}"
                    except (IndexError, ValueError):
                        pass

                # Find existing View Targets tab with the same number
                existing_tab = main_window.find_tab_by_title(view_targets_title)
                
                loading_dialog.set_message("Initializing view targets...", 25)
                QApplication.processEvents()
                
                if existing_tab:
                    view_targets_controller = main_window.tab_widgets['controllers'].get(view_targets_title)
                    if view_targets_controller:
                        loading_dialog.set_message("Loading guides...", 50)
                        QApplication.processEvents()
                        
                        # Pass the loading dialog to load_guides
                        view_targets_controller.load_guides(
                            selected_targets, 
                            self.organism, 
                            self.endonuclease,
                            loading_dialog=loading_dialog
                        )
                        
                        # Switch to the existing tab
                        main_window.view.tab_widget.setCurrentWidget(existing_tab)
                    else:
                        self.logger.error(f"View Targets controller not found for existing tab {view_targets_title}")
                else:
                    loading_dialog.set_message("Creating view targets...", 25)
                    QApplication.processEvents()
                    
                    view_targets_controller = self.global_settings.get_view_targets_window()
                    
                    # Pass the loading dialog to load_guides
                    view_targets_controller.load_guides(
                        selected_targets, 
                        self.organism, 
                        self.endonuclease,
                        loading_dialog=loading_dialog
                    )
                    
                    main_window.open_new_tab(view_targets_title, view_targets_controller)
                    
            finally:
                loading_dialog.close()
                QApplication.processEvents()
                
        except Exception as e:
            self.logger.error(f"Error in view_targets: {str(e)}")
            if self.view:
                QMessageBox.critical(self.view, "Error", f"An error occurred while viewing targets: {str(e)}")

    def gather_settings(self):
        """Process input data and direct to appropriate view"""
        try:
            input_data = self.view.get_find_targets_input()
            
            # For position-based searches, go directly to view targets
            if input_data['search_type'] == 'position':
                self.open_view_targets_directly(input_data)
            else:
                # For other search types, show find targets view first
                self.find_targets(input_data)
                
        except Exception as e:
            show_error(self.global_settings, "Error in find_targets", str(e))

    def open_view_targets_directly(self, input_data):
        """Open view targets directly for position-based searches"""
        try:
            # Create loading dialog with the main window as parent
            main_window = self.global_settings.main_window
            loading_dialog = LoadingDialog(main_window.view)
            loading_dialog.show()
            loading_dialog.set_progress(0)
            QApplication.processEvents()

            try:
                loading_dialog.set_message("Finding targets...", 25)
                QApplication.processEvents()
                
                # Get targets using the model
                targets = self.model.find_targets_by_position(
                    self.model._get_parser(self.model.get_cspr_file_path(input_data)), 
                    input_data
                )
                
                if targets:
                    loading_dialog.set_message("Creating view targets...", 50)
                    QApplication.processEvents()
                    
                    # Create view targets controller
                    view_targets_controller = self.global_settings.get_view_targets_window()
                    
                    loading_dialog.set_message("Loading guides...", 75)
                    QApplication.processEvents()
                    
                    # Load targets directly
                    view_targets_controller.load_targets(
                        targets,
                        input_data['organism'],
                        input_data['endonuclease']
                    )
                    
                    # Open view targets tab
                    self.global_settings.main_window.open_new_tab(
                        "View Targets", 
                        view_targets_controller
                    )
                else:
                    QMessageBox.warning(
                        self.view,
                        "No Targets Found",
                        "No targets were found in the specified position range."
                    )
                    
            finally:
                loading_dialog.close()
                QApplication.processEvents()
                
        except Exception as e:
            self.global_settings.logger.error(f"Error opening view targets directly: {str(e)}")
            show_error(self.global_settings, "Error", f"Could not open view targets: {str(e)}")
