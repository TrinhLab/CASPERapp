import os
from PyQt6 import QtWidgets, QtCore, uic
from PyQt6.QtWidgets import QMainWindow, QMessageBox
from views.HomeWindowView import HomeWindowView
from models.HomeWindowModel import HomeWindowModel
from utils.ui import show_error, show_message
from PyQt6.QtCore import QObject
from controllers.FindTargetsController import FindTargetsController
from models.DatabaseManager import FileChangeType

class HomeWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        try:
            self.model = HomeWindowModel(global_settings)
            self.view = HomeWindowView(global_settings)
            self.init_ui()
            self.setup_connections()
            self.model.load_data()
            self.global_settings.db_manager.db_files_changed.connect(self._handle_db_files_changed)
            self.global_settings.db_manager.db_validation_changed.connect(self._handle_db_validation_changed)
            self.global_settings.db_manager.db_state_changed.connect(self._handle_db_state_changed)
        except Exception as e:
            show_error(self.global_settings, "Error initializing HomeWindowController", str(e))

    def init_ui(self):
        try:
            self.load_combo_box_data()
            self.handle_search_type_change()
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in HomeWindowController", str(e))

    def load_combo_box_data(self):
        """Reload all combo box data"""
        try:
            self.model.load_data()
            
            organism_to_endonuclease = self.model.get_organism_to_endonuclease()
            annotation_files = self.model.get_annotation_files()
            
            # Update organisms combo box
            self.view.update_combo_box_organism(sorted(organism_to_endonuclease.keys()))
            
            # Update endonuclease combo box
            self.update_combo_box_endonuclease()
            
            # Update annotation files combo box
            self.view.update_combo_box_annotation_files(annotation_files)
            
        except Exception as e:
            show_error(self.global_settings, "Error loading dropdown data", str(e))

    def update_combo_box_endonuclease(self):
        selected_organism = self.view.combo_box_organism.currentText()
        endonuclease = self.model.get_organism_to_endonuclease().get(selected_organism, [])
        self.logger.debug(f"Updating endonuclease combo box for organism {selected_organism} with endonuclease: {endonuclease} in Main window")
        self.view.update_combo_box_endonuclease(endonuclease)

    def setup_connections(self):
        try:
            # grpNavigationMenu
            self.view.push_button_new_genome.clicked.connect(self.open_new_genome_module)
            self.view.push_button_new_endonuclease.clicked.connect(self.open_new_endonuclease_module)
            self.view.push_button_multitargeting_analysis.clicked.connect(self.open_multitargeting_analysis_module)
            self.view.push_button_population_analysis.clicked.connect(self.open_population_analysis_module)

            # grpStep1
            self.view.combo_box_organism.currentIndexChanged.connect(self.update_combo_box_endonuclease)

            # grpStep2
            self.view.push_button_ncbi_file_search.clicked.connect(self.open_ncbi_window)

            # grpStep3
            # self.view.radio_button_feature.clicked.connect(self.toggle_annotation)
            # self.view.radio_button_position.clicked.connect(self.toggle_annotation)
            self.view.radio_button_feature.clicked.connect(self.handle_search_type_change)
            self.view.radio_button_position.clicked.connect(self.handle_search_type_change)
            self.view.radio_button_sequence.clicked.connect(self.handle_search_type_change)
            self.view.push_button_find_view_targets.clicked.connect(self.gather_settings)

            # Add connection for annotation file changes
            self.view.combo_box_local_annotation_files.currentTextChanged.connect(self._on_annotation_file_changed)
        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in HomeWindowController", str(e))
    

    # Event Handlers
    def gather_settings(self):
        """Process input data and direct to appropriate view"""
        try:
            input_data = self.view.get_find_targets_input()
            
            if input_data['search_type'] == 'sequence':
                sequence = input_data['search_query'].strip()
                if len(sequence) < 100:
                    QMessageBox.warning(
                        self.view,
                        "Sequence Too Short",
                        "The sequence given is too small. At least 100 characters are required."
                    )
                    return
                self.open_view_targets(input_data)
            elif input_data['search_type'] == 'position':
                self.open_view_targets(input_data)
            else:
                self.open_find_targets_module()
                
        except Exception as e:
            show_error(self.global_settings, "Error in gather_settings", str(e))

    def open_view_targets(self, input_data):
        try:
            # Create find targets controller to use its model
            find_targets_controller = self.global_settings.get_find_targets_window()
            
            # Get targets using the model
            targets = find_targets_controller.model.find_targets(input_data)
            
            if targets:
                self.logger.debug(f"Found {len(targets)} targets")
                
                # Close existing View Targets tab if it exists
                main_window = self.global_settings.main_window
                existing_tab = main_window.find_tab_by_title("View Targets")
                if existing_tab:
                    tab_index = main_window.view.tab_widget.indexOf(existing_tab)
                    main_window._close_tab(tab_index)
                    self.logger.debug("Closed existing View Targets tab")
                
                # Create view targets controller
                view_targets_controller = self.global_settings.get_view_targets_window()
                
                view_targets_controller.load_guides(
                    targets,  # Pass the targets directly
                    input_data['organism'],
                    input_data['endonuclease']
                )
                
                # Open new view targets tab
                main_window.open_new_tab(
                    "View Targets", 
                    view_targets_controller
                )
                
            else:
                QMessageBox.warning(
                    self.view,
                    "No Targets Found",
                    "No targets were found for the specified search."
                )
                
        except Exception as e:
            self.global_settings.logger.error(f"Error opening view targets directly: {str(e)}")
            show_error(self.global_settings, "Error", f"Could not open view targets: {str(e)}")

    def open_find_targets_module(self):
        """Open find targets module for non-position searches"""
        try:
            # Close existing Find Targets tab if it exists
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Find Targets")
            if existing_tab:
                tab_index = main_window.view.tab_widget.indexOf(existing_tab)
                main_window._close_tab(tab_index)
                self.logger.debug("Closed existing Find Targets tab")
            
            # Create new find targets controller and load data
            find_targets_controller = self.global_settings.get_find_targets_window()
            input_data = self.view.get_find_targets_input()
            find_targets_controller.find_targets(input_data)
            
            # Open new Find Targets tab
            self.global_settings.main_window.open_new_tab("Find Targets", find_targets_controller)
            
        except Exception as e:
            show_error(self.global_settings, "Error in open_find_targets_module() in Home", str(e))

    def toggle_annotation(self):
        # Implementation for toggling annotation
        pass

    def open_new_genome_module(self):
        try:
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("New Genome")

            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("New Genome")
            else:
                new_genome_controller = self.global_settings.get_new_genome_window()
                main_window.open_new_tab("New Genome", new_genome_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_new_genome_widget() in Home", str(e))

    def open_new_endonuclease_module(self):
        try:
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Define New Endonuclease")
            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("Define New Endonuclease")
            else:
                new_endonuclease_controller = self.global_settings.get_new_endonuclease_window()
                main_window.open_new_tab("Define New Endonuclease", new_endonuclease_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_new_endonuclease_widget() in main", str(e))

    def open_multitargeting_analysis_module(self):
        try:
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Multitargeting Analysis")
            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("Multitargeting Analysis")
            else:
                multitargeting_controller = self.global_settings.get_multitargeting_window()
                main_window.open_new_tab("Multitargeting Analysis", multitargeting_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_multitargeting_analysis_widget() in Home", str(e))

    def open_population_analysis_module(self):
        try:
            main_window = self.global_settings.main_window
            existing_tab = main_window.find_tab_by_title("Population Analysis")
            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("Population Analysis")
            else:
                population_analysis_controller = self.global_settings.get_population_analysis_window()
                main_window.open_new_tab("Population Analysis", population_analysis_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_population_analysis_widget() in Home", str(e))

    def launch_populate_fna_files(self):
        # Implementation for launching populate FNA files
        pass

    def open_ncbi_window(self):
        try:
            ncbi_controller = self.global_settings.get_ncbi_window()
            self.global_settings.main_window.open_new_tab("NCBI Download Tool", ncbi_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_ncbi_window() in main", str(e))

    def _handle_db_files_changed(self, changes):
        """Handle database file changes"""
        try:
            # Reload model data if necessary
            self.model.update_for_file_changes(changes)
            
            # Update UI if needed
            if (FileChangeType.CSPR_ADDED in changes or 
                FileChangeType.CSPR_REMOVED in changes):
                # Update both organism and endonuclease combo boxes
                organism_to_endonuclease = self.model.get_organism_to_endonuclease()
                self.view.update_combo_box_organism(sorted(organism_to_endonuclease.keys()))
                self.update_combo_box_endonuclease()
                
            if (FileChangeType.GBFF_ADDED in changes or 
                FileChangeType.GBFF_REMOVED in changes):
                self.view.update_combo_box_annotation_files(self.model.get_annotation_files())
                
        except Exception as e:
            show_error(self.global_settings, "Error handling database changes", str(e))

    def _handle_db_validation_changed(self, is_valid, message):
        """Handle database validation state changes"""
        if not is_valid:
            self.view.show_warning("Database Warning", message)
        self._update_validation_state(is_valid)

    def _handle_db_state_changed(self, is_valid, message, changes):
        """Handle database state changes"""
        try:
            if not is_valid:
                show_error(self.global_settings, "Database Warning", message)
                return
            
            # Always reload model data when database state changes
            self.model.load_data()
            
            # Update all combo boxes
            organism_to_endonuclease = self.model.get_organism_to_endonuclease()
            self.view.update_combo_box_organism(sorted(organism_to_endonuclease.keys()))
            self.update_combo_box_endonuclease()
            self.view.update_combo_box_annotation_files(self.model.get_annotation_files())
                
        except Exception as e:
            show_error(self.global_settings, "Error handling database state change", str(e))

    def _check_and_update_home_tab(self, index):
        if self.global_settings.main_window.view.tab_widget.tabText(index) == "Home":
            self.load_combo_box_data()
            # Disconnect after updating to avoid unnecessary updates
            self.global_settings.main_window.view.tab_widget.currentChanged.disconnect(self._check_and_update_home_tab)

    def get_organism_to_endonuclease(self):
        return self.model.get_organism_to_endonuclease()

    def get_annotation_files(self):
        return self.model.get_annotation_files()
    
    def get_annotation_file(self):
        return self.view.get_annotation_file()

    def _on_annotation_file_changed(self, new_file):
        """Handle changes to the annotation file selection"""
        self.global_settings.set_current_annotation_file(new_file)

    def _update_cspr_related_ui(self):
        # Implementation to update UI elements that depend on CSPR files
        pass

    def _update_gbff_related_ui(self):
        # Implementation to update UI elements that depend on GBFF files
        pass

    def _update_validation_state(self, is_valid):
        # Implementation to update UI elements based on validation state
        pass

    def handle_search_type_change(self):
        """Update UI elements based on search type"""
        try:
            search_type = self.view.get_search_type()
            
            # Update button text
            if search_type in ['position', 'sequence']:
                self.view.push_button_find_view_targets.setText("View Targets")
            else:  # 'feature'
                self.view.push_button_find_view_targets.setText("Find Targets")

        except Exception as e:
            self.logger.error(f"Error updating search type UI: {str(e)}")

