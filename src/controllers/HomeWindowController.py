from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QMessageBox
from views.HomeWindowView import HomeWindowView
from models.HomeWindowModel import HomeWindowModel
from utils.ui import show_error
from models.DatabaseManager import FileChangeType
import time
from views.LoadingDialog import LoadingDialog
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt, QSize

class HomeWindowController:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self.is_active = True
        
        try:
            self.view = HomeWindowView(self.settings)
            self.model = HomeWindowModel(self.settings)
            self._setup_connections()
            self._init_ui()
        except Exception as e:
            show_error(self.settings, "Error initializing HomeWindowController", str(e))

    def deactivate(self):
        """Cleanup controller when deactivated"""
        try:
            self.is_active = False
            if hasattr(self, 'view'):
                # Safely disconnect signals
                try:
                    self.settings.db_manager.db_state_changed.disconnect(self._on_db_state_changed)
                except (TypeError, RuntimeError):
                    # Signal wasn't connected or already disconnected
                    pass
                    
                try:
                    self.settings.db_manager.db_files_changed.disconnect(self._on_db_files_changed)
                except (TypeError, RuntimeError):
                    # Signal wasn't connected or already disconnected
                    pass
                    
                # Delete view reference
                delattr(self, 'view')
                
        except Exception as e:
            self.logger.error(f"Error in deactivate: {str(e)}")

    def _init_ui(self):
        """Initialize UI with current data"""
        try:
            self._update_ui_with_model_data()
            self.handle_search_type_change()
        except Exception as e:
            show_error(self.settings, "Error initializing UI in HomeWindowController", str(e))

    def _update_ui_with_model_data(self):
        """Update all UI elements with current model data"""
        try:
            # Get all required data at once
            organism_data = self.model.get_organism_to_endonuclease()
            annotation_files = self.model.get_annotation_files()
            
            # Update UI elements
            self._update_organism_selection(organism_data)
            self._update_annotation_files(annotation_files)
        except Exception as e:
            show_error(self.settings, "Error updating UI with model data", str(e))

    def _update_organism_selection(self, organism_data, preserve_selection=True):
        """Update organism and its dependent endonuclease selection"""
        try:
            # Store current selection if needed
            current_organism = self.view.combo_box_organism.currentText() if preserve_selection else ""
            
            # Block signals during update
            with self._block_signals(self.view.combo_box_organism):
                # Update organisms
                self.view.update_combo_box_organism(sorted(organism_data.keys()))
                
                # Restore or select first item
                if preserve_selection and current_organism in organism_data:
                    self.view.combo_box_organism.setCurrentText(current_organism)
            
            # Update endonuclease based on current organism
            self._update_endonuclease_for_organism(organism_data)
        except Exception as e:
            self.logger.error(f"Error updating organism selection: {str(e)}")

    def _update_endonuclease_for_organism(self, organism_data):
        """Update endonuclease combo box based on current organism"""
        try:
            selected_organism = self.view.combo_box_organism.currentText()
            endonucleases = organism_data.get(selected_organism, [])
            self.logger.debug(f"Updating endonuclease combo box for organism {selected_organism} with endonuclease: {endonucleases} in Main window")
            self.view.update_combo_box_endonuclease(endonucleases)
        except Exception as e:
            self.logger.error(f"Error updating endonuclease selection: {str(e)}")

    def _update_annotation_files(self, annotation_files):
        """Update annotation files combo box"""
        try:
            self.view.update_combo_box_annotation_files(annotation_files)
        except Exception as e:
            self.logger.error(f"Error updating annotation files: {str(e)}")

    def _on_organism_changed(self, _):
        """Handle organism combo box changes"""
        self._update_endonuclease_for_organism(self.model.get_organism_to_endonuclease())

    def _setup_connections(self):
        try:
            # grpNavigationMenu
            self.view.push_button_new_genome.clicked.connect(self.open_new_genome)
            self.view.push_button_new_endonuclease.clicked.connect(self.open_new_endonuclease)
            self.view.push_button_multitargeting_analysis.clicked.connect(self.open_multitargeting_analysis)
            self.view.push_button_population_analysis.clicked.connect(self.open_population_analysis)

            # grpStep1
            self.view.combo_box_organism.currentIndexChanged.connect(self._on_organism_changed)

            # grpStep2
            self.view.push_button_ncbi_file_search.clicked.connect(self.open_ncbi)

            # grpStep3
            self.view.radio_button_feature.clicked.connect(self.handle_search_type_change)
            self.view.radio_button_position.clicked.connect(self.handle_search_type_change)
            self.view.radio_button_sequence.clicked.connect(self.handle_search_type_change)
            self.view.push_button_find_view_targets.clicked.connect(self.gather_settings)

            # Add connection for annotation file changes
            self.view.combo_box_local_annotation_files.currentTextChanged.connect(self._on_annotation_file_changed)

            # Add connections for database changes
            self.settings.db_manager.db_validation_changed.connect(self._on_db_validation_changed)
            self.settings.db_manager.db_state_changed.connect(self._on_db_state_changed)
            self.settings.db_manager.db_files_changed.connect(self._on_db_files_changed)

        except Exception as e:
            show_error(self.settings, "Error setting up connections in HomeWindowController", str(e))

    def _on_db_files_changed(self, changes):
        """Handle database file changes"""
        try:
            # Reload model data if necessary
            self.model.update_for_file_changes(changes)
            
            # Update UI based on change type
            if self._should_update_organisms(changes):
                self._update_organism_selection(self.model.get_organism_to_endonuclease())
                
            if self._should_update_annotations(changes):
                self._update_annotation_files(self.model.get_annotation_files())
                
        except Exception as e:
            show_error(self.settings, "Error handling database changes", str(e))

    @staticmethod
    def _should_update_organisms(changes):
        """Check if organisms need to be updated based on changes"""
        return (FileChangeType.CSPR_ADDED in changes or 
                FileChangeType.CSPR_REMOVED in changes)

    @staticmethod
    def _should_update_annotations(changes):
        """Check if annotations need to be updated based on changes"""
        return (FileChangeType.GBFF_ADDED in changes or 
                FileChangeType.GBFF_REMOVED in changes)

    class _block_signals:
        """Context manager for blocking Qt signals"""
        def __init__(self, widget):
            self.widget = widget

        def __enter__(self):
            self.widget.blockSignals(True)
            return self.widget

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.widget.blockSignals(False)

    def refresh_data(self):
        """Refresh all data and update UI"""
        try:
            if not self.is_active or not hasattr(self, 'view'):
                return
                
            self.logger.debug("Refreshing home window data")
            self.model.load_data()
            self._update_ui_with_model_data()
            
        except Exception as e:
            self.logger.error(f"Error refreshing home window data: {str(e)}")

    def _on_db_state_changed(self, is_valid, message, changes):
        """Handle database state changes"""
        try:
            if not self.is_active or not hasattr(self, 'view'):
                return
                
            self.logger.debug(f"Database state changed - Valid: {is_valid}, Message: {message}")
            if is_valid:
                self.refresh_data()
        except Exception as e:
            self.logger.error(f"Error handling database state change: {str(e)}")

    def _check_and_update_home_tab(self, index):
        if self.settings.main_window.view.tab_widget.tabText(index) == "Home":
            self.load_combo_box_data()
            # Disconnect after updating to avoid unnecessary updates
            self.settings.main_window.view.tab_widget.currentChanged.disconnect(self._check_and_update_home_tab)

    def get_organism_to_endonuclease(self):
        return self.model.get_organism_to_endonuclease()

    def get_annotation_files(self):
        return self.model.get_annotation_files()
    
    def get_annotation_file(self):
        return self.view.get_annotation_file()

    def _on_annotation_file_changed(self, new_file):
        """Handle changes to the annotation file selection"""
        self.logger.debug(f"Current annotation file changed to: {new_file}")
        self.settings.set_current_annotation_file(new_file)

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

    def _on_db_validation_changed(self, is_valid, message):
        """Handle database validation changes"""
        try:
            if self.is_active and hasattr(self, 'view'):
                if is_valid:
                    self.refresh_data()
        except Exception as e:
            self.logger.error(f"Error handling database validation change: {str(e)}")

    def open_view_targets(self, input_data):
        try:
            # Create and show loading dialog
            loading_dialog = LoadingDialog(self.view)
            loading_dialog.show()
            loading_dialog.set_progress(0)
            QApplication.processEvents()

            try:
                # Create find targets controller to use its model
                find_targets_controller = self.settings.get_find_targets_window()
                
                # For position searches, handle each query separately
                if input_data['search_type'] == 'position':
                    queries = input_data['search_query'].strip().split('\n')
                    all_targets = []
                    total_queries = len(queries)
                    
                    for i, query in enumerate(queries):
                        # Update loading progress for each query
                        progress = int((i / total_queries) * 80)  # Leave room for final steps
                        loading_dialog.set_message(f"Processing position {i+1} of {total_queries}...", progress)
                        QApplication.processEvents()
                        
                        # Create a copy of input data with single query
                        query_data = input_data.copy()
                        query_data['search_query'] = query.strip()
                        
                        # Get targets for this query
                        targets = find_targets_controller.model.find_targets(query_data)
                        if targets:
                            # Add query information to each target
                            for target in targets:
                                target['original_query'] = query.strip()
                            all_targets.extend(targets)
                    
                    targets = all_targets  # Use combined results
                    self.logger.debug(f"Processed {len(queries)} queries, found total {len(targets)} targets")
                else:
                    # For non-position searches, process normally
                    loading_dialog.set_message("Finding targets...", 20)
                    QApplication.processEvents()
                    targets = find_targets_controller.model.find_targets(input_data)
                
                if targets:
                    self.logger.debug(f"Found {len(targets)} targets")
                    loading_dialog.set_message("Preparing view targets...", 80)
                    QApplication.processEvents()
                    
                    # Close existing View Targets tab if it exists
                    main_window = self.settings.main_window
                    existing_tab = main_window.find_tab_by_title("View Targets")
                    if existing_tab:
                        tab_index = main_window.view.tab_widget.indexOf(existing_tab)
                        main_window._close_tab(tab_index)
                        self.logger.debug("Closed existing View Targets tab")
                    
                    # Create view targets controller
                    loading_dialog.set_message("Creating view targets...", 90)
                    QApplication.processEvents()
                    view_targets_controller = self.settings.get_view_targets_window()
                    
                    view_targets_controller.load_guides(
                        targets,
                        input_data['organism'],
                        input_data['endonuclease'],
                        loading_dialog=loading_dialog
                    )
                    
                    # Open new view targets tab
                    main_window.open_new_tab(
                        "View Targets", 
                        view_targets_controller
                    )
                    
                else:
                    QtWidgets.QMessageBox.warning(
                        self.view,
                        "No Targets Found",
                        "No targets were found for the specified search."
                    )
                    
            finally:
                loading_dialog.close()
                
        except Exception as e:
            self.settings.logger.error(f"Error opening view targets directly: {str(e)}")
            show_error(self.settings, "Error", f"Could not open view targets: {str(e)}")

    def open_find_targets(self):
        """Open find targets module for non-position searches"""
        try:
            # Show loading dialog
            loading_dialog = LoadingDialog(self.view)
            loading_dialog.show()
            loading_dialog.set_progress(0)
            QApplication.processEvents()
            
            try:
                # Find all existing Find Targets tabs
                main_window = self.settings.main_window
                existing_tabs = []
                tab_numbers = []
                
                # Get all tab titles
                for i in range(main_window.view.tab_widget.count()):
                    tab_title = main_window.view.tab_widget.tabText(i)
                    if tab_title.startswith("Find Targets"):
                        existing_tabs.append(tab_title)
                        # Extract number if it exists
                        if tab_title != "Find Targets":
                            try:
                                num = int(tab_title.split()[-1])
                                tab_numbers.append(num)
                            except ValueError:
                                continue
                
                # Determine new tab number
                new_tab_number = 1
                if tab_numbers:
                    new_tab_number = max(tab_numbers) + 1
                
                loading_dialog.set_progress(40)
                
                # Create new find targets controller and load data
                find_targets_controller = self.settings.get_find_targets_window()
                input_data = self.view.get_find_targets_input()
                loading_dialog.set_progress(60)
                
                find_targets_controller.find_targets(input_data)
                loading_dialog.set_progress(80)
                
                # Open new Find Targets tab with number if not the first one
                tab_title = "Find Targets" if not existing_tabs else f"Find Targets {new_tab_number}"
                self.settings.main_window.open_new_tab(tab_title, find_targets_controller)
                loading_dialog.set_progress(100)
                
            finally:
                loading_dialog.close()
                
        except Exception as e:
            show_error(self.settings, "Error in open_find_targets() in Home", str(e))

    def open_new_genome(self):
        try:
            main_window = self.settings.main_window
            existing_tab = main_window.find_tab_by_title("New Genome")

            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("New Genome")
            else:
                new_genome_controller = self.settings.get_new_genome_window()
                main_window.open_new_tab("New Genome", new_genome_controller)
        except Exception as e:
            show_error(self.settings, "Error in open_new_genome() in Home", str(e))

    def open_new_endonuclease(self):
        try:
            # Create new endonuclease controller
            new_endonuclease_controller = self.settings.get_new_endonuclease_window()
            
            # Get the window from the controller
            window = new_endonuclease_controller.view
            
            # Set window properties
            window.setWindowModality(Qt.WindowModality.ApplicationModal)  # Make it modal
            window.setMinimumSize(QSize(500, 650))  # Smaller minimum size
            window.resize(QSize(600, 650))  # Set initial size
            
            # Get the screen where the main window is
            main_window = self.settings.main_window.view
            screen = main_window.screen()
            if not screen:
                screen = QtWidgets.QApplication.primaryScreen()
            
            # Get the available geometry of the screen (accounts for taskbars/docks)
            screen_geometry = screen.availableGeometry()
            
            # Calculate the center point of the screen
            center_point = screen_geometry.center()
            
            # Center the window on screen
            window_geometry = window.frameGeometry()
            window_geometry.moveCenter(center_point)
            window.move(window_geometry.topLeft())
            
            # Show the window
            window.show()
            window.raise_()
            window.activateWindow()
            
            # Store reference to prevent garbage collection
            self._current_new_endonuclease_window = new_endonuclease_controller
            
            self.logger.debug("New Endonuclease window opened successfully")
        except Exception as e:
            show_error(self.settings, "Error opening new endonuclease window", str(e))

    def open_multitargeting_analysis(self):
        try:
            start_time = time.time()
            self.logger.debug("Starting multitargeting analysis module launch")
            
            main_window = self.settings.main_window
            existing_tab = main_window.find_tab_by_title("Multitargeting Analysis")
            
            tab_check_time = time.time()
            self.logger.debug(f"Tab check took: {tab_check_time - start_time:.2f} seconds")

            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("Multitargeting Analysis")
                self.logger.debug(f"Switched to existing tab: {time.time() - tab_check_time:.2f} seconds")
            else:
                controller_start = time.time()
                multitargeting_controller = self.settings.get_multitargeting_window()
                self.logger.debug(f"Controller creation took: {time.time() - controller_start:.2f} seconds")
                
                tab_open_start = time.time()
                main_window.open_new_tab("Multitargeting Analysis", multitargeting_controller)
                self.logger.debug(f"Tab opening took: {time.time() - tab_open_start:.2f} seconds")

            self.logger.debug(f"Total multitargeting module launch took: {time.time() - start_time:.2f} seconds")
        except Exception as e:
            show_error(self.settings, "Error in open_multitargeting_analysis() in Home", str(e))

    def open_population_analysis(self):
        try:
            main_window = self.settings.main_window
            existing_tab = main_window.find_tab_by_title("Population Analysis")
            if existing_tab:
                main_window.view.tab_widget.setCurrentWidget(existing_tab)
                main_window._resize_for_tab("Population Analysis")
            else:
                population_analysis_controller = self.settings.get_population_analysis_window()
                main_window.open_new_tab("Population Analysis", population_analysis_controller)
        except Exception as e:
            show_error(self.settings, "Error in open_population_analysis() in Home", str(e))

    def open_ncbi(self):
        try:
            ncbi_controller = self.settings.get_ncbi_window()
            self.settings.main_window.open_new_tab("NCBI Download Tool", ncbi_controller)
        except Exception as e:
            show_error(self.settings, "Error in open_ncbi() in main", str(e))

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
                if len(sequence) > 10000:
                    QMessageBox.warning(
                        self.view,
                        "Sequence Too Long",
                        "The sequence given is too large. Maximum allowed length is 10,000 base pairs."
                    )
                    return
                self.open_view_targets(input_data)
            elif input_data['search_type'] == 'position':
                self.open_view_targets(input_data)
            else:
                self.open_find_targets()
        except Exception as e:
            show_error(self.settings, "Error in gather_settings", str(e))


