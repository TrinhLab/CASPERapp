from PyQt6.QtWidgets import QMainWindow
from views.MultitargetingWindowView import MultitargetingWindowView
from models.MultitargetingWindowModel import MultitargetingWindowModel
from utils.ui import show_error, show_message

class MultitargetingWindowController(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        
        try:
            self._model = MultitargetingWindowModel(global_settings)
            self._view = MultitargetingWindowView(global_settings)
            self.setCentralWidget(self._view)
            
            self._init_ui()
            self._setup_connections()
        except Exception as e:
            show_error(self.settings, "Error initializing MultitargetingWindowController", str(e))

    def _init_ui(self):
        """Initialize UI components"""
        try:
            # Set up initial data
            organisms = self._model.get_organisms()
            
            self._view.combo_box_organism.clear()
            self._view.combo_box_organism.addItems(organisms)
            
            # If we have organisms, trigger the first one
            if organisms:
                self._on_organism_changed(0)
            
            # Initialize plots
            self._view.setup_plots()
            
            # Connect max results line edit
            self._view.line_edit_max_results.textChanged.connect(self._on_max_results_changed)
            
        except Exception as e:
            self.logger.error(f"Error in _init_ui: {str(e)}")
            show_error(self.settings, "Error", f"Failed to initialize UI: {str(e)}")

    def _setup_connections(self):
        """Set up signal-slot connections"""
        # Organism and endonuclease selection
        self._view.combo_box_organism.currentIndexChanged.connect(self._on_organism_changed)
        self._view.combo_box_endonuclease.currentIndexChanged.connect(self._on_endonuclease_changed)
        
        # Buttons
        self._view.push_button_analyze.clicked.connect(self._on_analyze_clicked)
        # self._view.push_button_statistics_overview.clicked.connect(self._on_statistics_overview_clicked)
        # self._view.tool_button_sql_settings.clicked.connect(self._on_sql_settings_clicked)
        
        # Table selection
        self._view.table_seeds.itemSelectionChanged.connect(self._on_seed_selected)
        self._view.check_box_select_all.stateChanged.connect(self._on_select_all_changed)

        self._view.push_button_export_selected_gRNAs.clicked.connect(self._handle_export)

    def _on_organism_changed(self, index):
        """Handle organism selection change"""
        try:
            organism = self._view.combo_box_organism.currentText()
            if not organism:
                return
                
            # Get endos for selected organism
            endos = self._model.get_endos_for_organism(organism)
            
            # Update endonuclease combo box
            self._view.combo_box_endonuclease.clear()
            self._view.combo_box_endonuclease.addItems(endos)
            
            # Update file paths
            if endos:
                self._model.set_files(organism, endos[0])
                
        except Exception as e:
            self.logger.error(f"Error in _on_organism_changed: {str(e)}")
            show_error(self.settings, "Error", f"Failed to update endonucleases: {str(e)}")

    def _on_endonuclease_changed(self, index):
        """Handle endonuclease selection change"""
        self._update_analysis_button_state()

    def _on_analyze_clicked(self):
        """Handle analyze button click"""
        try:
            organism = self._view.combo_box_organism.currentText()
            endo = self._view.combo_box_endonuclease.currentText()
            
            if not organism or not endo:
                show_error(self.settings, "Analysis Error", "Please select both an organism and an endonuclease.")
                return
                
            # Load data
            try:
                self._model.set_files(organism, endo)
            except FileNotFoundError as e:
                show_error(self.settings, "File Error", 
                          f"Could not find required files for {organism} with {endo}. Please ensure the files exist.")
                return
            except ValueError as e:
                show_error(self.settings, "Input Error", str(e))
                return
                
            seeds_data = self._model.get_repeats_data()
            
            # Update UI
            self._view.update_seeds_table(seeds_data)
            self._update_plots()
            
        except Exception as e:
            show_error(self.settings, "Analysis Error", str(e))

    def _on_statistics_overview_clicked(self):
        """Handle statistics overview button click"""
        try:
            stats = self._model.calculate_statistics()
            self._show_statistics_dialog(stats)
        except Exception as e:
            show_error(self.settings, "Statistics Error", str(e))

    def _on_sql_settings_clicked(self):
        """Handle SQL settings button click"""
        try:
            current_settings = self._model.get_sql_settings()
            if self._show_sql_settings_dialog(current_settings):
                new_settings = self._get_sql_settings_from_dialog()
                self._model.update_sql_settings(new_settings)
        except Exception as e:
            show_error(self.settings, "SQL Settings Error", str(e))

    def _on_seed_selected(self):
        """Handle seed selection in table"""
        try:
            selected_items = self._view.table_seeds.selectedItems()
            if selected_items:
                row = selected_items[0].row()
                seed = self._view.table_seeds.item(row, 0).text()
                
                # Get seed data
                seed_data = self._model.get_seed_data(seed)
                if not seed_data:
                    return

                # Process seed data for visualization
                kstats = self._model.get_kstats()
                seed_data_processed, event_data = self._process_seed_data(seed_data, kstats)
                
                # Update chromosome viewer
                self._view.fill_chromosome_viewer(seed_data_processed, event_data)
                
                # Update only the chromosome bar plot
                chromosome_data = self._model.get_chro_bar_data(seed)
                # Update only the chromosome plot, keep other plots unchanged
                self._view._update_repeat_vs_chromosome_plot(chromosome_data)
                
        except Exception as e:
            self.logger.error(f"Error handling seed selection: {str(e)}")
            show_error(self.settings, "Error", f"Failed to display seed data: {str(e)}")

    def _process_seed_data(self, seed_data, kstats):
        """Process seed data for visualization"""
        try:
            seed_data_processed = {}
            event_data = {}
            
            for data in seed_data:
                # Split chromosome and location strings into lists
                chromos = [int(x) for x in data[0].split(',')]
                locs = [int(x) for x in data[1].split(',')]
                pams = data[2].split(',')
                scores = data[3].split(',')
                fives = data[4].split(',')
                threes = data[5].split(',')
                
                # Process each chromosome location
                for i in range(len(chromos)):
                    chromo = chromos[i]
                    pos = locs[i]
                    
                    # Normalize location
                    dir = "+" if pos >= 0 else "-"
                    normalized_location = abs(float(pos) / float(kstats[chromo - 1]))
                    
                    # Store data
                    if chromo in seed_data_processed:
                        seed_data_processed[chromo].append(normalized_location)
                        event_data[chromo].append([
                            normalized_location, 
                            pos, 
                            fives[i] + threes[i], 
                            pams[i], 
                            scores[i], 
                            dir
                        ])
                    else:
                        seed_data_processed[chromo] = [normalized_location]
                        event_data[chromo] = [[
                            normalized_location, 
                            pos, 
                            fives[i] + threes[i], 
                            pams[i], 
                            scores[i], 
                            dir
                        ]]
                        
            return seed_data_processed, event_data
            
        except Exception as e:
            self.logger.error(f"Error processing seed data: {str(e)}")
            raise

    def _on_select_all_changed(self, state):
        """Handle select all checkbox state change"""
        self._view.table_seeds.selectAll() if state else self._view.table_seeds.clearSelection()

    def _update_analysis_button_state(self):
        """Update analyze button enabled state"""
        has_organism = bool(self._view.combo_box_organism.currentText())
        has_endo = bool(self._view.combo_box_endonuclease.currentText())
        self._view.push_button_analyze.setEnabled(has_organism and has_endo)
        
        # Clear any existing data if selection changes
        if not (has_organism and has_endo):
            self._view.table_seeds.setRowCount(0)
            self._view.update_plots(None, None, None)

    def _update_plots(self):
        try:
            repeats_data = self._model.get_repeats_vs_seeds_data()
            sequences_data = self._model.get_seeds_vs_repeats_data()

            # Get statistics for overview tab
            stats = self._model.calculate_statistics()
            
            # Update statistics labels
            if stats:
                self._view.update_statistics_labels(
                    total_repeats=stats.get('repeat_count', 0),
                    avg_repeats=stats.get('average', 0),
                    median_repeats=stats.get('median', 0),
                    mode_repeats=stats.get('mode', 0)
                )

            # Update all plots at once
            self._view.update_plots(repeats_data, sequences_data, None)  # chromosome_data will be updated on seed selection
            
        except Exception as e:
            self.logger.error(f"Error in _update_plots: {str(e)}")
            show_error(self.settings, "Plot Update Error", str(e))

    def _show_statistics_dialog(self, stats):
        """Show statistics overview dialog"""
        # Implement statistics dialog display
        pass

    def _show_sql_settings_dialog(self, current_settings):
        """Show SQL settings dialog"""
        # Implement SQL settings dialog display
        return False

    def _get_sql_settings_from_dialog(self):
        """Get settings from SQL settings dialog"""
        # Implement getting settings from dialog
        return {}

    def _on_max_results_changed(self, value):
        """Handle changes to max results setting"""
        try:
            if value == "":  # Handle empty input
                self._model.set_row_limit(1000)  # Reset to default
                return
            
            # Convert to int and update model
            limit = int(value)
            if limit <= 0:  # Handle negative or zero values
                limit = -1  # Use -1 to indicate no limit
            self._model.set_row_limit(limit)
            
        except ValueError:
            # Reset to default if invalid input
            self._model.set_row_limit(1000)
            self._view.line_edit_max_results.setText("1000")

    def _handle_export(self):
        """Handle export button click"""
        try:
            selected_items = []
            selected_rows = self._view.table_seeds.selectedItems()
            
            if not selected_rows:
                show_message(
                    "Warning",
                    "Please select at least one row to export."
                )
                return

            # Get unique rows (since selecting one row selects all its columns)
            selected_row_numbers = set()
            for item in selected_rows:
                selected_row_numbers.add(item.row())

            # For each selected row, create a dictionary with the row data
            for row in selected_row_numbers:
                item_data = {
                    'seed': self._view.table_seeds.item(row, 0).text(),
                    'total_repeats': self._view.table_seeds.item(row, 1).text(),
                    'avg_repeats_or_scaffold': self._view.table_seeds.item(row, 2).text(),
                    'consensus_sequence': self._view.table_seeds.item(row, 3).text(),
                    'percent_consensus': self._view.table_seeds.item(row, 4).text(),
                    'score': self._view.table_seeds.item(row, 5).text(),
                    'pam': self._view.table_seeds.item(row, 6).text(),
                    'strand': self._view.table_seeds.item(row, 7).text()
                }
                selected_items.append(item_data)

            # Get export window controller and show dialog
            export_controller = self.settings.get_export_selected_grnas_window()
            export_controller.show_dialog(selected_items, "Multitargeting")

        except Exception as e:
            self.logger.error(f"Error handling export: {str(e)}")
            show_error(self.settings, "Export Error", str(e))
