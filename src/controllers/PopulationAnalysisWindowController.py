from PyQt6 import QtWidgets
from utils.ui import show_error, show_message
from views.PopulationAnalysisWindowView import PopulationAnalysisWindowView
from models.PopulationAnalysisWindowModel import PopulationAnalysisWindowModel
import logging

class PopulationAnalysisWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = PopulationAnalysisWindowModel(self.global_settings)
        self.view = PopulationAnalysisWindowView(self.global_settings)
        self._init_ui()
        self.setup_connections()

    def setup_connections(self):
        try:
            # Group: Select Organisms
            self.view.combo_box_endonuclease.currentIndexChanged.connect(self.change_endo)
            self.view.push_button_analyze_organism.clicked.connect(self.pre_analyze)
            
            # Group: Seed Analysis
            self.view.push_button_query_seed.clicked.connect(self.custom_seed_search)
            self.view.push_button_clear_seeds.clicked.connect(self.clear)
            self.view.push_button_find_locations.clicked.connect(self.find_locations)
            self.view.push_button_clear_locations.clicked.connect(self.clear_loc_table)
            
            # Tables sorting
            self.view.table_seed.horizontalHeader().sectionClicked.connect(self.seed_table_sorting)
            self.view.table_locations.horizontalHeader().sectionClicked.connect(self.loc_table_sorter)
            
            # Add new connection for export button
            self.view.push_button_export_selected_gRNAs.clicked.connect(self.export_selected_seeds)
        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in population analysis.", str(e))

    def _init_ui(self):
        self.launch()

    def launch(self):
        try:
            self.logger = self.global_settings.get_logger()
            self.logger.info("Launching Population Analysis Window")
            self.get_data()
        except Exception as e:
            self.logger.error(f"Error in launch(): {str(e)}")
            show_error(self.global_settings, "Error in launch() in population analysis.", str(e))

    def get_data(self):
        try:
            self.logger.info("Getting data for Population Analysis")
            self.fillEndo()
        except Exception as e:
            self.logger.error(f"Error in get_data(): {str(e)}")
            show_error(self.global_settings, "Error in get_data() in population analysis.", str(e))

    def fillEndo(self):
        try:
            self.logger.info("Starting fillEndo()")
            endos = self.model.load_endonucleases()
            self.logger.debug(f"Loaded endonucleases: {endos}")
            
            if not endos:
                self.logger.warning("No endonucleases found")
                show_error(self.global_settings, "Error", "No endonucleases found")
                return
            
            self.logger.info(f"Updating dropdown with {len(endos)} endonucleases")
            self.view.update_endo_dropdown(endos.keys())
            self.change_endo()
        except Exception as e:
            self.logger.error(f"Error in fillEndo(): {str(e)}")
            show_error(self.global_settings, "Error in fillEndo() in population analysis.", str(e))

    def change_endo(self):
        try:
            selected_endo = self.view.get_selected_endo()
            if not selected_endo:
                return
            org_files = self.model.get_organism_files(selected_endo)
            self.view.update_org_table(org_files)
        except Exception as e:
            show_error(self.settings, "Error in change_endo() in population analysis.", str(e))

    def pre_analyze(self):
        try:
            selected_indexes = [index.row() for index in self.view.table_organism.selectionModel().selectedRows()]
            if len(selected_indexes) == 0:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Error",
                    message="Please select CSPR file(s) for analysis."
                )
                return

            self.model.cspr_files = [self.model.index_to_cspr[index] for index in selected_indexes]
            self.model.db_files = [self.model.index_to_db[index] for index in selected_indexes]

            self.model.get_org_names()
            self.fill_data()
        except Exception as e:
            show_error(self.global_settings, "Error in pre_analyze() in population analysis.", str(e))

    def fill_data(self):
        try:
            # Get seeds shared between ALL organisms
            self.model.seeds = self.model.get_shared_seeds(self.model.db_files, True)
            
            if len(self.model.seeds) > 0:
                # Process seed data for the table
                seed_data = []
                for seed in self.model.seeds:
                    data = self.model.get_seed_data(seed, self.model.db_files)
                    processed_data = self.process_seed_data(seed, data)
                    if processed_data:  # Only add if data was processed successfully
                        seed_data.append(processed_data)

                if seed_data:  # Only update table if we have data
                    self.view.update_shared_seeds_table(seed_data)
            
            # Always generate and display heatmap for 2 or more organisms
            if len(self.model.db_files) > 1:
                heatmap_data = self.model.get_heatmap_data(self.model.db_files)
                self.view.plot_heatmap(heatmap_data, self.model.org_names)
            else:
                self.logger.warning("Not enough organisms selected for heatmap")

        except Exception as e:
            show_error(self.global_settings, "Error in fill_data() in population analysis.", str(e))

    def process_seed_data(self, seed, data):
        """Process seed data and return a tuple of values for the table"""
        try:
            # self.logger.debug(f"Processing seed data: {data}")
            
            if not data or data['org_count'] == 0:
                self.logger.warning(f"No data found for seed {seed}")
                return None

            # Calculate coverage percentage
            coverage = (data['org_count'] / len(self.model.db_files)) * 100
            coverage = float("%.2f" % coverage)

            # Calculate average repeats per scaffold
            avg_rep_per_scaff = data['total_count'] / data['org_count']
            avg_rep_per_scaff = float("%.2f" % avg_rep_per_scaff)

            # Handle missing data in threes/fives
            threes = data['threes']
            fives = data['fives']
            if len(threes) < len(fives):
                threes.extend([''] * (len(fives) - len(threes)))
            elif len(fives) < len(threes):
                fives.extend([''] * (len(threes) - len(fives)))

            # Find majority sequence
            majority_index = 0
            if not threes or threes[0] == '':
                majority = max(set(fives), key=fives.count)
                majority_index = fives.index(majority)
                consensus_seq = fives[majority_index] + seed
                percent_consensus = (fives.count(fives[majority_index]) / len(fives)) * 100
            elif not fives or fives[0] == '':
                majority = max(set(threes), key=threes.count)
                majority_index = threes.index(majority)
                consensus_seq = seed + threes[majority_index]
                percent_consensus = (threes.count(threes[majority_index]) / len(threes)) * 100
            else:
                # Both threes and fives present
                combined = [f"{f}{t}" for f, t in zip(fives, threes)]
                majority = max(set(combined), key=combined.count)
                majority_index = combined.index(majority)
                consensus_seq = fives[majority_index] + seed + threes[majority_index]
                percent_consensus = (combined.count(majority) / len(combined)) * 100

            percent_consensus = float("%.2f" % percent_consensus)

            # Determine strand
            strand = "+" if int(data['locs'][majority_index]) >= 0 else "-"

            # Create the row data
            row_data = (
                seed,                   # Seed
                coverage,              # % Coverage
                data['total_count'],   # Total Repeats
                avg_rep_per_scaff,     # Avg. Repeats/Scaffold
                consensus_seq,         # Consensus Sequence
                percent_consensus,     # % Consensus
                data['scores'][majority_index],  # Score
                data['pams'][majority_index],    # PAM
                strand                 # Strand
            )

            self.logger.debug(f"Processed seed data: {row_data}")
            return row_data

        except Exception as e:
            self.logger.error(f"Error processing seed data: {str(e)}")
            show_error(self.global_settings, f"Error processing seed {seed}", str(e))
            return None

    def custom_seed_search(self):
        try:
            seeds = self.view.line_edit_seed.text().split(',')
            seeds = [seed.strip().upper() for seed in seeds if seed.strip()]

            if not seeds:
                self.pre_analyze()
                return

            seed_data = []
            for seed in seeds:
                data = self.model.get_seed_data(seed, self.model.db_files)
                if data['org_count'] > 0:
                    seed_data.append(self.process_seed_data(seed, data))
                else:
                    show_message(
                        fontSize=12,
                        icon=QtWidgets.QMessageBox.Icon.Critical,
                        title="Seed Error",
                        message=f"{seed}: No such seed exists in the repeats section of any organism selected."
                    )
                    return

            self.view.update_seed_table(seed_data)
        except Exception as e:
            show_error(self.global_settings, "Error in custom_seed_search() in population analysis.", str(e))

    def find_locations(self):
        try:
            selected_seeds = self.view.get_selected_seeds()
            if not selected_seeds:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Error",
                    message="Please select at least 1 seed to find locations of."
                )
                return

            # Clear the locations table before adding new entries
            self.view.table_locations.setRowCount(0)

            # Get and display new locations
            locations = self.model.get_seed_locations(selected_seeds, self.model.db_files)
            self.view.update_loc_finder_table(locations)
        except Exception as e:
            show_error(self.global_settings, "Error in find_locations() in population analysis.", str(e))

    def clear_loc_table(self):
        try:
            self.view.clear_loc_finder_table()
        except Exception as e:
            show_error(self.global_settings, "Error in clear_loc_table() in population analysis.", str(e))

    def seed_table_sorting(self, logicalIndex):
        try:
            self.view.sort_table2(logicalIndex)
        except Exception as e:
            show_error(self.global_settings, "Error in seed_table_sorting() in population analysis.", str(e))

    def loc_table_sorter(self, logicalIndex):
        try:
            self.view.sort_loc_finder_table(logicalIndex)
        except Exception as e:
            show_error(self.global_settings, "Error in loc_table_sorter() in population analysis.", str(e))

    def clear(self):
        try:
            self.view.clear_shared_seeds_table()
        except Exception as e:
            show_error(self.global_settings, "Error in clear() in population analysis.", str(e))

    def go_back(self):
        try:
            self.global_settings.main_window.show()
            self.view.hide()
        except Exception as e:
            show_error(self.global_settings, "Error in go_back() in population analysis.", str(e))

    def export_tool(self):
        try:
            selected_items = self.view.get_selected_seeds_for_export()
            if not selected_items:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Nothing Selected",
                    message="No targets were highlighted. Please highlight the targets you want to be exported to a CSV File!"
                )
                return

            self.global_settings.main_window.export_tool_window.launch(selected_items, "pa")
        except Exception as e:
            show_error(self.global_settings, "Error in export_tool() in population analysis.", str(e))

    def show(self):
        self.view.show()

    def hide(self):
        self.view.hide()

    def closeEvent(self, event):
        try:
            self.global_settings.main_window.closeFunction()
            event.accept()
        except Exception as e:
            show_error(self.global_settings, "Error in closeEvent() in population analysis.", str(e))

    def export_selected_seeds(self):
        try:
            selected_seeds = []
            selected_rows = self.view.table_seed.selectionModel().selectedRows()
            self.logger.debug(f"Selected rows: {selected_rows}")
            
            if not selected_rows:
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Nothing Selected",
                    message="No seeds were selected. Please select seeds to export."
                )
                return
            
            for row_idx in selected_rows:
                seed_data = {
                    # Seed, % Coverage, Total Repeats, Avg. Repeats/Scaffold, Consensus Sequence, Full Sequence, % Consensus, Score, PAM, Strand
                    # ('TCCCTGGTTCGAATCC', 100.0, 4, 2.0, 'TTGGTCCCTGGTTCGAATCC', 50.0, '55', 'GGG', '-')
                    'seed': self.view.table_seed.item(row_idx.row(), 0).text(),  # Seed column
                    'percent_coverage': self.view.table_seed.item(row_idx.row(), 1).text(),  # % Coverage column
                    'total_repeats': self.view.table_seed.item(row_idx.row(), 2).text(),  # Total Repeats column
                    'avg_repeats_or_scaffold': self.view.table_seed.item(row_idx.row(), 3).text(),  # Avg. Repeats/Scaffold column
                    'consensus_sequence': self.view.table_seed.item(row_idx.row(), 4).text(),  # Consensus Sequence column
                    'full_sequence': self.view.table_seed.item(row_idx.row(), 4).text(),  # Full Sequence column
                    'percent_consensus': self.view.table_seed.item(row_idx.row(), 5).text(),  # % Consensus column
                    'score': self.view.table_seed.item(row_idx.row(), 6).text(),    # Score column
                    'pam': self.view.table_seed.item(row_idx.row(), 7).text(),      # PAM column
                    'strand': self.view.table_seed.item(row_idx.row(), 8).text(),   # Strand column
                }
                selected_seeds.append(seed_data)
            
            # Get export window from global settings and show dialog
            export_window = self.global_settings.get_export_selected_grnas_window()
            print(f"Selected seeds: {selected_seeds}")
            export_window.show_dialog(selected_seeds, "Population Analysis")
            
        except Exception as e:
            show_error(self.global_settings, "Error exporting selected seeds", str(e))
