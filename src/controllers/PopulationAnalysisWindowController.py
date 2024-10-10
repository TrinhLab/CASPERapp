from PyQt6 import QtWidgets
from utils.ui import show_error, show_message, position_window
from views.PopulationAnalysisWindowView import PopulationAnalysisWindowView
from models.PopulationAnalysisWindowModel import PopulationAnalysisWindowModel

class PopulationAnalysisWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = PopulationAnalysisWindowModel(global_settings)
        self.view = PopulationAnalysisWindowView(global_settings)
        self.setup_connections()

    def setup_connections(self):
        self.view.goBackButton.clicked.connect(self.go_back)
        self.view.analyze_button.clicked.connect(self.pre_analyze)
        self.view.clear_Button.clicked.connect(self.clear)
        self.view.export_button.clicked.connect(self.export_tool)
        self.view.find_locs_button.clicked.connect(self.find_locations)
        self.view.clear_loc_button.clicked.connect(self.clear_loc_table)
        self.view.query_seed_button.clicked.connect(self.custom_seed_search)
        self.view.endoBox.currentIndexChanged.connect(self.change_endo)
        self.view.table2.horizontalHeader().sectionClicked.connect(self.table2_sorting)
        self.view.loc_finder_table.horizontalHeader().sectionClicked.connect(self.loc_table_sorter)

    def launch(self):
        try:
            self.get_data()
            self.view.show()
        except Exception as e:
            show_error(self.global_settings, "Error in launch() in population analysis.", str(e))

    def get_data(self):
        try:
            self.fillEndo()
        except Exception as e:
            show_error(self.global_settings, "Error in get_data() in population analysis.", str(e))

    def fillEndo(self):
        try:
            endos = self.model.load_endonucleases()
            self.view.update_endo_dropdown(endos.keys())
            self.change_endo()
        except Exception as e:
            show_error(self.global_settings, "Error in fillEndo() in population analysis.", str(e))

    def change_endo(self):
        try:
            selected_endo = self.view.get_selected_endo()
            org_files = self.model.get_organism_files(selected_endo)
            self.view.update_org_table(org_files)
        except Exception as e:
            show_error(self.global_settings, "Error in change_endo() in population analysis.", str(e))

    def pre_analyze(self):
        try:
            selected_indexes = self.view.get_selected_organisms()
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
            self.view.show_loading_window(5)
            self.model.seeds = self.model.get_shared_seeds(self.model.db_files, True)
            
            if len(self.model.seeds) == 0:
                self.view.hide_loading_window()
                return

            seed_data = []
            for seed in self.model.seeds:
                data = self.model.get_seed_data(seed, self.model.db_files)
                seed_data.append(self.process_seed_data(seed, data))

            self.view.update_shared_seeds_table(seed_data)
            
            if len(self.model.db_files) > 1:
                heatmap_data = self.model.get_heatmap_data(self.model.db_files)
                self.view.plot_heatmap(heatmap_data, self.model.org_names)

            self.view.hide_loading_window()
        except Exception as e:
            show_error(self.global_settings, "Error in fill_data() in population analysis.", str(e))

    def process_seed_data(self, seed, data):
        # Process the seed data and return a list of values for the table
        # This method should contain the logic to calculate percentages, averages, etc.
        # Return a list that matches the columns in the shared seeds table
        pass

    def custom_seed_search(self):
        try:
            seeds = self.view.get_seed_input().split(',')
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

            self.view.update_shared_seeds_table(seed_data)
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

            locations = self.model.get_seed_locations(selected_seeds, self.model.db_files)
            self.view.update_loc_finder_table(locations)
        except Exception as e:
            show_error(self.global_settings, "Error in find_locations() in population analysis.", str(e))

    def clear_loc_table(self):
        try:
            self.view.clear_loc_finder_table()
        except Exception as e:
            show_error(self.global_settings, "Error in clear_loc_table() in population analysis.", str(e))

    def table2_sorting(self, logicalIndex):
        try:
            self.view.sort_table2(logicalIndex)
        except Exception as e:
            show_error(self.global_settings, "Error in table2_sorting() in population analysis.", str(e))

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
