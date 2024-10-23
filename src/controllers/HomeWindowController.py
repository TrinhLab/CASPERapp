import os
from PyQt6 import QtWidgets, QtCore, uic
from PyQt6.QtWidgets import QMainWindow
from views.HomeWindowView import HomeWindowView
from models.HomeWindowModel import HomeWindowModel
from utils.ui import show_error, show_message, scale_ui, center_ui, position_window
from PyQt6.QtCore import QObject
from controllers.FindTargetsController import FindTargetsController

class HomeWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        try:
            self.model = HomeWindowModel(global_settings)
            self.view = HomeWindowView(global_settings)
            self.init_ui()
            self.setup_connections()
            self.global_settings.db_state_updated.connect(self.refresh_data)

            main_window = self.global_settings.main_window
            
            # Ensure the model data is loaded
            self.model.load_data()
            # self.find_targets_controller = FindTargetsController(global_settings)
        except Exception as e:
            show_error(self.global_settings, "Error initializing HomeWindowController", str(e))

    def init_ui(self):
        try:
            self.view.push_button_view_targets.setEnabled(False)
            self.view.push_button_generate_library.setEnabled(False)
            self.load_combo_box_data()
            self.view.reset_progress_bar()
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in HomeWindowController", str(e))

    def load_combo_box_data(self):
        try:
            self.model.load_data()

            organism_to_endonuclease = self.model.get_organism_to_endonuclease()
            annotation_files = self.model.get_annotation_files()

            self.logger.debug(f"Updating Organisms combo box with organisms: {organism_to_endonuclease.keys()} in Main window")
            self.view.update_combo_box_organism(list(organism_to_endonuclease.keys()))

            self.update_combo_box_endonuclease()

            self.logger.debug(f"Updating Annotation files combo box with annotation files: {annotation_files} in Main window")
            self.view.update_combo_box_annotation_files(annotation_files)
        except Exception as e:
            show_error(self.global_settings, "Error loading dropdown data in HomeWindowController", str(e))

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
            self.view.radio_button_feature.clicked.connect(self.toggle_annotation)
            self.view.radio_button_position.clicked.connect(self.toggle_annotation)
            self.view.push_button_find_targets.clicked.connect(self.gather_settings)
            self.view.push_button_view_targets.clicked.connect(self.view_results)
            self.view.push_button_generate_library.clicked.connect(self.prep_gen_lib)
        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in HomeWindowController", str(e))
    

    # Event Handlers
    def gather_settings(self):
        try:
            # input_data = self.view.get_find_targets_input()
            # self.find_targets_controller.find_targets(input_data)
            # self.global_settings.main_window.open_new_tab("Find Targets", self.find_targets_controller)

            self.open_find_targets_module()
        except Exception as e:
            show_error(self.global_settings, "Error in find_targets", str(e))

    def view_results(self):
        # Implementation for viewing results
        pass

    def toggle_annotation(self):
        # Implementation for toggling annotation
        pass

    def prep_gen_lib(self):
        # Implementation for preparing gene library
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
            multitargeting_controller = self.global_settings.get_multitargeting_window()
            multitargeting_window = multitargeting_controller.view
            self.global_settings.main_window.open_new_tab("Multitargeting Analysis", multitargeting_window)
        except Exception as e:
            show_error(self.global_settings, "Error in open_multitargeting_analysis_widget() in Home", str(e))

    def open_population_analysis_module(self):
        try:
            population_analysis_controller = self.global_settings.get_population_analysis_window()
            population_analysis_window = population_analysis_controller.view
            self.global_settings.main_window.open_new_tab("Population Analysis", population_analysis_window)
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

    def open_find_targets_module(self):
        try:
            find_targets_controller = self.global_settings.get_find_targets_window()
            input_data = self.view.get_find_targets_input()
            find_targets_controller.find_targets(input_data)
            self.global_settings.main_window.open_new_tab("Find Targets", find_targets_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_find_targets_module() in Home", str(e))

    def open_view_targets_module(self):
        try:
            view_targets_controller = self.global_settings.get_view_targets_window()
            self.global_settings.main_window.open_new_tab("View Targets", view_targets_controller)
        except Exception as e:
            show_error(self.global_settings, "Error in open_view_targets_module() in Home", str(e))

    def refresh_data(self, is_valid, message, cspr_files):
        self.logger.info(f"Refreshing Home Window data after database state update. Valid: {is_valid}, Message: {message}")
        if is_valid:
            self.load_combo_box_data()
            # If the current tab is not Home, we need to update it when it becomes visible
            main_window = self.global_settings.main_window
            current_tab_text = main_window.view.tab_widget.tabText(main_window.view.tab_widget.currentIndex())
            if current_tab_text != "Home":
                main_window.view.tab_widget.currentChanged.connect(self._check_and_update_home_tab)
        else:
            self.logger.warning(f"Database state update received, but it's not valid. Message: {message}")

    def _check_and_update_home_tab(self, index):
        if self.global_settings.main_window.view.tab_widget.tabText(index) == "Home":
            self.load_combo_box_data()
            # Disconnect after updating to avoid unnecessary updates
            self.global_settings.main_window.view.tab_widget.currentChanged.disconnect(self._check_and_update_home_tab)

    def get_organism_to_endonuclease(self):
        return self.model.get_organism_to_endonuclease()

    def get_annotation_files(self):
        return self.model.get_annotation_files()

