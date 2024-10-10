import os
from PyQt6 import QtWidgets, QtCore, uic
from PyQt6.QtWidgets import QMainWindow
from views.HomeWindowView import HomeWindowView
from models.HomeWindowModel import HomeWindowModel
from utils.ui import show_error, show_message, scale_ui, center_ui, position_window
from PyQt6.QtCore import QObject

class HomeWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        try:
            self.model = HomeWindowModel(global_settings)
            self.view = HomeWindowView(global_settings)
            self.tab_widgets = {}  # Store references to tab widgets
            self.init_ui()
            self.setup_connections()
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
        # Implementation for gathering settings
        pass

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
            new_genome_controller = self.global_settings.get_new_genome_window()
            self.global_settings.main_window.open_new_tab("New Genome", new_genome_controller)
        except Exception as e:
            self.logger.error(f"Error in open_new_genome_module: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error in open_new_genome_widget() in Home", str(e))

    def open_new_endonuclease_module(self):
        try:
            new_endonuclease_controller = self.global_settings.get_new_endonuclease_window()
            self.global_settings.main_window.open_new_tab("Define New Endonuclease", new_endonuclease_controller)
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
            ncbi_window = ncbi_controller.view
            self.global_settings.main_window.open_new_tab("NCBI Download Tool", ncbi_window)
        except Exception as e:
            show_error(self.global_settings, "Error in open_ncbi_window() in main", str(e))