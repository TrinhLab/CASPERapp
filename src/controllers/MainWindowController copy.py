import os
from PyQt6 import QtWidgets, QtCore
from PyQt6.QtWidgets import QMainWindow
from views.MainWindowView import MainWindowView
from models.MainWindowModel import MainWindowModel
from controllers.MultitargetingWindowController import MultitargetingWindowController
from utils.ui import show_error, show_message, scale_ui, center_ui, position_window
from utils.web import ncbi_page, repo_page, ncbi_blast_page
from PyQt6.QtCore import QObject

class MainWindowController(QObject):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        try:
            self.model = MainWindowModel(global_settings)
            self.view = MainWindowView(global_settings)
            # self.setup_connections()
            # self.init_ui()
            self.show()
        except Exception as e:
            show_error(global_settings, "Error initializing MainWindowController", str(e))
            raise

    def setup_connections(self):
        try:
            # menuBar
            self.view.action_change_directory.triggered.connect(self.change_database_directory)
            # self.view.action_exit.triggered.connect(self.close)
            # self.view.action_open_genome_browser.triggered.connect(self.open_genome_browser)
            self.view.action_open_repository.triggered.connect(self.open_repository_website)
            self.view.action_open_NCBI_BLAST.triggered.connect(self.open_ncbi_blast_website)
            self.view.action_open_NCBI.triggered.connect(self.open_ncbi_website)

            # grpNavigationMenu
            self.view.push_button_new_genome.clicked.connect(self.open_new_genome_widget)
            self.view.push_button_new_endonuclease.clicked.connect(self.open_new_endonuclease_widget)
            self.view.push_button_multitargeting_analysis.clicked.connect(self.open_multitargeting_analysis_widget)
            self.view.push_button_population_analysis.clicked.connect(self.open_population_analysis_widget)

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
            # self.view.theme_toggle_button.clicked.connect(self.toggle_theme)

            # Custom title bar connections
            self.view.theme_toggle_button.clicked.connect(self.toggle_theme)
            self.view.minimize_button.clicked.connect(self.view.showMinimized)
            self.view.maximize_button.clicked.connect(self.toggle_maximize)
            self.view.close_button.clicked.connect(self.view.close)
        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in MainWindowController", str(e))
    
    def init_ui(self):
        try:
            self.view.push_button_view_targets.setEnabled(False)
            self.view.push_button_generate_library.setEnabled(False)
            self.load_combo_box_data()
            self.view.reset_progress_bar()
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in MainWindowController", str(e))

    def load_combo_box_data(self):
        try:
            self.model.load_data()
            self.update_combo_box_data()
        except Exception as e:
            show_error(self.global_settings, "Error loading dropdown data in MainWindowController", str(e))

    def update_combo_box_data(self):
        organism_to_endonuclease = self.model.get_organism_to_endonuclease()
        annotation_files = self.model.get_annotation_files()

        self.logger.debug(f"Updating Organisms combo box with organisms: {organism_to_endonuclease.keys()} in Main window")
        self.view.update_combo_box_organism(list(organism_to_endonuclease.keys()))

        self.update_combo_box_endonuclease()

        self.logger.debug(f"Updating Annotation files combo box with annotation files: {annotation_files} in Main window")
        self.view.update_combo_box_annotation_files(annotation_files)

    def update_combo_box_endonuclease(self):
        selected_organism = self.view.combo_box_organism.currentText()
        endonuclease = self.model.get_organism_to_endonuclease().get(selected_organism, [])
        self.logger.debug(f"Updating endonuclease combo box for organism {selected_organism} with endonuclease: {endonuclease} in Main window")
        self.view.update_combo_box_endonuclease(endonuclease)

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

    def open_new_genome_widget(self):
        try:
            new_genome_window = self.global_settings.new_genome_window
            self._open_widget("New Genome", new_genome_window.view)
        except Exception as e:
            show_error(self.global_settings, "Error in open_new_genome_widget() in main", e)

    def open_new_endonuclease_widget(self):
        try:
            new_endo_window = self.global_settings.new_endonuclease_window
            self._open_widget("New Endonuclease", new_endo_window.view)
        except Exception as e:
            show_error(self.global_settings, "Error in open_new_endonuclease_widget() in main", str(e))

    def open_multitargeting_analysis_widget(self):
        try:
            multitargeting_window = self.global_settings.multitargeting_window
            self._open_widget("Multitargeting Analysis", multitargeting_window)
        except Exception as e:
            show_error(self.global_settings, "Error in open_multitargeting_analysis_widget() in main", str(e))

    def open_population_analysis_widget(self):
        try:
            population_analysis_window = self.global_settings.population_analysis_window
            self._open_widget("Population Analysis", population_analysis_window)
        except Exception as e:
            show_error(self.global_settings, "Error in open_population_analysis_widget() in main", str(e))

    def launch_populate_fna_files(self):
        # Implementation for launching populate FNA files
        pass

    def open_ncbi_window(self):
        try:
            ncbi_window = self.global_settings.ncbi_window
            self._open_widget("NCBI Window", ncbi_window.view)
        except Exception as e:
            show_error(self.global_settings, "Error in open_ncbi_window() in main", str(e))

    def toggle_theme(self):
        try:
            self.global_settings.toggle_dark_mode()
            self.view.update_theme_icon()
            self.apply_theme()
        except Exception as e:
            self.logger.error(f"Error toggling theme: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error toggling theme", str(e))

    def toggle_maximize(self):
        if self.view.isMaximized():
            self.view.showNormal()
        else:
            self.view.showMaximized()

    def apply_theme(self):
        # Apply theme to all widgets
        if self.global_settings.is_dark_mode():
            # Apply dark theme
            self.view.setStyleSheet("""
                QWidget { background-color: #2b2b2b; color: #ffffff; }
                QPushButton { background-color: #3a3a3a; border: 1px solid #5a5a5a; }
                QPushButton:hover { background-color: #4a4a4a; }
                QLineEdit, QTextEdit, QPlainTextEdit { background-color: #3a3a3a; border: 1px solid #5a5a5a; }
                QComboBox { background-color: #3a3a3a; border: 1px solid #5a5a5a; }
                QMenuBar { background-color: #2b2b2b; }
                QMenuBar::item:selected { background-color: #3a3a3a; }
                QMenu { background-color: #2b2b2b; }
                QMenu::item:selected { background-color: #3a3a3a; }
            """)
        else:
            # Apply light theme
            self.view.setStyleSheet("""
                QWidget { background-color: #f0f0f0; color: #000000; }
                QPushButton { background-color: #e0e0e0; border: 1px solid #c0c0c0; }
                QPushButton:hover { background-color: #d0d0d0; }
                QLineEdit, QTextEdit, QPlainTextEdit { background-color: #ffffff; border: 1px solid #c0c0c0; }
                QComboBox { background-color: #ffffff; border: 1px solid #c0c0c0; }
                QMenuBar { background-color: #f0f0f0; }
                QMenuBar::item:selected { background-color: #e0e0e0; }
                QMenu { background-color: #f0f0f0; }
                QMenu::item:selected { background-color: #e0e0e0; }
            """)

    # Utility Functions
    def some_long_running_task(self):
        if self.view.progress_bar:
            self.view.reset_progress()
        # ... (perform task steps)
        if self.view.progress_bar:
            self.view.set_progress(50)
        # ... (more task steps)
        if self.view.progress_bar:
            self.view.set_progress(100)

    def show(self):
        try:
            saved_position = self.global_settings.load_window_position("main_window")
            if saved_position:
                self.view.move(saved_position)
            else:
                center_ui(self.view)
            self.view.show()
            self.apply_theme()
        except Exception as e:
            self.global_settings.logger.error(f"Error showing main window: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error showing main window", e)

    def closeEvent(self, event):
        self.global_settings.save_window_position("main_window", self.view.pos())
        event.accept()

    def change_database_directory(self):
        try:
            new_directory = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Select Database Directory", self.global_settings.CSPR_DB,
                QtWidgets.QFileDialog.Option.ShowDirsOnly
            )
            if new_directory:
                if self.global_settings.validate_db_path(new_directory):
                    self.global_settings.save_database_path(new_directory)
                    self.global_settings.initialize_app_directories()
                    self.load_combo_box_data()
                    show_message(12, QtWidgets.QMessageBox.Icon.Information,
                                 "Success", "Database directory changed successfully.")
                else:
                    show_message(12, QtWidgets.QMessageBox.Icon.Warning,
                                 "Invalid Directory", "The selected directory does not contain valid CSPR files.")
        except Exception as e:
            self.logger.error(f"Error changing database directory: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error changing database directory", str(e))
    
    def open_ncbi_website(self):
        ncbi_page()

    def open_repository_website(self):
        repo_page()

    def open_ncbi_blast_website(self):
        ncbi_blast_page()

    def _open_widget(self, title, widget) -> None:
        index = self.view.stacked_widget.indexOf(widget)
        if index == -1:
            self.view.stacked_widget.addWidget(widget)
            index = self.view.stacked_widget.indexOf(widget)
        self.view.stacked_widget.setCurrentIndex(index)
        # Make sure the stacked widget is visible
        self.view.stacked_widget.show()