import os
from PyQt6 import QtWidgets 
from models.StartupWindowModel import StartupWindowModel
from utils.ui import show_message, show_error
from views.StartupWindowView import StartupWindowView

class StartupWindowController:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = self.settings.get_logger()

        try:
            self.view = StartupWindowView(self.settings)
            self.model = StartupWindowModel(self.settings)

            self._setup_connections()
            self._init_ui()
        except Exception as e:
            show_error(self.settings, "Error initializing StartupWindowController", str(e))

    def _setup_connections(self):
        self.view.push_button_change_directory.clicked.connect(self._set_database_directory)
        self.view.push_button_go_to_home_or_new_genome.clicked.connect(self._handle_go_to_home_or_new_genome)
        self.view.db_path_text_changed.connect(self._on_db_path_text_changed)
        self.model.db_state_updated.connect(self._on_db_state_updated)
        self.view.open_new_genome_requested.connect(self.open_new_genome_tab)

    def _on_db_path_text_changed(self, new_path):
        self.model.save_db_path(new_path)

    def _on_db_state_updated(self, is_valid, message, cspr_files):
        self.view.set_db_status(is_valid, message)

    def _init_ui(self):
        db_path = self.model.get_db_path()
        self.logger.debug(f"Initial database path: {db_path}")
        self._init_db_state(db_path)

    def _init_db_state(self, db_path):
        self.view.set_db_path(db_path)
        is_valid, message = self.settings.validate_db_path(db_path)
        self.view.set_db_status(is_valid, message)

    def _set_database_directory(self):
        try:
            directory_path = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Open a folder...", self.settings.get_db_path(), QtWidgets.QFileDialog.Option.ShowDirsOnly)

            if directory_path:
                directory_path = self.settings.adjust_path_for_os(directory_path)
                self.logger.debug(f"New directory path selected: {directory_path}")
                self.model.save_db_path(directory_path)
                self.view.set_db_path(directory_path)
        except Exception as e:
            self.logger.error(f"Error in _set_database_directory: {str(e)}", exc_info=True)
            show_error(self.settings, "Error in change_directory() in startup window", e)

    def _handle_go_to_home_or_new_genome(self):
        self.logger.debug(f"Handle go to home or new genome: {self.model.get_db_path()}")
        is_valid, message = self.settings.validate_db_path(self.model.get_db_path())
        if is_valid:
            self.settings.set_first_time_startup_completed()  # New method call
            self.settings.main_window._switch_to_home_from_startup()
        else:
            self.logger.warning(f"Invalid database path: {message}")
            self.open_new_genome_tab()

    def open_new_genome_tab(self):
        try:
            self.logger.debug("Opening New Genome tab")
            self.settings.main_window.open_new_genome_tab()
        except Exception as e:
            self.logger.error(f"Error opening New Genome tab: {str(e)}", exc_info=True)
            show_error(self.settings, "Error opening New Genome module", str(e))
