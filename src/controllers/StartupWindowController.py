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

    def _init_ui(self):
        default_db_path = self.settings.get_db_path()
        self.logger.debug(f"Initial database path: {default_db_path}")
        self.model.set_db_path(default_db_path)
        self.view.set_db_path(default_db_path)
        self._check_cspr_files()

    def _check_cspr_files(self):
        db_path = self.model.get_db_path()
        self.logger.debug(f"Checking CSPR files in: {db_path}")
        if self.settings.validate_db_path(db_path):
            self.logger.debug("Valid DB path found")
            self.view.set_valid_db_state()
        else:
            self.logger.debug("Invalid DB path")
            self.view.set_invalid_db_state()

    def _set_database_directory(self):
        try:
            directory_path = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Open a folder...", self.settings.get_db_path(), QtWidgets.QFileDialog.Option.ShowDirsOnly)

            if not directory_path: 
                self.logger.debug("User cancelled directory selection")
                return

            if not os.path.isdir(directory_path):
                self.logger.debug(f"Selected path is not a directory: {directory_path}")
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Not a directory",
                    message="The directory you selected does not exist.",
                )
                return

            directory_path = self.settings.adjust_path_for_os(directory_path)
            self.logger.debug(f"New directory path selected: {directory_path}")
            self.view.set_db_path(directory_path)
            self.model.set_db_path(directory_path)
            self._check_cspr_files()  # This will update the UI based on the new directory
        except Exception as e:
            self.logger.error(f"Error in _set_database_directory: {str(e)}", exc_info=True)
            show_error(self.settings, "Error in change_directory() in startup window", e)

    def _handle_go_to_home_or_new_genome(self):
        if self.settings.validate_db_path(self.model.get_db_path()):
            self.settings.main_window._switch_to_home_from_startup()
        else:
            self.open_new_genome_tab()

    def open_new_genome_tab(self):
        try:
            new_genome_controller = self.settings.get_new_genome_window()
            self.settings.main_window.open_new_tab("New Genome", new_genome_controller)
        except Exception as e:
            show_error(self.settings, "Error opening New Genome module", str(e))
