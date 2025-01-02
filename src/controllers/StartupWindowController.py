import os
from PyQt6 import QtWidgets 
from models.StartupWindowModel import StartupWindowModel
from utils.ui import show_message, show_error
from views.StartupWindowView import StartupWindowView
import sys

class StartupWindowController:
    def __init__(self, global_settings, keep_db_path=False):
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self.is_active = True
        self.keep_db_path = keep_db_path

        try:
            self.view = StartupWindowView(self.settings)
            self.model = StartupWindowModel(self.settings)

            self._setup_connections()
            self._init_ui()

            self.view.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Fixed,
                QtWidgets.QSizePolicy.Policy.Fixed
            )
        except Exception as e:
            show_error(self.settings, "Error initializing StartupWindowController", str(e))

    def _setup_connections(self):
        self.view.push_button_change_directory.clicked.connect(self._set_database_directory)
        self.view.push_button_go_to_home_or_new_genome.clicked.connect(self._handle_go_to_home_or_new_genome)
        self.view.db_path_text_changed.connect(self._on_db_path_text_changed)
        self.model.db_state_updated.connect(self._on_db_state_updated)
        self.settings.db_manager.db_validation_changed.connect(self._on_db_validation_changed)
        self.settings.db_manager.db_files_changed.connect(self._on_db_files_changed)
        self.view.open_new_genome_requested.connect(self.open_new_genome_tab)

    def _on_db_path_text_changed(self, new_path):
        """Handle database path text changes"""
        try:
            # Avoid triggering save if the path hasn't actually changed
            if new_path == self.model.get_db_path():
                return
            
            self.model.save_db_path(new_path)
        except Exception as e:
            self.logger.error(f"Error handling db path change: {str(e)}")

    def _on_db_state_updated(self, is_valid, message, cspr_files):
        if self.is_active and hasattr(self, 'view'):
            self.view.set_db_status(is_valid, message)

    def _on_db_validation_changed(self, is_valid, message):
        """Handle database validation state changes"""
        if self.is_active and hasattr(self, 'view'):
            self.view.set_db_status(is_valid, message)

    def _on_db_files_changed(self, changes):
        """Handle database file changes"""
        if self.is_active and hasattr(self, 'view'):
            # Re-validate the current path
            db_path = self.model.get_db_path()
            is_valid, message = self.settings.validate_db_path(db_path)
            self.view.set_db_status(is_valid, message)
            
            # If path is now valid and we have CSPR files, update button text
            if is_valid:
                self.view.push_button_go_to_home_or_new_genome.setText("Go to Home")

    def _init_ui(self):
        """Initialize the UI with the correct database path"""
        try:
            # Always get the current path from settings
            db_path = self.model.get_db_path()
            
            # For true first time startup (no previous path), show default path
            if self.settings.is_first_time_startup and not db_path:
                db_path = self.settings.db_manager.get_default_database_path()
            # For invalid path case, keep the existing path
            elif not self.keep_db_path and not self.settings.is_first_time_startup:
                db_path = ''
                
            self.logger.debug(f"Initial database path: {db_path}, keep_db_path: {self.keep_db_path}, "
                             f"first_time_startup: {self.settings.is_first_time_startup}")
            
            # Set the path in the view first
            self.view.set_db_path(db_path)
            
            # Then initialize the state
            if db_path:
                is_valid, message = self.settings.validate_db_path(db_path)
                self.view.set_db_status(is_valid, message)
            else:
                self.view.set_db_status(False, "No directory selected")
                
        except Exception as e:
            self.logger.error(f"Error in _init_ui: {str(e)}")
            raise

    def _init_db_state(self, db_path):
        """Initialize database state"""
        try:
            # Only update the view's path if it's different from current
            current_view_path = self.view.get_db_path()
            if current_view_path != db_path:
                self.view.set_db_path(db_path)
                
            # Validate and update status
            is_valid, message = self.settings.validate_db_path(db_path)
            self.view.set_db_status(is_valid, message)
            
        except Exception as e:
            self.logger.error(f"Error in _init_db_state: {str(e)}")
            raise

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
            self.settings.set_first_time_startup_completed()
            self.restart_application()
        else:
            self.logger.warning(f"Invalid database path: {message}")
            self.open_new_genome_tab()

    def restart_application(self):
        try:
            self.logger.info("Restarting application...")
            # Get the current application instance
            app = QtWidgets.QApplication.instance()
            app.exit(1000)  
        except Exception as e:
            self.logger.error(f"Error restarting application: {str(e)}", exc_info=True)
            show_error(self.settings, "Error restarting application", str(e))

    def open_new_genome_tab(self):
        try:
            self.logger.debug("Opening New Genome tab")
            if hasattr(self.settings, 'main_window'):
                self.settings.main_window.open_new_genome_tab()
        except Exception as e:
            self.logger.error(f"Error opening New Genome tab: {str(e)}", exc_info=True)
            show_error(self.settings, "Error opening New Genome module", str(e))

    # Add this method
    def deactivate(self):
        self.is_active = False
        if hasattr(self, 'view'):
            delattr(self, 'view')
