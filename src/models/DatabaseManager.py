import os
from PyQt6.QtCore import QObject, pyqtSignal, QFileSystemWatcher

class DatabaseManager(QObject):
    db_state_updated = pyqtSignal(bool, str, list)  # Combined signal

    def __init__(self, logger, config_manager):
        super().__init__()
        self.logger = logger
        self.config_manager = config_manager
        self.db_path = None
        self.file_watcher = QFileSystemWatcher()  # Initialize file_watcher here
        self.file_watcher.directoryChanged.connect(self._on_directory_changed)
        self.load_database_path()
        self._update_watched_directory()

    def load_database_path(self):
        """Load the database path from .env file or set default if empty."""
        db_path = self.config_manager.get_env_value('CSPR_DB', '')
        # Remove both single and double quotes if present
        db_path = db_path.strip("'\"")
        if not db_path:
            db_path = self.get_default_database_path()
            self.save_db_path(db_path)
        self.db_path = db_path
        return self.db_path

    def validate_db_path(self, path):
        """Validate that the given path exists and contains CSPR files."""
        self.logger.debug(f"Validating DB path: {path}")
        if not os.path.isdir(path):
            self.logger.debug(f"Path is not a directory: {path}")
            return False, "The selected path is not a directory."
        has_cspr_files = any(file.endswith(".cspr") for file in os.listdir(path))
        if not has_cspr_files:
            self.logger.debug(f"Path {path} does not contain CSPR files")
            return False, "The selected directory does not contain any CSPR files."
        self.logger.debug(f"Path {path} is valid and contains CSPR files")
        return True, "Valid database path selected."

    def save_db_path(self, path):
        """Set and save the database path."""
        if not path:
            self.logger.warning("Attempting to save an empty database path")
            return False, "Empty database path is not allowed."

        # Ensure the path is a string and properly quoted
        path = str(path).strip("'\"")

        # Validate the database path
        is_valid, message = self.validate_db_path(path)
        if not is_valid:
            self.logger.warning(f"Invalid database path: {path}")
            self.db_state_updated.emit(False, message, [])
            self.db_path = path
            self.config_manager.set_env_value('CSPR_DB', path)
            self._update_watched_directory()
            return False, message

        # Set the db_path attribute
        self.db_path = path

        try:
            self.config_manager.set_env_value('CSPR_DB', path)
            self.logger.info(f"Database path set and saved: {path}")
            self.db_state_updated.emit(True, "Database path saved successfully.", [])
            self._update_watched_directory()
            return True, "Database path saved successfully."
        except Exception as e:
            error_message = f"Error saving database path: {str(e)}"
            self.logger.error(error_message)
            self.db_state_updated.emit(False, error_message, [])
            return False, error_message

    def get_db_path(self):
        return self.db_path

    def ensure_db_path_exists(self):
        """Ensure that the database path exists, creating it if necessary."""
        if not os.path.exists(self.db_path):
            try:
                os.makedirs(self.db_path)
                self.logger.info(f"Created database directory: {self.db_path}")
            except Exception as e:
                self.logger.error(f"Failed to create database directory: {self.db_path}. Error: {str(e)}")
                raise

    def get_default_database_path(self):
        """Get the default database path."""
        documents_path = os.path.expanduser("~/Documents")
        default_path = os.path.join(documents_path, 'CASPERdb')
        return self.adjust_path_for_os(default_path)

    def adjust_path_for_os(self, path):
        """Adjust the file path based on the operating system."""
        if os.name == "nt":  # Windows
            adjusted_path = path.replace("/", "\\")
        else:  # Unix-like systems
            adjusted_path = path.replace("\\", "/")
        
        # Add trailing slash for directories on non-Windows systems
        # if os.name != "nt" and os.path.isdir(adjusted_path):
        #     adjusted_path = os.path.join(adjusted_path, '')
        
        return adjusted_path

    def _update_watched_directory(self):
        """Update the directory being watched by QFileSystemWatcher."""
        if self.file_watcher.directories():
            self.file_watcher.removePaths(self.file_watcher.directories())
        if self.db_path and os.path.isdir(self.db_path):
            self.file_watcher.addPath(self.db_path)
            self.logger.debug(f"Now watching directory: {self.db_path}")

    def _on_directory_changed(self, path):
        """Handle changes in the watched directory."""
        self.logger.debug(f"Detected change in directory: {path}")
        cspr_files = self._get_cspr_files()
        is_valid, message = self.validate_db_path(path)
        self.db_state_updated.emit(is_valid, message, cspr_files)

    def _get_cspr_files(self):
        """Get a list of CSPR files in the current database directory."""
        if not self.db_path or not os.path.isdir(self.db_path):
            return []
        return [f for f in os.listdir(self.db_path) if f.endswith('.cspr')]

    def check_db_state(self):
        """Check the current state of the database and emit signals if changed."""
        self.logger.debug("Checking database state")
        if not self.db_path:
            self.load_database_path()

        is_valid, message = self.validate_db_path(self.db_path)
        self.logger.debug(f"Database state: valid={is_valid}, message={message}")
        
        cspr_files = self._get_cspr_files()
        self.logger.debug(f"Found {len(cspr_files)} CSPR files in the database")
        
        message = f"Database is valid. Contains {len(cspr_files)} CSPR files."
        
        self.db_state_updated.emit(is_valid, message, cspr_files)
