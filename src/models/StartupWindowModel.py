from PyQt6.QtCore import QObject, pyqtSignal

class StartupWindowModel(QObject):
    db_state_updated = pyqtSignal(bool, str, list)

    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        
        # Connect to the DatabaseManager's signals instead of GlobalSettings
        self.settings.db_manager.db_state_changed.connect(self.on_db_state_updated)

    def get_db_path(self):
        return self.settings.get_db_path()

    def save_db_path(self, directory_path):
        """Save the database path and trigger validation"""
        self.logger.debug(f"Saving database path: {directory_path}")
        # The db_manager will emit its own signals that we're now listening to
        success, message = self.settings.save_db_path(directory_path)
        return success, message

    def on_db_state_updated(self, is_valid, message, cspr_files):
        """Handle database state updates"""
        self.logger.debug(f"StartupWindowModel received db state update: valid={is_valid}, message={message}, cspr_files_count={len(cspr_files)}")
        self.db_state_updated.emit(is_valid, message, cspr_files)
