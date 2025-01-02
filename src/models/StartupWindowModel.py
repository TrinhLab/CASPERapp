from PyQt6.QtCore import QObject, pyqtSignal

class StartupWindowModel(QObject):
    db_state_updated = pyqtSignal(bool, str, list)
    _is_saving = False  # Add flag to prevent recursion

    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        
        # Connect to the DatabaseManager's signals instead of GlobalSettings
        self.settings.db_manager.db_state_changed.connect(self.on_db_state_updated)

    def get_db_path(self):
        """Get the current database path without modifying it"""
        return self.settings.get_db_path()

    def save_db_path(self, path):
        """Save the database path"""
        try:
            if self._is_saving:  # Prevent recursive saves
                return
                
            self._is_saving = True
            try:
                # Don't clear the path if it's invalid - let the controller handle that
                self.settings.save_db_path(path)
                self.settings.update_db_state()
            finally:
                self._is_saving = False
                
        except Exception as e:
            self._is_saving = False
            self.logger.error(f"Error saving database path: {str(e)}")
            raise

    def on_db_state_updated(self, is_valid, message, cspr_files):
        """Handle database state updates"""
        if self._is_saving:  # Don't emit signals during save operation
            return
            
        self.logger.debug(f"StartupWindowModel received db state update: valid={is_valid}, message={message}")
        self.db_state_updated.emit(is_valid, message, cspr_files)
