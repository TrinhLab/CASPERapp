from PyQt6.QtCore import QObject, pyqtSignal

class StartupWindowModel(QObject):
    db_state_updated = pyqtSignal(bool, str, list)

    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        
        # Connect to the combined signal from GlobalSettings
        self.settings.db_state_updated.connect(self.on_db_state_updated)

    def get_db_path(self):
        return self.settings.get_db_path()

    def save_db_path(self, directory_path):
        success, message = self.settings.save_db_path(directory_path)
        # Note: The actual db_state_updated signal will be emitted by the DatabaseManager

    def on_db_state_updated(self, is_valid, message, cspr_files):
        self.logger.debug(f"StartupWindowModel received db state update: valid={is_valid}, message={message}, cspr_files_count={len(cspr_files)}")
        self.db_state_updated.emit(is_valid, message, cspr_files)
