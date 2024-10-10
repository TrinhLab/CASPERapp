class StartupWindowModel:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = global_settings.get_logger()

    def get_db_path(self):
        return self.settings.get_db_path()

    def set_db_path(self, directory_path):
        print(f"Setting database directory to: {directory_path}")
        self.settings.set_db_path(directory_path)

