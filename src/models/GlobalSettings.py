import os
import logging
import sys
import platform
from functools import lru_cache
import importlib
from PyQt6.QtCore import QSettings, QObject, pyqtSignal
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import QApplication

from models.DatabaseManager import DatabaseManager
from models.ConfigManager import ConfigManager

class GlobalSettings(QObject):
    db_state_updated = pyqtSignal(bool, str, list)  # Combined signal
    first_time_startup = pyqtSignal()  # New signal
    endonuclease_updated = pyqtSignal()

    def __init__(self, app_dir_path):
        super().__init__()
        
        self.app_dir_path = app_dir_path
        
        self.logger = self._setup_logging()
        
        self.config_manager = ConfigManager(app_dir_path=self.app_dir_path, logger=self.logger)
        self.config_manager.load_env()
        
        self.is_first_time_startup = self.config_manager.get_env_value('FIRST_TIME_START', 'TRUE').upper() == 'TRUE'
        
        self._initialize_directories()  # Add this line
        
        self.db_manager = DatabaseManager(self.logger, self.config_manager)
        self.db_manager.db_state_updated.connect(self._on_db_state_updated)
        
        self.CSPR_DB = self.db_manager.get_db_path()
        self.algorithms = self.config_manager.get_config_value('algorithms', ["Azimuth 2.0"])

        self.settings = QSettings("TrinhLab-UTK", "CASPER")
        self.theme = self.settings.value("theme", "light")

        self.light_palette = None
        self.dark_palette = None
        self.initialize_palettes()

        self.main_window = None 

    def _initialize_directories(self):
        self.src_dir_path = os.path.join(self.app_dir_path, 'src')
        self.ui_dir_path = os.path.join(self.src_dir_path, self.config_manager.get_config_value('paths.ui'))
        self.controllers_dir_path = os.path.join(self.src_dir_path, self.config_manager.get_config_value('paths.controllers'))
        self.assets_dir_path = os.path.join(self.app_dir_path, self.config_manager.get_config_value('paths.assets'))
        self.models_dir_path = os.path.join(self.src_dir_path, 'models')
        self.views_dir_path = os.path.join(self.src_dir_path, 'views')
        self.utils_dir_path = os.path.join(self.src_dir_path, 'utils')
        self.SeqFinder_dir_path = os.path.join(self.src_dir_path, 'SeqFinder')
        self.casper_info_path = os.path.join(self.config_manager.config_dir_path, 'CASPERinfo')

    def _setup_logging(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        
        log_dir_path = os.path.join(self.app_dir_path, 'logs')
        os.makedirs(log_dir_path, exist_ok=True)
        log_file_path = os.path.join(log_dir_path, 'app.log')
        
        fh = logging.FileHandler(log_file_path, mode='w')
        fh_formatter = logging.Formatter('%(asctime)s %(levelname)s %(lineno)d:%(filename)s(%(process)d) - %(message)s')
        fh.setFormatter(fh_formatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

        logger.info(f"System OS: {platform.system()}")
        if hasattr(sys, 'frozen'):
            logger.info("Running a packaged version of CASPER.")
        else:
            logger.info("Running a non-packaged version of CASPER.")

        return logger

    def get_db_path(self):
        return self.db_manager.get_db_path()

    def validate_db_path(self, path):
        return self.db_manager.validate_db_path(path)

    def save_db_path(self, path):
        return self.db_manager.save_db_path(path)

    def _on_db_state_updated(self, is_valid, message, cspr_files):
        self.db_state_updated.emit(is_valid, message, cspr_files)

    def ensure_db_path_exists(self):
        self.db_manager.ensure_db_path_exists()

    def adjust_path_for_os(self, path):
        return self.db_manager.adjust_path_for_os(path)

    def get_app_dir_path(self):
        return self.app_dir_path
    
    def get_src_dir_path(self):
        return self.src_dir_path
    
    def get_ui_dir_path(self):
        return self.ui_dir_path
    
    def get_assets_dir_path(self):
        return self.assets_dir_path
    
    def get_controllers_dir_path(self):
        return self.controllers_dir_path

    def get_SeqFinder_dir_path(self):
        return self.SeqFinder_dir_path

    def get_logger(self):
        return self.logger
    
    def get_casper_info_path(self):
        return self.casper_info_path
    
    def get_theme(self):
        return self.theme
    
    def set_theme(self, theme):
        self.theme = theme
        self.settings.setValue("theme", theme)
        self.apply_theme()

    def apply_theme(self):
        app = QApplication.instance()
        if self.theme == "dark":
            app.setPalette(self.dark_palette)
        else:
            app.setPalette(self.light_palette)

    def initialize_palettes(self):
        self.light_palette = QPalette()  # Use default Qt light palette
        self.dark_palette = QPalette()

        # Set up dark palette
        self.dark_palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
        self.dark_palette.setColor(QPalette.ColorRole.WindowText, QColor(255, 255, 255))
        self.dark_palette.setColor(QPalette.ColorRole.Base, QColor(25, 25, 25))
        self.dark_palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
        self.dark_palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(255, 255, 255))
        self.dark_palette.setColor(QPalette.ColorRole.ToolTipText, QColor(255, 255, 255))
        self.dark_palette.setColor(QPalette.ColorRole.Text, QColor(255, 255, 255))
        self.dark_palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
        self.dark_palette.setColor(QPalette.ColorRole.ButtonText, QColor(255, 255, 255))
        self.dark_palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 0, 0))
        self.dark_palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
        self.dark_palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
        self.dark_palette.setColor(QPalette.ColorRole.HighlightedText, QColor(0, 0, 0))

    def save_config(self):
        self.config_manager.save_config()

    def get_config_value(self, key, default=None):
        return self.config_manager.get_config_value(key, default)

    def set_config_value(self, key, value):
        self.config_manager.set_config_value(key, value)

    def save_window_position(self, window_name, position):
        self.settings.setValue(f"{window_name}_pos", position)

    def load_window_position(self, window_name):
        return self.settings.value(f"{window_name}_pos")

    def save_window_size(self, window_name, size):
        self.settings.setValue(f"{window_name}_size", size)

    def load_window_size(self, window_name):
        return self.settings.value(f"{window_name}_size")

    @lru_cache(maxsize=None)
    def _get_window_class(self, window_name):
        module_name = f"controllers.{window_name}Controller"
        class_name = f"{window_name}Controller"
        module = importlib.import_module(module_name)
        return getattr(module, class_name)

    def _create_window(self, window_name):
        WindowClass = self._get_window_class(window_name)
        controller = WindowClass(self)
        return controller
    
    def get_startup_window(self):
        if not hasattr(self, '_startup_window'):
            from controllers.StartupWindowController import StartupWindowController
            self._startup_window = StartupWindowController(self)
        return self._startup_window
    
    def get_home_window(self):
        controller = self._create_window("HomeWindow")
        self._current_home_window = controller
        return controller

    def get_new_genome_window(self):
        controller = self._create_window("NewGenomeWindow")
        self._current_new_genome_window = controller
        return controller

    def get_new_endonuclease_window(self):
        controller = self._create_window("NewEndonuclease")
        self._current_new_endonuclease_window = controller
        return controller

    def get_ncbi_window(self):
        controller = self._create_window("NCBIWindow")
        self._current_ncbi_window = controller
        return controller

    def get_multitargeting_window(self):
        controller = self._create_window("MultitargetingWindow")
        self._current_multitargeting_window = controller
        return controller

    def get_population_analysis_window(self):
        controller = self._create_window("PopulationAnalysisWindow")
        self._current_population_analysis_window = controller
        return controller

    def get_find_targets_window(self):
        controller = self._create_window("FindTargets")
        self._current_find_targets_window = controller
        return controller
    
    def get_view_targets_window(self):
        controller = self._create_window("ViewTargets")
        self._current_view_targets_window = controller  
        return controller

    def set_main_window(self, main_window):
        self.main_window = main_window

    def update_db_state(self):
        self.db_manager.check_db_state()

    def _on_env_file_created(self):
        self.logger.info("GlobalSettings: _on_env_file_created")
        self.is_first_time_startup = True
        self.first_time_startup.emit()

    def set_first_time_startup_completed(self):
        self.logger.info("First time startup completed")
        self.config_manager.set_env_value('FIRST_TIME_START', 'FALSE')
        self.is_first_time_startup = False

    def check_and_emit_first_time_startup(self):
        if self.is_first_time_startup:
            self.logger.info("Emitting first_time_startup signal")
            self.first_time_startup.emit()
            # We no longer set FIRST_TIME_START to FALSE here
    
    def get_organism_to_endonuclease(self):
        if hasattr(self, '_current_home_window'):
            return self._current_home_window.model.get_organism_to_endonuclease()
        else:
            self.logger.warning("Home window not initialized when trying to get organism_to_endonuclease")
            return {}

    def get_annotation_files(self):
        return self.home_window_model.get_annotation_files()
    
    def get_endonucleases(self):
        return self.config_manager.get_endonucleases()

# Global instance
global_settings = None
