import os
import logging
import sys
import platform
from functools import cached_property, lru_cache
import yaml
from collections import OrderedDict
from dotenv import load_dotenv, set_key
from controllers.NCBIRenameWindowController import NCBIRenameWindowController
from PyQt6.QtCore import QSettings, QObject, QEvent, Qt
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import QApplication
import time
import importlib

def represent_ordereddict(dumper, data):
    return dumper.represent_mapping('tag:yaml.org,2002:map', data.items())

yaml.add_representer(OrderedDict, represent_ordereddict)

class GlobalSettings:
    def __init__(self, application_directory):
        start_time = time.time()
        
        self.app_dir = application_directory
        self.src_dir = os.path.join(self.app_dir, 'src')
        self.config_dir = os.path.join(self.app_dir, 'config')
        self.env_path = os.path.join(self.app_dir, '.env')
        self._load_env()
        self.config = self._initialize_config()  # Initialize config before using it
        
        # Set up logging first
        self.logger = self._setup_logging()
        
        self.CSPR_DB = self.load_database_path()
        self.algorithms = self.config.get('algorithms', ["Azimuth 2.0"])
        self.ui_dir = os.path.join(self.src_dir, self.config['paths']['ui'])
        self.controllers_dir = os.path.join(self.src_dir, self.config['paths']['controllers'])
        self.assets_dir = os.path.join(self.app_dir, self.config['paths']['assets'])
        self.models_dir = os.path.join(self.src_dir, 'models')
        self.views_dir = os.path.join(self.src_dir, 'views')
        self.utils_dir = os.path.join(self.src_dir, 'utils')
        self.SeqFinder_dir = os.path.join(self.src_dir, 'SeqFinder')
        self.casper_info_path = os.path.join(self.config_dir, 'CASPERinfo')

        self.settings = QSettings("TrinhLab-UTK", "CASPER")

        self.logger.info(f"Initialized CSPR_DB: {self.CSPR_DB}")

        self.update_env_file()

        self.theme = self.settings.value("theme", "light")

        self.light_palette = None
        self.dark_palette = None
        self.initialize_palettes()

        self.save_config()

        end_time = time.time()
        initialization_time = end_time - start_time
        self.logger.info(f"GlobalSettings initialization time: {initialization_time:.4f} seconds")

        self.main_window = None  # Initialize as None

    def _load_env(self):
        if not os.path.exists(self.env_path):
            self._create_default_env()
        load_dotenv(self.env_path)

    def _create_default_env(self):
        with open(self.env_path, 'w') as f:
            f.write(f'APP_DIR="{self.app_dir}"\n')
            f.write('CSPR_DB=""\n')  # Empty string, but with quotes
            f.write('LOG_LEVEL=DEBUG\n')

    def _setup_logging(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        
        log_dir = os.path.join(self.app_dir, 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_file_path = os.path.join(log_dir, 'app.log')
        
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

    def load_database_path(self):
        """Load the database path from .env file."""
        db_path = os.getenv('CSPR_DB', '')
        # Remove both single and double quotes if present
        db_path = db_path.strip("'\"")
        if db_path and self.validate_db_path(db_path):
            return db_path
        return ''

    def validate_db_path(self, path):
        """Validate that the given path exists and contains CSPR files."""
        self.logger.debug(f"Validating DB path: {path}")
        if not os.path.isdir(path):
            self.logger.debug(f"Path is not a directory: {path}")
            return False
        has_cspr_files = any(file.endswith(".cspr") for file in os.listdir(path))
        self.logger.debug(f"Path {path} contains CSPR files: {has_cspr_files}")
        return has_cspr_files

    def save_database_path(self, path):
        """Save the database path to .env file."""
        if self.validate_db_path(path):
            self._write_to_env('CSPR_DB', self.adjust_path_for_os(path))
            self.CSPR_DB = self.adjust_path_for_os(path)
            print(f"Saved database path: {self.CSPR_DB}")
        else:
            print(f"Invalid database path: {path}")

    def _write_to_env(self, key, value):
        """Write a key-value pair to the .env file, ensuring proper quoting."""
        with open(self.env_path, 'r') as f:
            lines = f.readlines()

        with open(self.env_path, 'w') as f:
            for line in lines:
                if line.startswith(f'{key}='):
                    if key == 'CSPR_DB' and platform.system() == "Darwin":
                        # Add a trailing slash if it's not already there
                        value = value.rstrip('/') + '/'
                    f.write(f'{key}="{value}"\n')  # Always use double quotes
                else:
                    f.write(line)

    def ensure_correct_env_format(self):
        """Ensure the .env file has the correct format for CSPR_DB."""
        with open(self.env_path, 'r') as f:
            lines = f.readlines()

        with open(self.env_path, 'w') as f:
            for line in lines:
                if line.startswith('CSPR_DB='):
                    value = line.split('=', 1)[1].strip().strip("'\"")
                    f.write(f'CSPR_DB="{value}"\n')
                else:
                    f.write(line)

    def get_app_dir(self):
        return self.app_dir
    
    def get_src_dir(self):
        return self.src_dir
    
    def get_ui_dir(self):
        return self.ui_dir
    
    def get_assets_dir(self):
        return self.assets_dir
    
    def get_controllers_dir(self):
        return self.controllers_dir

    def get_SeqFinder_dir(self):
        return self.SeqFinder_dir

    def get_logger(self):
        return self.logger
    
    def get_casper_info_path(self):
        return self.casper_info_path
    
    def get_db_path(self):
        return self.adjust_path_for_os(self.CSPR_DB)
    
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

    @cached_property
    def main_window(self):
        from controllers.MainWindowController import MainWindowController
        if not hasattr(self, '_main_window'):
            self._main_window = MainWindowController(self)
        return self._main_window
    
    @cached_property
    def home_window(self):
        from controllers.HomeWindowController import HomeWindowController
        if not hasattr(self, '_home_window'):
            self._home_window = HomeWindowController(self)
        return self._home_window

    # @cached_property
    # def startup_window(self):
    #     from views.StartupWindowView_old import StartupWindowView
    #     from controllers.StartupWindowController import StartupWindowController
    #     if not hasattr(self, '_startup_window'):
    #         ui = StartupWindowView(self)
    #         self._startup_window = StartupWindowController(ui, self)
    #     return self._startup_window

    @cached_property
    def multitargeting_window(self):
        from controllers.MultitargetingWindowController import MultitargetingWindowController
        if not hasattr(self, '_multitargeting_window'):
            self._multitargeting_window = MultitargetingWindowController(self)
        return self._multitargeting_window

    @cached_property
    def population_analysis_window(self):
        from controllers.PopulationAnalysisWindowController import PopulationAnalysisWindowController
        if not hasattr(self, '_population_analysis_window'):
            self._population_analysis_window = PopulationAnalysisWindowController(self)
        return self._population_analysis_window
    
    @cached_property
    def new_endonuclease_window(self):
        from controllers.NewEndonucleaseController import NewEndonucleaseController
        if not hasattr(self, '_new_endonuclease_window'):
            self._new_endonuclease_window = NewEndonucleaseController(self)
        return self._new_endonuclease_window
    
    @cached_property
    def ncbi_window(self):
        from controllers.NCBIWindowController import NCBIWindowController
        if not hasattr(self, '_ncbi_window'):
            self._ncbi_window = NCBIWindowController(self)
        return self._ncbi_window

    @cached_property
    def ncbi_rename_window(self):
        if not hasattr(self, '_ncbi_rename_window'):
            self._ncbi_rename_window = NCBIRenameWindowController(self, [], None)
        return self._ncbi_rename_window

    def get_default_database_path(self):
        default_path = os.path.join(self.app_dir, 'CSPR_DB')
        return self.adjust_path_for_os(default_path)

    def set_db_path(self, directory_path):
        self.CSPR_DB = self.adjust_path_for_os(directory_path)
        self.save_database_path(self.CSPR_DB)

    def adjust_path_for_os(self, path):
        """Adjust the file path based on the operating system."""
        if platform.system() == "Windows":
            return path.replace("/", "\\")
        else:
            path = path.replace("\\", "/")
            if not path.endswith("/"):
                path += "/"
        return path

    def initialize_app_directories(self):
        """Create required directories if they don't exist."""
        required_dirs = ["FNA", "GBFF"]
        for directory in required_dirs:
            path = os.path.join(self.CSPR_DB, directory)
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)
                self.logger.debug(f"Directory created: {path}")

    def _initialize_config(self):
        config_path = os.path.join(self.config_dir, 'config.yml')
        if not os.path.exists(config_path):
            print(f"Config file not found. Creating default config at {config_path}")
            self._create_default_config(config_path)
        return self._load_config(config_path)

    def _create_default_config(self, config_path):
        config = OrderedDict([
            ('paths', OrderedDict([
                ('ui', 'ui'),
                ('controllers', 'controllers'),
                ('assets', 'assets'),
                ('casper_info', 'CASPERinfo')
            ])),
            ('algorithms', ['Azimuth 2.0'])
        ])
        
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        print(f"Default config created at {config_path}")

    def _load_config(self, config_path):
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            print(f"Error loading config from {config_path}: {e}")
            return None

    def save_config(self):
        config_path = os.path.join(self.config_dir, 'config.yml')
        config = OrderedDict([
            ('paths', OrderedDict([
                ('ui', self.config['paths']['ui']),
                ('controllers', self.config['paths']['controllers']),
                ('assets', self.config['paths']['assets']),
                ('casper_info', self.config['paths']['casper_info'])
            ])),
            ('algorithms', self.algorithms)
        ])
        
        with open(config_path, 'w') as f:
            yaml.dump({'paths': config['paths']}, f, default_flow_style=False)
            f.write("\n")  # Add an extra newline here
            yaml.dump({'algorithms': config['algorithms']}, f, default_flow_style=False)
        print(f"Config saved to {config_path}")

    def save_window_position(self, window_name, position):
        self.settings.setValue(f"{window_name}_pos", position)

    def load_window_position(self, window_name):
        return self.settings.value(f"{window_name}_pos")

    def save_window_size(self, window_name, size):
        self.settings.setValue(f"{window_name}_size", size)

    def load_window_size(self, window_name):
        return self.settings.value(f"{window_name}_size")

    def update_env_file(self):
        """Update the .env file to ensure correct formatting."""
        with open(self.env_path, 'r') as f:
            lines = f.readlines()

        with open(self.env_path, 'w') as f:
            for line in lines:
                if line.startswith('CSPR_DB='):
                    value = line.split('=', 1)[1].strip().strip("'\"")
                    f.write(f'CSPR_DB="{value}"\n')
                else:
                    f.write(line)


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


    @cached_property
    def new_endonuclease_controller(self):
        from controllers.NewEndonucleaseController import NewEndonucleaseController
        from views.NewEndonucleaseView import NewEndonucleaseView
        from models.NewEndonucleaseModel import NewEndonucleaseModel
        
        view = NewEndonucleaseView(self)
        model = NewEndonucleaseModel(self)
        return NewEndonucleaseController(self, view, model)

    @cached_property
    def population_analysis_window(self):
        from controllers.PopulationAnalysisWindowController import PopulationAnalysisWindowController
        if not hasattr(self, '_population_analysis_window') or self._population_analysis_window is None:
            self._population_analysis_window = PopulationAnalysisWindowController(self)
        return self._population_analysis_window

    def _initialize_new_genome_window(self):
        from controllers.NewGenomeWindowController import NewGenomeWindowController
        return NewGenomeWindowController(self)

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
        self._current_home_window = controller  # Store reference
        return controller

    def get_new_genome_window(self):
        controller = self._create_window("NewGenomeWindow")
        self._current_new_genome_window = controller  # Store reference
        return controller

    def get_new_endonuclease_window(self):
        controller = self._create_window("NewEndonuclease")
        self._current_new_endonuclease_window = controller  # Store reference
        return controller

    def get_ncbi_window(self):
        controller = self._create_window("NCBIWindow")
        self._current_ncbi_window = controller  # Store reference
        return controller

    def get_multitargeting_window(self):
        controller = self._create_window("MultitargetingWindow")
        self._current_multitargeting_window = controller  # Store reference
        return controller

    def get_population_analysis_window(self):
        controller = self._create_window("PopulationAnalysisWindow")
        self._current_population_analysis_window = controller  # Store reference
        return controller

    def set_main_window(self, main_window):
        self.main_window = main_window

# Global instance
global_settings = None