import os
import logging
import sys
import platform
from functools import lru_cache
import importlib
from PyQt6.QtCore import QSettings, QObject, pyqtSignal, QThread
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import QApplication
import time

from models.DatabaseManager import DatabaseManager, FileChangeType
from models.ConfigManager import ConfigManager

class ModulePreloader(QThread):
    finished = pyqtSignal(str, object)
    
    def __init__(self, global_settings, module_name):
        super().__init__()
        self.global_settings = global_settings
        self.module_name = module_name
        self.module = None  # Store the loaded module
        
    def run(self):
        try:
            module_path = f"controllers.{self.module_name}Controller"
            if module_path not in self.global_settings._module_cache:
                # Get root directory
                if hasattr(sys, 'frozen'):
                    root_dir = os.path.join(os.path.dirname(sys.executable), 'src')
                    if platform.system() == 'Darwin':
                        root_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(sys.executable))), 
                                              'Contents', 'Resources', 'src')
                else:
                    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

                if root_dir not in sys.path:
                    sys.path.insert(0, root_dir)

                controller_file = os.path.join(root_dir, 'controllers', f"{self.module_name}Controller.py")
                
                if os.path.exists(controller_file):
                    spec = importlib.util.spec_from_file_location(module_path, controller_file)
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                    sys.modules[module_path] = module
                    self.module = module
                    self.finished.emit(self.module_name, module)
                    
        except Exception as e:
            self.global_settings.logger.error(f"Error preloading module {self.module_name}: {str(e)}")

class GlobalSettings(QObject):
    first_time_startup = pyqtSignal()
    endonuclease_updated = pyqtSignal()
    annotation_file_changed = pyqtSignal(str)
    theme_changed = pyqtSignal(str)

    def __init__(self, app_dir_path):
        super().__init__()
        
        self.app_dir_path = app_dir_path
        self.logger = self._setup_logging()
        
        # Initialize important attributes
        self._current_annotation_file = None
        self._module_cache = {}
        self._preloading_modules = {}
        self.main_window = None
        
        # Only preload essential controllers for startup
        self._preload_essential_controllers()
        
        # Start background loading of commonly used modules
        self._background_load_common_modules()
        
        self.config_manager = ConfigManager(app_dir_path=self.app_dir_path, logger=self.logger)
        self.config_manager.load_env()
        
        self.is_first_time_startup = self.config_manager.get_env_value('FIRST_TIME_START', 'TRUE').upper() == 'TRUE'
        
        # Defer database initialization until needed
        self._initialize_directories()
        self._init_db_manager()
        
        # Defer theme initialization
        self._init_theme_settings()

    def _init_db_manager(self):
        """Initialize database manager lazily"""
        if not hasattr(self, 'db_manager'):
            self.db_manager = DatabaseManager(self.logger, self.config_manager)
            self.db_manager.db_files_changed.connect(self._on_db_files_changed)
            self.db_manager.db_validation_changed.connect(self._on_db_validation_changed)
            self.db_manager.db_state_changed.connect(self._on_db_state_changed)
            self.CSPR_DB = self.db_manager.get_db_path()
            self.algorithms = self.config_manager.get_config_value('algorithms', ["Azimuth 2.0"])

    def _init_theme_settings(self):
        """Initialize theme settings lazily"""
        if not hasattr(self, 'settings'):
            self.settings = QSettings("TrinhLab-UTK", "CASPER")
            self.theme = self.settings.value("theme", "light")
            self.light_palette = None
            self.dark_palette = None

    def _preload_essential_controllers(self):
        """Preload only the essential controllers needed for startup"""
        try:
            essential_controllers = [
                "StartupWindow",
                "HomeWindow"
            ] if self.is_first_time_startup else ["HomeWindow"]
            
            for controller_name in essential_controllers:
                self._preload_controller(controller_name)
                
        except Exception as e:
            self.logger.warning(f"Essential controller preloading failed: {str(e)}")

    def _preload_controller(self, window_name):
        """Preload a single controller with optimized imports"""
        try:
            module_path = f"controllers.{window_name}Controller"
            if module_path not in self._module_cache:
                # Get root directory only once
                if not hasattr(self, '_root_dir'):
                    if hasattr(sys, 'frozen'):
                        self._root_dir = os.path.join(os.path.dirname(sys.executable), 'src')
                        if platform.system() == 'Darwin':
                            self._root_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(sys.executable))), 
                                                      'Contents', 'Resources', 'src')
                    else:
                        self._root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

                    if self._root_dir not in sys.path:
                        sys.path.insert(0, self._root_dir)

                controller_file = os.path.join(self._root_dir, 'controllers', f"{window_name}Controller.py")
                
                if os.path.exists(controller_file):
                    spec = importlib.util.spec_from_file_location(module_path, controller_file)
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                    sys.modules[module_path] = module
                    self._module_cache[module_path] = module
                    
        except Exception as e:
            self.logger.warning(f"Failed to preload controller {window_name}: {str(e)}")

    def _on_db_files_changed(self, changes):
        """Handle database file changes"""
        self.logger.debug(f"Database files changed: {changes}")
        # Components should connect directly to db_manager signals
        # This method is for global-level handling if needed

    def _on_db_validation_changed(self, is_valid, message):
        """Handle database validation state changes"""
        self.logger.debug(f"Database validation changed - Valid: {is_valid}, Message: {message}")
        # Handle any global-level validation state changes

    def _on_db_state_changed(self, is_valid, message, changes):
        """Handle combined database state changes"""
        self.logger.debug(f"Database state changed - Valid: {is_valid}, Message: {message}, Changes: {changes}")
        # Handle any global-level state changes

    def get_db_path(self):
        """Get the current database path"""
        return self.db_manager.get_db_path()

    def validate_db_path(self, path):
        """Validate the given database path"""
        return self.db_manager.validate_db_path(path)

    def save_db_path(self, path):
        """Save the database path"""
        return self.db_manager.save_db_path(path)

    def ensure_db_path_exists(self):
        """Ensure the database path exists"""
        self.db_manager.ensure_db_path_exists()

    def update_db_state(self):
        """Check and update the database state"""
        self.db_manager.check_db_state()

    def _initialize_directories(self):
        """Initialize application directories"""
        # app_dir_path is already set in __init__
        
        # Set up source directory paths
        self.src_dir_path = os.path.join(self.app_dir_path, 'src')
        self.ui_dir_path = os.path.join(self.src_dir_path, 'ui')
        self.controllers_dir_path = os.path.join(self.src_dir_path, 'controllers')
        self.models_dir_path = os.path.join(self.src_dir_path, 'models')
        self.views_dir_path = os.path.join(self.src_dir_path, 'views')
        self.utils_dir_path = os.path.join(self.src_dir_path, 'utils')
        self.SeqFinder_dir_path = os.path.join(self.src_dir_path, 'SeqFinder')
        
        # Set up other resource paths
        self.assets_dir_path = os.path.join(self.app_dir_path, 'assets')
        self.config_dir_path = os.path.join(self.app_dir_path, 'config')
        self.casper_info_path = os.path.join(self.config_dir_path, 'CASPERinfo')
        self.off_target_dir_path = os.path.join(self.models_dir_path, 'OffTarget')

        # Ensure critical directories exist
        os.makedirs(self.config_dir_path, exist_ok=True)

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
    
    def get_off_target_dir_path(self):
        return self.off_target_dir_path
    
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
        """Get the controller class with optimized loading"""
        try:
            start_time = time.time()
            
            # Check if module is already cached
            module_path = f"controllers.{window_name}Controller"
            if module_path in self._module_cache:
                controller_module = self._module_cache[module_path]
            else:
                # Fall back to regular import if not cached
                if hasattr(sys, 'frozen'):
                    root_dir = os.path.join(os.path.dirname(sys.executable), 'src')
                    if platform.system() == 'Darwin':
                        root_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(sys.executable))), 
                                              'Contents', 'Resources', 'src')
                else:
                    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

                if root_dir not in sys.path:
                    sys.path.insert(0, root_dir)

                controller_file = os.path.join(root_dir, 'controllers', f"{window_name}Controller.py")
                
                if not os.path.exists(controller_file):
                    raise ImportError(f"Controller file not found: {controller_file}")

                spec = importlib.util.spec_from_file_location(module_path, controller_file)
                controller_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(controller_module)
                sys.modules[module_path] = controller_module
                self._module_cache[module_path] = controller_module

            class_name = f"{window_name}Controller"
            if not hasattr(controller_module, class_name):
                raise AttributeError(f"Controller module does not contain class {class_name}")

            self.logger.debug(f"Window class retrieval took: {time.time() - start_time:.2f} seconds")
            return getattr(controller_module, class_name)

        except Exception as e:
            self.logger.error(f"Failed to load controller {window_name}: {str(e)}")
            raise ImportError(f"Could not load controller for {window_name}") from e

    def _create_window(self, window_name):
        """Create a window instance with dynamic module loading"""
        try:
            WindowClass = self._get_window_class(window_name)
            controller = WindowClass(self)
            
            # Store the reference to prevent garbage collection
            setattr(self, f'_current_{window_name.lower()}_window', controller)
            
            return controller
        except Exception as e:
            self.logger.error(f"Error creating window {window_name}: {str(e)}")
            raise
    
    def get_startup_window(self):
        if not hasattr(self, '_startup_window'):
            from controllers.StartupWindowController import StartupWindowController
            self._startup_window = StartupWindowController(self)
        return self._startup_window
    
    def get_home_window(self):
        """Get or create home window with proper initialization"""
        try:
            controller = self._create_window("HomeWindow")
            self._current_home_window = controller
            
            # Initialize annotation file if needed
            if not hasattr(self, '_current_annotation_file'):
                self._current_annotation_file = None
                if hasattr(controller, 'view'):
                    self._current_annotation_file = controller.view.get_annotation_file()
                
            return controller
        except Exception as e:
            self.logger.error(f"Error creating home window: {str(e)}")
            raise

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

    def _background_load_common_modules(self):
        """Start background loading of commonly used modules"""
        try:
            common_modules = ["MultitargetingWindow", "PopulationAnalysisWindow"]
            for module_name in common_modules:
                if (module_name not in self._module_cache and 
                    module_name not in self._preloading_modules):
                    preloader = ModulePreloader(self, module_name)
                    preloader.finished.connect(self._on_module_preloaded)
                    self._preloading_modules[module_name] = preloader
                    preloader.start()
        except Exception as e:
            self.logger.warning(f"Error starting background module loading: {str(e)}")

    def _on_module_preloaded(self, module_name, module):
        """Handle completion of module preloading"""
        try:
            module_path = f"controllers.{module_name}Controller"
            self._module_cache[module_path] = module
            if module_name in self._preloading_modules:
                preloader = self._preloading_modules[module_name]
                if not preloader.isRunning():  # Only remove if thread is finished
                    del self._preloading_modules[module_name]
            self.logger.debug(f"Module {module_name} preloaded successfully")
        except Exception as e:
            self.logger.error(f"Error handling preloaded module: {str(e)}")

    def get_multitargeting_window(self):
        """Create and return MultitargetingController instance with optimized loading"""
        try:
            start_time = time.time()
            self.logger.debug("Starting multitargeting window creation")
            
            # Check if module is being preloaded
            if "MultitargetingWindow" in self._preloading_modules:
                preloader = self._preloading_modules["MultitargetingWindow"]
                if preloader.isRunning():
                    self.logger.debug("Waiting for preloader to complete...")
                    preloader.wait()
                    if preloader.module:  # Use the stored module
                        WindowClass = getattr(preloader.module, "MultitargetingWindowController")
                    else:
                        WindowClass = self._get_window_class("MultitargetingWindow")
                else:
                    WindowClass = self._get_window_class("MultitargetingWindow")
            else:
                WindowClass = self._get_window_class("MultitargetingWindow")
            
            # Create controller instance
            controller_start = time.time()
            controller = WindowClass(self)
            self.logger.debug(f"Controller instantiation took: {time.time() - controller_start:.2f} seconds")
            
            # Store the reference
            self._current_multitargeting_window = controller
            
            self.logger.debug(f"Total multitargeting window creation took: {time.time() - start_time:.2f} seconds")
            return controller
            
        except Exception as e:
            self.logger.error(f"Error creating multitargeting window: {str(e)}")
            raise

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

    def set_current_annotation_file(self, annotation_file):
        """Set the current annotation file and notify listeners"""
        try:
            if not hasattr(self, '_current_annotation_file'):
                self._current_annotation_file = None
            
            if self._current_annotation_file != annotation_file:
                self._current_annotation_file = annotation_file
                self.logger.debug(f"Current annotation file changed to: {annotation_file}")
                self.annotation_file_changed.emit(annotation_file)
        except Exception as e:
            self.logger.error(f"Error setting current annotation file: {str(e)}")

    def get_current_annotation_file(self):
        """Get the currently selected annotation file"""
        try:
            if not self._current_annotation_file and hasattr(self, '_current_home_window'):
                # Try to get from home window if not set
                home_controller = self._current_home_window
                if hasattr(home_controller, 'view'):
                    self._current_annotation_file = home_controller.view.get_annotation_file()
                    self.logger.debug(f"Got annotation file from home window: {self._current_annotation_file}")
            return self._current_annotation_file
        except Exception as e:
            self.logger.error(f"Error getting current annotation file: {str(e)}")
            return None

    def get_scoring_options_window(self, view_targets_controller):
        """Create and return ScoringOptionsController instance"""
        from controllers.ScoringOptionsController import ScoringOptionsController
        return ScoringOptionsController(self, view_targets_controller)

    def set_theme(self, theme):
        """Set the current theme and notify listeners"""
        self.theme = theme
        self.theme_changed.emit(theme)

    def get_theme(self):
        """Get the current theme"""
        return self.theme

    def get_cotargeting_window(self, view_targets_controller=None):
        """Create and return CoTargetingController instance"""
        if not hasattr(self, '_cotargeting_controller'):
            from controllers.CoTargetingController import CoTargetingController
            self._cotargeting_controller = CoTargetingController(self, view_targets_controller)
        return self._cotargeting_controller

    def get_export_selected_grnas_window(self):
        """Get or create ExportSelectedgRNAs window"""
        if not hasattr(self, '_export_selected_grnas_controller'):
            from controllers.ExportSelectedgRNAsController import ExportSelectedgRNAsController
            self._export_selected_grnas_controller = ExportSelectedgRNAsController(self)
        return self._export_selected_grnas_controller

    def adjust_path_for_os(self, path):
        """
        Adjust file path based on operating system
        """
        try:
            # Convert path separators to match the current OS
            adjusted_path = os.path.normpath(path)
            
            # For Windows, ensure the path uses backslashes
            if platform.system() == 'Windows':
                adjusted_path = adjusted_path.replace('/', '\\')
            # For Unix-like systems (Linux, macOS), ensure the path uses forward slashes
            else:
                adjusted_path = adjusted_path.replace('\\', '/')
            
            self.logger.debug(f"Adjusted path from '{path}' to '{adjusted_path}'")
            return adjusted_path
            
        except Exception as e:
            self.logger.error(f"Error adjusting path: {str(e)}")
            return path  # Return original path if adjustment fails

# Global instance
global_settings = None
