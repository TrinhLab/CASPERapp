import sys
import os
import platform
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt, QCoreApplication
from models.GlobalSettings import GlobalSettings
from utils.ui import show_error
import importlib
import concurrent.futures

def preload_modules():
    """Preload commonly used modules in parallel"""
    modules_to_load = [
        'PyQt6.QtWidgets',
        'PyQt6.QtCore',
        'PyQt6.QtGui',
        'controllers.MainWindowController',
        'views.MainWindowView',
        'models.MainWindowModel',
        'models.DatabaseManager',
        'models.ConfigManager'
    ]
    
    def import_module(module_name):
        try:
            importlib.import_module(module_name)
            return True, module_name
        except Exception as e:
            return False, f"Failed to load {module_name}: {str(e)}"
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(import_module, modules_to_load)

def get_app_directory():
    """Determine the application root directory based on whether we're frozen or not"""
    if hasattr(sys, 'frozen'):
        if platform.system() == 'Darwin':  # macOS
            bundle_dir = os.path.abspath(os.path.dirname(sys.executable))
            return os.path.join(os.path.dirname(os.path.dirname(bundle_dir)), 'Contents', 'Resources')
        else:
            return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    RESTART_CODE = 1000
    
    # Preload modules in parallel
    preload_modules()
    
    while True:
        app = QApplication(sys.argv)
        app.setOrganizationName("TrinhLab-UTK")
        app.setApplicationName("CASPER")

        # Enable high DPI scaling
        if hasattr(Qt.ApplicationAttribute, 'AA_UseHighDpiPixmaps'):
            app.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps)
            app.setAttribute(Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)

        # Enable Qt's built-in caching mechanisms
        QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
        
        # Get the application directory
        app_dir_path = get_app_directory()

        try:
            global_settings = GlobalSettings(app_dir_path)
            
            # Import here after preloading
            from controllers.MainWindowController import MainWindowController
            main_window_controller = MainWindowController(global_settings)
            global_settings.set_main_window(main_window_controller)
            main_window_controller.show()

            exit_code = app.exec()

            if exit_code != RESTART_CODE:
                sys.exit(exit_code)
                break

            main_window_controller = None
            global_settings = None
            app = None
        except Exception as e:
            show_error(global_settings, "An error occurred during application initialization", e)

if __name__ == '__main__':
    main()
