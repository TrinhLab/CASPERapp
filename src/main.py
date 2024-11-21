import sys
import os
import platform
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt
from models.GlobalSettings import GlobalSettings
from utils.ui import show_error

def get_app_directory():
    """Determine the application root directory based on whether we're frozen or not"""
    if hasattr(sys, 'frozen'):
        if platform.system() == 'Darwin':  # macOS
            # Get the path to the executable inside the .app bundle
            bundle_dir = os.path.abspath(os.path.dirname(sys.executable))
            # Navigate up to Contents directory and set Resources as app_dir
            return os.path.join(os.path.dirname(os.path.dirname(bundle_dir)), 'Contents', 'Resources')
        else:
            # For other platforms when frozen
            return os.path.dirname(sys.executable)
    else:
        # Development environment - go up one level from src directory
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    RESTART_CODE = 1000  # Define restart code constant
    
    while True:
        app = QApplication(sys.argv)
        app.setOrganizationName("TrinhLab-UTK")
        app.setApplicationName("CASPER")

        if hasattr(Qt.ApplicationAttribute, 'AA_UseHighDpiPixmaps'):
            app.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps)

        # Get the application directory
        app_dir_path = get_app_directory()

        try:
            global_settings = GlobalSettings(app_dir_path)
            
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
