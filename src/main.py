import sys
import os
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt
from models.GlobalSettings import GlobalSettings
from utils.ui import show_error
from controllers.MainWindowController import MainWindowController

def main():
    app = QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")

    if hasattr(Qt.ApplicationAttribute, 'AA_UseHighDpiPixmaps'):
        app.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps)

    app_dir_path = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))

    global_settings = GlobalSettings(app_dir_path)
    global_settings.apply_theme()

    try:
        main_window_controller = MainWindowController(global_settings)
        global_settings.set_main_window(main_window_controller)
        main_window_controller.show()

        app.exec()
    except Exception as e:
        show_error(global_settings, "An error occurred during application initialization", e)

if __name__ == '__main__':
    main()
