import sys
import os
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QObject, QEvent, Qt
from models.GlobalSettings import GlobalSettings
from utils.ui import show_error
from controllers.MainWindowController import MainWindowController

def main():
    app = QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")

    if hasattr(Qt.ApplicationAttribute, 'AA_UseHighDpiPixmaps'):
        app.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps)

    app_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    global_settings = GlobalSettings(app_dir)

    # Apply the saved theme (or default to light if it's the first startup)
    global_settings.apply_theme()

    print(f"App directory: {global_settings.get_app_dir()}")
    print(f"Src directory: {global_settings.get_src_dir()}")
    print(f"UI directory: {global_settings.get_ui_dir()}")
    print(f"Controllers directory: {global_settings.get_controllers_dir()}")
    print(f"CSPR_DB path: {global_settings.CSPR_DB}")
    print(f"CASPERinfo path: {global_settings.get_casper_info_path()}")
    print(f"Assets directory: {global_settings.get_assets_dir()}")

    try:
        main_window_controller = MainWindowController(global_settings)
        global_settings.set_main_window(main_window_controller)  # Set the main window
        main_window_controller.show()

        sys.exit(app.exec())
    except Exception as e:
        show_error(global_settings, "An error occurred during application initialization", e)

if __name__ == '__main__':
    main()
