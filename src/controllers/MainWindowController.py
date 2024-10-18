import os
from PyQt6 import QtWidgets, QtCore, uic
from PyQt6.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout
from views.MainWindowView import MainWindowView
from models.MainWindowModel import MainWindowModel
from controllers.MultitargetingWindowController import MultitargetingWindowController
from utils.ui import show_error, show_message, scale_ui, center_ui, position_window
from utils.web import ncbi_page, repo_page, ncbi_blast_page 
from PyQt6.QtCore import QObject, Qt
import qdarktheme

class MainWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.tab_widgets = {}  # Store references to tab widgets
        self.startup_controller = None
        self.is_first_time_startup = self.global_settings.is_first_time_startup

        try:
            self.view = MainWindowView(global_settings)
            self._setup_connections()
            self._init_ui()
            
            # Check and emit first_time_startup signal after initialization
            self.global_settings.check_and_emit_first_time_startup()
        except Exception as e:
            show_error(self.global_settings, "Error initializing MainWindowController", str(e))

    def _setup_connections(self):
        try:
            # menuBar
            self.view.action_change_database_directory.triggered.connect(self._change_database_directory)
            # self.view.action_open_genome_browser.triggered.connect(self.open_genome_browser)
            self.view.action_open_repository.triggered.connect(self._open_repository_website)
            self.view.action_open_NCBI_BLAST.triggered.connect(self._open_ncbi_blast_website)
            self.view.action_open_NCBI.triggered.connect(self._open_ncbi_website)

            # Title Bar
            self.view.close_window_button.clicked.connect(self._close_window)
            self.view.minimize_window_button.clicked.connect(self._minimize_window)
            self.view.maximize_window_button.clicked.connect(self._maximize_window)
            self.view.theme_toggle_button.clicked.connect(self._toggle_theme)

            # Tab bar
            self.view.tab_widget.tab_closed.connect(self._on_tab_closed)
            self.view.tab_widget.tabCloseRequested.connect(self._close_tab)

            self.global_settings.first_time_startup.connect(self._handle_first_time_startup)

        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in MainWindowController", str(e))

    def _init_ui(self):
        try:
            if self.is_first_time_startup:
                self.logger.info("First time startup detected in _init_ui. Opening startup tab.")
                self._open_startup_tab()
            else:
                db_path = self.global_settings.get_db_path()
                is_valid, message = self.global_settings.validate_db_path(db_path)
                if db_path and is_valid:
                    self.logger.info(f"Database path is set and valid: {db_path}")
                    self._open_home_tab()
                else:
                    self.logger.info(f"Database path is not set or invalid: {db_path}. {message}")
                    self._open_startup_tab()
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in MainWindowController", str(e))

    def _handle_first_time_startup(self):
        self.logger.info("First time startup signal received. Opening startup tab.")
        self.is_first_time_startup = True
        self._open_startup_tab()

    def _open_startup_tab(self):
        try:
            self.startup_controller = self.global_settings.get_startup_window()
            self.open_new_tab("Startup", self.startup_controller)
            self.view.setFixedSize(750, 550)
            self.view.setWindowFlags(self.view.windowFlags() & ~Qt.WindowType.WindowMaximizeButtonHint)
        except Exception as e:
            show_error(self.global_settings, "Error opening startup tab", str(e))

    def _switch_to_home_from_startup(self):
        try:
            self.logger.debug("Switching to home from startup")
            # Find the startup tab
            startup_tab = None
            for i in range(self.view.tab_widget.count()):
                if self.view.tab_widget.tabText(i) == "Startup":
                    startup_tab = self.view.tab_widget.widget(i)
                    break

            if startup_tab:
                index = self.view.tab_widget.indexOf(startup_tab)
                self._close_tab(index)
                self.logger.debug(f"Closed startup tab at index {index}")
            else:
                self.logger.warning("Startup tab not found when trying to close it")

            self._open_home_tab()
        except Exception as e:
            self.logger.error(f"Error switching to home from startup: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error switching to home tab", str(e))

    def _open_home_tab(self):
        try:
            home_controller = self.global_settings.get_home_window()
            self.open_new_tab("Home", home_controller)
        except Exception as e:
            show_error(self.global_settings, "Error opening home tab", str(e))

    def _change_database_directory(self):
        try:
            new_directory = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Select Database Directory", self.global_settings.CSPR_DB,
                QtWidgets.QFileDialog.Option.ShowDirsOnly
            )
            if new_directory:
                if self.global_settings.validate_db_path(new_directory):
                    self.global_settings.save_database_path(new_directory)
                    self.global_settings.initialize_app_directories()
                    self.load_combo_box_data()
                    show_message(12, QtWidgets.QMessageBox.Icon.Information,
                                 "Success", "Database directory changed successfully.")
                else:
                    show_message(12, QtWidgets.QMessageBox.Icon.Warning,
                                 "Invalid Directory", "The selected directory does not contain valid CSPR files.")
        except Exception as e:
            self.logger.error(f"Error changing database directory: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error changing database directory", str(e))
    
    def _open_ncbi_website(self):
        ncbi_page()

    def _open_repository_website(self):
        repo_page()

    def _open_ncbi_blast_website(self):
        ncbi_blast_page()

    def _close_window(self):
        self.view.close()

    def _minimize_window(self):
        self.view.showMinimized()

    def _maximize_window(self):
        if self.view.isMaximized():
            self.view.showNormal()
        else:
            self.view.showMaximized()

    def _on_tab_closed(self, widget):
        """
        Handle the tab_closed signal by removing references to the deleted widget.
        """
        # Iterate through the tab_widgets to find and remove the closed widget
        for title, tab_widget in list(self.tab_widgets.items()):
            if tab_widget == widget:
                self.logger.info(f"Tab '{title}' closed. Dereferencing the widget.")
                del self.tab_widgets[title]
                break

    def open_new_tab(self, title, content):
        try:
            self.logger.debug(f"Attempting to open new tab: {title}")
            
            if title in self.tab_widgets:
                self.logger.debug(f"Tab '{title}' already exists, switching to it")
                existing_widget = self.tab_widgets[title]
                self.view.tab_widget.setCurrentWidget(existing_widget)
                return

            # Determine if content is a controller or a widget
            if hasattr(content, 'view'):
                self.logger.debug("Content is a controller, getting its view")
                widget = content.view
            else:
                self.logger.debug("Content is a widget")
                widget = content

            # Create a wrapper widget with padding
            self.logger.debug("Creating wrapper widget with padding")
            wrapper = QWidget()
            layout = QVBoxLayout(wrapper)
            layout.setContentsMargins(10, 10, 10, 10)
            layout.addWidget(widget)

            # Add the wrapper to the tab widget
            self.logger.debug("Adding wrapper to tab widget")
            index = self.view.tab_widget.addTab(wrapper, title)
            self.view.tab_widget.setCurrentIndex(index)
            self.tab_widgets[title] = wrapper

            # Resize the main window if it's the New Genome tab
            if title == "New Genome":
                self.logger.debug("Resizing main window for New Genome tab")
                new_size = widget.sizeHint()
                self.view.resize(new_size)

            # center_ui(self.view)  # Re-center the window after resizing
            self.logger.info(f"Opened new tab '{title}' at index {index}")
        except Exception as e:
            self.logger.error(f"Error opening tab '{title}': {str(e)}", exc_info=True)
            show_error(self.global_settings, f"Error opening tab '{title}'", str(e))

    def _close_tab(self, index):
        if 0 <= index < self.view.tab_widget.count():
            widget = self.view.tab_widget.widget(index)
            title = self.view.tab_widget.tabText(index)
            self.view.tab_widget.removeTab(index)
            if widget:
                widget.deleteLater()
            # Remove the tab from our tab_widgets dictionary
            if title in self.tab_widgets:
                del self.tab_widgets[title]
            self.logger.debug(f"Closed tab '{title}' at index {index}")
        else:
            self.logger.warning(f"Attempted to close non-existent tab at index {index}")

    def _toggle_theme(self):
        try:
            self.global_settings.set_theme("dark" if self.global_settings.get_theme() == "light" else "light")
            self.view.update_theme_icon()
            self.view.apply_theme()
        except Exception as e:
            show_error(self.global_settings, "Error toggling theme", str(e))

    def show(self):
        try:
            saved_position = self.global_settings.load_window_position("main_window")
            if saved_position:
                self.view.move(saved_position)
            else:
                center_ui(self.view)
            self.view.show()
            self.view.apply_theme()
        except Exception as e:
            self.global_settings.logger.error(f"Error showing main window: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error showing main window", e)

    def open_new_genome_tab(self):
        tab_title = "New Genome"
        existing_tab = self.find_tab_by_title(tab_title)
        if existing_tab:
            self.view.tab_widget.setCurrentWidget(existing_tab)
        else:
            new_genome_controller = self.global_settings.get_new_genome_window()
            self.open_new_tab(tab_title, new_genome_controller)

    def find_tab_by_title(self, title):
        for i in range(self.view.tab_widget.count()):
            if self.view.tab_widget.tabText(i) == title:
                return self.view.tab_widget.widget(i)
        return None
