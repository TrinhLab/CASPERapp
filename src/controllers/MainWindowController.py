import os
from PyQt6 import QtWidgets, QtCore, uic, QtGui
from PyQt6.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout
from views.MainWindowView import MainWindowView
from models.MainWindowModel import MainWindowModel
from controllers.MultitargetingWindowController import MultitargetingWindowController
from utils.ui import show_error, show_message, scale_ui, center_ui, position_window
from utils.web import ncbi_page, repo_page, ncbi_blast_page 
from PyQt6.QtCore import QObject, Qt
import qdarktheme
from PyQt6.QtCore import QSize

class MainWindowController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.tab_widgets = {}  # Store references to tab widgets
        self.startup_controller = None
        self.is_first_time_startup = self.global_settings.is_first_time_startup
        self.tab_sizes = {
            "Startup": QSize(750, 550),
            "New Genome": QSize(575, 700),
            "Home": QSize(1000, 700),  # Add a default size for Home tab
            "Define New Endonuclease": QSize(450, 600),
            "NCBI Download Tool": QSize(1000, 700),
            "View Targets": QSize(1500, 700),
        }
        self.current_tab = None

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
        except Exception as e:
            show_error(self.global_settings, "Error opening startup tab", str(e))

    def _switch_to_home_from_startup(self):
        try:
            self.logger.debug("Switching to home from startup")
            # Find the startup tab
            startup_tab = self.find_tab_by_title("Startup")
            if startup_tab:
                index = self.view.tab_widget.indexOf(startup_tab)
                self._close_tab(index)
                self.logger.debug(f"Closed startup tab at index {index}")
                
                # Deactivate the startup controller
                if self.startup_controller:
                    self.startup_controller.deactivate()
                    self.startup_controller = None
            else:
                self.logger.warning("Startup tab not found when trying to close it")

            self.close_new_genome_and_switch_to_home()
            self._center_window()  # Center the window after switching to home
        except Exception as e:
            self.logger.error(f"Error switching to home from startup: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error switching to home tab", str(e))

    def _center_window(self):
        try:
            center_point = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
            frame_geometry = self.view.frameGeometry()
            frame_geometry.moveCenter(center_point)
            self.view.move(frame_geometry.topLeft())
            self.logger.debug(f"Centered window. New position: {self.view.pos()}")
        except Exception as e:
            self.logger.error(f"Error centering window: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error centering window", str(e))

    def _open_home_tab(self):
        try:
            home_controller = self.global_settings.get_home_window()
            self.open_new_tab("Home", home_controller)
        except Exception as e:
            show_error(self.global_settings, "Error opening home tab", str(e))

    def _change_database_directory(self):
        try:
            new_directory = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Select Database Directory", self.global_settings.get_db_path(),
                QtWidgets.QFileDialog.Option.ShowDirsOnly
            )
            if new_directory:
                is_valid, message = self.global_settings.validate_db_path(new_directory)
                if is_valid:
                    self._process_valid_directory(new_directory)
                else:
                    self._handle_invalid_directory(new_directory, message)
        except Exception as e:
            self.logger.error(f"Error changing database directory: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error changing database directory", str(e))

    def _process_valid_directory(self, new_directory):
        self.global_settings.save_db_path(new_directory)
        self.global_settings.update_db_state()
        show_message("Success", "Database directory changed successfully.")
        
        # If we're currently on the startup tab, switch to the home tab
        if self.startup_controller and self.view.tab_widget.currentWidget() == self.startup_controller.view:
            self._switch_to_home_from_startup()

    def _handle_invalid_directory(self, new_directory, message):
        reply = QtWidgets.QMessageBox.question(
            self.view,
            "Invalid Directory",
            f"The selected directory does not contain valid CSPR files: {message}\n\n"
            "Would you like to analyze a new genome in this directory?",
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No
        )
        
        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            self.global_settings.save_db_path(new_directory)
            self.global_settings.update_db_state()
            self.open_new_genome_tab()
        else:
            show_message("Operation Cancelled", "Database directory change cancelled.")

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
            
            # Check if the tab already exists
            existing_tab = self.find_tab_by_title(title)
            if existing_tab:
                self.logger.debug(f"Tab '{title}' already exists, switching to it")
                self.view.tab_widget.setCurrentWidget(existing_tab)
                self._resize_for_tab(title)
                return

            # If the tab doesn't exist, create a new one
            if hasattr(content, 'view'):
                widget = content.view
            else:
                widget = content

            # Create a wrapper widget with padding
            wrapper = QWidget()
            layout = QVBoxLayout(wrapper)
            layout.setContentsMargins(10, 10, 10, 10)
            layout.addWidget(widget)

            # Add the wrapper to the tab widget
            index = self.view.tab_widget.addTab(wrapper, title)
            self.view.tab_widget.setCurrentIndex(index)
            self.tab_widgets[title] = wrapper

            self._resize_for_tab(title)

            self.logger.info(f"Opened new tab '{title}' at index {index}")
        except Exception as e:
            self.logger.error(f"Error opening tab '{title}': {str(e)}", exc_info=True)
            show_error(self.global_settings, f"Error opening tab '{title}'", str(e))

        self.view.tab_widget.currentChanged.connect(self._on_tab_changed)

    def _resize_for_tab(self, title):
        if title in self.tab_sizes:
            new_size = self.tab_sizes[title]
            if title == "Startup":
                # For Startup tab, set fixed size and disable maximize button
                self.view.setFixedSize(new_size)
                self.view.setWindowFlags(self.view.windowFlags() & ~Qt.WindowType.WindowMaximizeButtonHint)
            else:
                # For other tabs, allow resizing but set a minimum size
                self.view.setMinimumSize(QSize(400, 300))  # Set a reasonable minimum size
                self.view.setMaximumSize(QtCore.QSize(16777215, 16777215))
                self.view.setWindowFlags(self.view.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint)
                
                # Resize to the specified size for the tab
                self.view.resize(new_size)
        else:
            # Default behavior for unknown tabs
            self.view.setMinimumSize(QSize(400, 300))  # Set a reasonable minimum size
            self.view.setMaximumSize(QtCore.QSize(16777215, 16777215))
            self.view.setWindowFlags(self.view.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint)
        
        # Ensure window flags are updated
        self.view.show()
        
        # Update the current tab
        self.current_tab = title

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

            # If we're closing the New Genome tab and Home tab exists, refresh it
            if title == "New Genome":
                home_tab = self.find_tab_by_title("Home")
                if home_tab:
                    home_controller = self.global_settings.get_home_window()
                    home_controller.refresh_data()
        else:
            self.logger.warning(f"Attempted to close non-existent tab at index {index}")

        if self.view.tab_widget.count() > 0:
            new_index = self.view.tab_widget.currentIndex()
            new_tab_title = self.view.tab_widget.tabText(new_index)
            self._resize_for_tab(new_tab_title)

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
        # Check if the New Genome tab already exists
        existing_tab = self.find_tab_by_title("New Genome")
        if existing_tab:
            # If it exists, just switch to it
            self.view.tab_widget.setCurrentWidget(existing_tab)
        else:
            # If it doesn't exist, create a new one
            new_genome_controller = self.global_settings.get_new_genome_window()
            new_genome_view = new_genome_controller.view
            tab_index = self.view.tab_widget.addTab(new_genome_view, "New Genome")
            self.view.tab_widget.setCurrentIndex(tab_index)
            self.tab_widgets["New Genome"] = new_genome_view
        
        self._resize_for_tab("New Genome")
        
        # Ensure the window is visible and brought to front
        self.view.show()
        self.view.raise_()
        self.view.activateWindow()

        # Log the current state
        self.logger.debug(f"Window visibility after opening New Genome tab: {self.view.isVisible()}")
        self.logger.debug(f"Window geometry after opening New Genome tab: {self.view.geometry()}")

    def find_tab_by_title(self, title):
        for i in range(self.view.tab_widget.count()):
            if self.view.tab_widget.tabText(i) == title:
                return self.view.tab_widget.widget(i)
        return None

    def _on_tab_changed(self, index):
        # Save the current tab size before switching
        if self.current_tab:
            current_size = self.view.size()
            if current_size.width() >= 400 and current_size.height() >= 300:
                self.tab_sizes[self.current_tab] = current_size

        # Get the new tab title
        new_tab_title = self.view.tab_widget.tabText(index)

        # Resize for the new tab
        self._resize_for_tab(new_tab_title)

    def close_new_genome_and_switch_to_home(self):
        try:
            self.logger.debug("Attempting to close New Genome tab and switch to Home")
            
            # Find and close the New Genome tab if it's open
            new_genome_tab = self.find_tab_by_title("New Genome")
            if new_genome_tab:
                index = self.view.tab_widget.indexOf(new_genome_tab)
                self._close_tab(index)
                self.logger.debug("Closed New Genome tab")
            else:
                self.logger.debug("New Genome tab not found")

            # Switch to the Home tab or create it if it doesn't exist
            home_tab = self.find_tab_by_title("Home")
            if home_tab:
                self.view.tab_widget.setCurrentWidget(home_tab)
                self.logger.debug("Switched to existing Home tab")
            else:
                self._open_home_tab()
                self.logger.debug("Opened new Home tab")

            # Resize for the Home tab
            self._resize_for_tab("Home")

        except Exception as e:
            self.logger.error(f"Error in close_new_genome_and_switch_to_home: {str(e)}", exc_info=True)
            show_error(self.global_settings, "Error switching to Home tab", str(e))


