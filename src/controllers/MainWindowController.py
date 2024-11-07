from PyQt6 import QtWidgets, QtCore, QtGui
from PyQt6.QtWidgets import QWidget, QVBoxLayout
from views.MainWindowView import MainWindowView
from models.MainWindowModel import MainWindowModel
from utils.ui import show_error, show_message 
from utils.web import ncbi_page, repo_page, ncbi_blast_page 
from PyQt6.QtCore import Qt
from PyQt6.QtCore import QSize
from utils.LoggingMixin import LoggingMixin

class MainWindowController(LoggingMixin):
    def __init__(self, global_settings):
        LoggingMixin.__init__(self)
        self.global_settings = global_settings
        self.tab_widgets = {
            'widgets': {},  
            'controllers': {}  
        }
        self.startup_controller = None
        self.is_first_time_startup = self.global_settings.is_first_time_startup
        self.shared_tab_size = QSize(850, 850)
        self.startup_size = QSize(750, 550)
        self.current_tab = None

        try:
            self.view = MainWindowView(global_settings)
            self._setup_connections()
            self._init_ui()
            self.global_settings.check_and_emit_first_time_startup()
        except Exception as e:
            self.log_error("__init__", e)
            show_error(self.global_settings, "Error initializing MainWindowController", str(e))

    def _setup_connections(self):
        self.log_method_call("_setup_connections")
        
        # menuBar
        self.view.action_change_database_directory.triggered.connect(self._change_database_directory)
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

    def _init_ui(self):
        self.log_method_call("_init_ui")
        
        if self.is_first_time_startup:
            self.log_info("First time startup detected. Opening startup tab.")
            self._open_startup_tab()
            return

        db_path = self.global_settings.get_db_path()
        is_valid, message = self.global_settings.validate_db_path(db_path)
        
        if db_path and is_valid:
            self.log_info(f"Database path is valid: {db_path}")
            self._open_home_tab()
        else:
            self.log_warning(f"Invalid database path: {db_path}. {message}")
            self._open_startup_tab()

    def _handle_first_time_startup(self):
        self.log_info("First time startup signal received")
        self.is_first_time_startup = True
        self._open_startup_tab()

    def _open_startup_tab(self):
        try:
            self.startup_controller = self.global_settings.get_startup_window()
            self.open_new_tab("Startup", self.startup_controller)
        except Exception as e:
            self.log_error("_open_startup_tab", e)
            show_error(self.global_settings, "Error opening startup tab", str(e))

    def _switch_to_home_from_startup(self):
        self.log_method_call("_switch_to_home_from_startup")
        
        startup_tab = self.find_tab_by_title("Startup")
        if startup_tab:
            index = self.view.tab_widget.indexOf(startup_tab)
            self._close_tab(index)
            
            if self.startup_controller:
                self.startup_controller.deactivate()
                self.startup_controller = None
        else:
            self.log_warning("Startup tab not found when trying to close it")

        self.close_new_genome_and_switch_to_home()
        self._center_window()

    def _center_window(self):
        try:
            center_point = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
            frame_geometry = self.view.frameGeometry()
            frame_geometry.moveCenter(center_point)
            self.view.move(frame_geometry.topLeft())
            self.log_debug(f"Window centered at {self.view.pos()}")
        except Exception as e:
            self.log_error("_center_window", e)
            show_error(self.global_settings, "Error centering window", str(e))

    def _change_database_directory(self):
        try:
            new_directory = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Select Database Directory", 
                self.global_settings.get_db_path(),
                QtWidgets.QFileDialog.Option.ShowDirsOnly
            )
            
            if not new_directory:
                return

            is_valid, message = self.global_settings.validate_db_path(new_directory)
            if is_valid:
                self._process_valid_directory(new_directory)
            else:
                self._handle_invalid_directory(new_directory, message)
                
        except Exception as e:
            self.log_error("_change_database_directory", e)
            show_error(self.global_settings, "Error changing database directory", str(e))

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

    def _process_valid_directory(self, new_directory):
        try:
            self.global_settings.save_db_path(new_directory)
            self.global_settings.update_db_state()
            show_message("Success", "Database directory changed successfully.")
            
            if (self.startup_controller and 
                self.view.tab_widget.currentWidget() == self.startup_controller.view):
                self._switch_to_home_from_startup()
        except Exception as e:
            self.log_error("_process_valid_directory", e)
            show_error(self.global_settings, "Error processing directory", str(e))

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
        Handle the tab_closed signal from CloseableTabWidget
        """
        try:
            # Remove references from both widgets and controllers dictionaries
            for title in list(self.tab_widgets['widgets'].keys()):
                if self.tab_widgets['widgets'][title] == widget:
                    self.logger.info(f"Tab '{title}' closed. Cleaning up references.")
                    del self.tab_widgets['widgets'][title]
                    if title in self.tab_widgets['controllers']:
                        del self.tab_widgets['controllers'][title]
                    break
        except Exception as e:
            self.logger.error(f"Error in _on_tab_closed: {str(e)}")

    def _open_home_tab(self):
        """Opens the home tab"""
        try:
            home_controller = self.global_settings.get_home_window()
            self.open_new_tab("Home", home_controller)
            self.log_info("Home tab opened successfully")
        except Exception as e:
            self.log_error("_open_home_tab", e)
            show_error(self.global_settings, "Error opening home tab", str(e))

    def open_new_tab(self, title, content):
        """Opens a new tab with the given title and content"""
        try:
            self.log_debug(f"Opening new tab: {title}")
            
            # Check if the tab already exists
            existing_tab = self.find_tab_by_title(title)
            if existing_tab:
                self.log_debug(f"Tab '{title}' already exists, switching to it")
                self.view.tab_widget.setCurrentWidget(existing_tab)
                self._resize_for_tab(title)
                return

            # Create widget from content
            if hasattr(content, 'view'):
                widget = content.view
                # Store controller reference
                self.tab_widgets['controllers'][title] = content
            else:
                widget = content

            # Create wrapper widget with padding
            wrapper = QWidget()
            layout = QVBoxLayout(wrapper)
            layout.setContentsMargins(10, 10, 10, 10)
            layout.addWidget(widget)

            # Add the wrapper to the tab widget and store reference
            index = self.view.tab_widget.addTab(wrapper, title)
            self.view.tab_widget.setCurrentIndex(index)
            self.tab_widgets['widgets'][title] = wrapper

            self._resize_for_tab(title)
            self.log_info(f"Tab '{title}' opened successfully at index {index}")
            
        except Exception as e:
            self.log_error("open_new_tab", e)
            show_error(self.global_settings, f"Error opening tab '{title}'", str(e))

    def _resize_for_tab(self, title):
        if title == "Startup":
            # For Startup tab, set fixed size and disable maximize button
            self.view.setFixedSize(self.startup_size)
            self.view.setWindowFlags(self.view.windowFlags() & ~Qt.WindowType.WindowMaximizeButtonHint)
        else:
            # For all other tabs, use the shared size and allow resizing
            self.view.setMinimumSize(QSize(400, 300))
            self.view.setMaximumSize(QtCore.QSize(16777215, 16777215))
            self.view.setWindowFlags(self.view.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint)
            
            # Only resize if coming from Startup tab or if no current size is set
            if self.current_tab == "Startup" or self.view.size() == self.startup_size:
                self.view.resize(self.shared_tab_size)
        
        # Ensure window flags are updated
        self.view.show()
        
        # Update the current tab
        self.current_tab = title

    def _close_tab(self, index):
        """
        Handle tab closure using CloseableTabWidget
        """
        if 0 <= index < self.view.tab_widget.count():
            title = self.view.tab_widget.tabText(index)
            
            # Let CloseableTabWidget handle the widget cleanup
            self.view.tab_widget.closeTab(index)
            
            # Clean up our references
            if title in self.tab_widgets['widgets']:
                del self.tab_widgets['widgets'][title]
            if title in self.tab_widgets['controllers']:
                del self.tab_widgets['controllers'][title]

            self.logger.debug(f"Closed tab '{title}' at index {index}")

            # Handle post-close operations
            if title == "New Genome":
                home_tab = self.find_tab_by_title("Home")
                if home_tab:
                    home_controller = self.global_settings.get_home_window()
                    home_controller.refresh_data()

            # Resize for the current tab
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
                # center_ui(self.view)
                pass
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
        """Find a tab by its title"""
        return self.tab_widgets['widgets'].get(title)

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




