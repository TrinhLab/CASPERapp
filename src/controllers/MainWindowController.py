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
        self.settings = global_settings
        self.tab_widgets = {
            'widgets': {},  
            'controllers': {}  
        }
        self.startup_controller = None
        self.is_first_time_startup = self.settings.is_first_time_startup
        self.shared_tab_size = QSize(850, 850)
        self.startup_size = QSize(750, 550)
        self.current_tab = None
        self.previous_size = None
        # Add flag to track dialog state
        self._cspr_deletion_dialog_shown = False
        
        # Define minimum sizes for different tab types
        self.tab_min_sizes = {
            "Startup": QSize(750, 550),  # Startup has a fixed size
            "View Targets": QSize(1300, 800),  # Large tabs
            "Multitargeting Analysis": QSize(1300, 800),
            "Population Analysis": QSize(1000, 700),  # Medium-large tabs
            "Find Targets": QSize(1000, 700),
            "New Genome": QSize(800, 600),  # Medium tabs
            "Home": QSize(800, 600),
            "default": QSize(600, 500)  # Default minimum size for any other tab
        }

        try:
            self.view = MainWindowView(self.settings)
            self._setup_connections()
            self._init_ui()
            self.settings.check_and_emit_first_time_startup()
            # Center the window after initialization
            self._center_window()
        except Exception as e:
            self.log_error("__init__", e)
            show_error(self.settings, "Error initializing MainWindowController", str(e))

    def _setup_connections(self):
        self.log_method_call("_setup_connections")
        
        # menuBar
        self.view.action_change_database_directory.triggered.connect(self._change_database_directory)
        self.view.action_open_repository.triggered.connect(self._open_repository_website)
        self.view.action_open_NCBI_BLAST.triggered.connect(self._open_ncbi_blast_website)
        self.view.action_open_NCBI.triggered.connect(self._open_ncbi_website)

        # Tab bar
        self.view.tab_widget.tab_closed.connect(self._on_tab_closed)
        self.view.tab_widget.tab_closing.connect(self._close_tab)
        self.view.tab_widget.currentChanged.connect(self._on_current_tab_changed)

        self.settings.first_time_startup.connect(self._handle_first_time_startup)

        # Add Button Menu
        self.view.action_new_endonuclease.triggered.connect(self.open_new_endonuclease_window)

        # Settings Menu
        self.view.action_toggle_theme.triggered.connect(self._toggle_theme)

        # Database state changes
        self.settings.db_manager.db_state_changed.connect(self._on_db_state_changed)
        self.settings.db_manager.db_validation_changed.connect(self._on_db_validation_changed)

    def _init_ui(self):
        self.log_method_call("_init_ui")
        
        # Check if it's first time startup
        if self.is_first_time_startup:
            self.log_info("First time startup detected. Opening startup tab.")
            self._open_startup_tab()
            return

        db_path = self.settings.get_db_path()
        is_valid, message = self.settings.validate_db_path(db_path)
        
        if not is_valid:
            self.log_warning(f"Invalid database path: {db_path}. {message}")
            # Always open startup tab for invalid paths, keeping the existing path
            self._open_startup_tab(keep_db_path=True)
            return  # Add explicit return to prevent further execution
        
        self.log_info(f"Database path is valid: {db_path}")
        self._open_home_tab()

    def _handle_first_time_startup(self):
        self.log_info("First time startup signal received")
        self.is_first_time_startup = True
        self._open_startup_tab()

    def _open_startup_tab(self, keep_db_path=False):
        try:
            # First safely deactivate any existing controllers
            for title, controller in list(self.tab_widgets['controllers'].items()):
                try:
                    if hasattr(controller, 'deactivate'):
                        controller.deactivate()
                except Exception as e:
                    self.logger.error(f"Error deactivating controller for {title}: {str(e)}")
            
            # Force close all tabs
            while self.view.tab_widget.count() > 0:
                try:
                    widget = self.view.tab_widget.widget(0)
                    self.view.tab_widget.removeTab(0)
                    if widget:
                        widget.deleteLater()
                except Exception as e:
                    self.logger.error(f"Error closing tab: {str(e)}")
            
            # Clear tab tracking dictionaries
            self.tab_widgets['widgets'].clear()
            self.tab_widgets['controllers'].clear()
            
            # Then open the startup tab
            self.startup_controller = self.settings.get_startup_window(keep_db_path=keep_db_path)
            self.open_new_tab("Startup", self.startup_controller)
        except Exception as e:
            self.log_error("_open_startup_tab", e)
            show_error(self.settings, "Error opening startup tab", str(e))

    def _switch_to_home_from_startup(self):
        self.log_method_call("_switch_to_home_from_startup")
        
        # First deactivate startup controller
        if self.startup_controller:
            self.startup_controller.deactivate()
            self.startup_controller = None

        # Close startup tab if it exists
        startup_tab = self.find_tab_by_title("Startup")
        if startup_tab:
            index = self.view.tab_widget.indexOf(startup_tab)
            self._close_tab(index)
        else:
            self.log_warning("Startup tab not found when trying to close it")

        # Open home tab and ensure it's properly initialized
        self._open_home_tab()
        
        # Center the window after all tab operations
        self._center_window()

    def _center_window(self):
        """Center the window on the current screen"""
        try:
            # Get the current screen where the window is or the primary screen
            window_screen = self.view.screen()
            if not window_screen:
                window_screen = QtGui.QGuiApplication.primaryScreen()
            
            # Get the geometry of the screen
            screen_geometry = window_screen.availableGeometry()
            
            # Calculate the center point
            center_point = screen_geometry.center()
            
            # Get the window geometry
            frame_geometry = self.view.frameGeometry()
            
            # Move the window's center to the screen's center
            frame_geometry.moveCenter(center_point)
            
            # Move the window to the calculated position
            self.view.move(frame_geometry.topLeft())
            
            self.log_debug(f"Window centered on screen at {self.view.pos()}")
        except Exception as e:
            self.log_error("_center_window", e)
            show_error(self.settings, "Error centering window", str(e))

    def _change_database_directory(self):
        try:
            new_directory = QtWidgets.QFileDialog.getExistingDirectory(
                self.view, "Select Database Directory", 
                self.settings.get_db_path(),
                QtWidgets.QFileDialog.Option.ShowDirsOnly
            )
            
            if not new_directory:
                return

            is_valid, message = self.settings.validate_db_path(new_directory)
            if is_valid:
                self._process_valid_directory(new_directory)
            else:
                self._handle_invalid_directory(new_directory, message)
                
        except Exception as e:
            self.log_error("_change_database_directory", e)
            show_error(self.settings, "Error changing database directory", str(e))

    def _handle_invalid_directory(self, new_directory, message):
        """Handle invalid directory selection with option to analyze new genomes"""
        self.logger.debug("Entering _handle_invalid_directory")  # Add entry log
        
        # Always show error for non-existent directories
        if message == "The selected directory does not exist.":
            self.logger.debug("Directory does not exist, showing error")  # Add debug log
            show_error(self.settings, "Invalid Directory", message)
            return
        
        self.logger.debug("Showing analyze new genome dialog")  # Add debug log
        # Show dialog for analyzing new genome
        reply = QtWidgets.QMessageBox.question(
            self.view,
            "Invalid Directory",
            "Would you like to analyze a new genome in this directory? testing",
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No
        )
        
        self.logger.debug(f"User reply to analyze new genome: {reply == QtWidgets.QMessageBox.StandardButton.Yes}")  # Add debug log
        
        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            self.logger.debug("User chose to analyze new genome")  # Add debug log
            # Temporarily disconnect validation signal to prevent warning
            self.settings.db_manager.db_validation_changed.disconnect()
            try:
                # Set directory change flag
                self.settings.db_manager.is_changing_directory = True
                self.logger.debug(f"Starting directory change process to: {new_directory}")
                
                # Save the new path without switching to startup tab
                self.settings.save_db_path(new_directory)
                self.settings.update_db_state()
                
                # Store the new directory for later use
                self._pending_directory_change = new_directory
                
                # Open new genome tab without closing other tabs
                self.logger.debug("Opening new genome tab")  # Add debug log
                self.open_new_genome_tab()
                
                # Connect to the new genome tab's completion signal
                new_genome_tab = self.find_tab_by_title("New Genome")
                if new_genome_tab and new_genome_tab in self.tab_widgets['controllers']:
                    new_genome_controller = self.tab_widgets['controllers']["New Genome"]
                    new_genome_controller.process_completed.connect(self._on_new_genome_completed)
                    self.logger.debug("Connected to new genome completion signal")  # Add debug log
                
            except Exception as e:
                self.logger.error(f"Error in _handle_invalid_directory: {str(e)}")  # Add error log
                raise
            finally:
                # Reconnect the signal after operation is complete
                self.settings.db_manager.db_validation_changed.connect(self._on_db_validation_changed)
                self.logger.debug("Reconnected validation signal")  # Add debug log
        else:
            self.logger.debug("User cancelled new genome analysis")  # Add debug log
            show_message("Operation Cancelled", "Database directory change cancelled.")

    def _on_new_genome_completed(self):
        """Handle completion of new genome analysis"""
        try:
            if hasattr(self, '_pending_directory_change'):
                # Show success message
                show_message(
                    "Database Directory Change Complete",
                    f"Successfully changed database directory to:\n{self._pending_directory_change}",
                    QtWidgets.QMessageBox.Icon.Information
                )
                
                # Refresh home tab
                home_tab = self.find_tab_by_title("Home")
                if home_tab and home_tab in self.tab_widgets['controllers']:
                    home_controller = self.tab_widgets['controllers']["Home"]
                    home_controller.refresh_data()
                
                # Clean up
                delattr(self, '_pending_directory_change')
        except Exception as e:
            self.logger.error(f"Error handling new genome completion: {str(e)}")

    def _process_valid_directory(self, new_directory):
        try:
            self.settings.save_db_path(new_directory)
            self.settings.update_db_state()
            show_message("Success", "Database directory changed successfully.")
            
            if (self.startup_controller and 
                self.view.tab_widget.currentWidget() == self.startup_controller.view):
                self._switch_to_home_from_startup()
        except Exception as e:
            self.log_error("_process_valid_directory", e)
            show_error(self.settings, "Error processing directory", str(e))

    def _open_ncbi_website(self):
        ncbi_page()

    def _open_repository_website(self):
        repo_page()
 
    def _open_ncbi_blast_website(self):
        ncbi_blast_page()

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
            home_controller = self.settings.get_home_window()
            self.open_new_tab("Home", home_controller)
            self.log_info("Home tab opened successfully")
        except Exception as e:
            self.log_error("_open_home_tab", e)
            show_error(self.settings, "Error opening home tab", str(e))

    def open_new_tab(self, title, content):
        """Opens a new tab with the given title and content"""
        try:
            self.log_debug(f"Opening new tab: {title}")
            
            # Check if the tab already exists
            existing_tab = self.find_tab_by_title(title)
            if existing_tab:
                self.log_debug(f"Tab '{title}' already exists, switching to it")
                self.view.tab_widget.setCurrentWidget(existing_tab)
                self._resize_for_tab(title, center_window=False)  # Don't center when switching to existing tab
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

            self._resize_for_tab(title, center_window=True)  # Center when opening new tab
            self.log_info(f"Tab '{title}' opened successfully at index {index}")
            
        except Exception as e:
            self.log_error("open_new_tab", e)
            show_error(self.settings, f"Error opening tab '{title}'", str(e))

    def _resize_for_tab(self, title, center_window=True):
        """Handle window resizing for different tab types"""
        try:
            # Get the minimum size for this tab type
            min_size = self.tab_min_sizes.get(title, self.tab_min_sizes["default"])
            
            if title == "Startup":
                # Startup tab has a fixed size
                self.view.setFixedSize(min_size)
            else:
                # Store current size before applying constraints if coming from a different size category
                current_min_size = self.tab_min_sizes.get(self.current_tab, self.tab_min_sizes["default"])
                if current_min_size != min_size:
                    self.previous_size = self.view.size()
                
                # Calculate new dimensions
                new_width = max(self.view.width(), min_size.width())
                new_height = max(self.view.height(), min_size.height())
                
                # Only resize if dimensions need to increase
                if new_width > self.view.width() or new_height > self.view.height():
                    self.view.resize(QSize(new_width, new_height))
                
                # Set size constraints
                self.view.setMinimumSize(min_size)
                self.view.setMaximumSize(QtCore.QSize(16777215, 16777215))
                
                # Restore previous size if available and coming from a larger minimum size tab
                if self.current_tab and self.previous_size:
                    current_min_size = self.tab_min_sizes.get(self.current_tab, self.tab_min_sizes["default"])
                    if current_min_size.width() > min_size.width() or current_min_size.height() > min_size.height():
                        # Only restore if the previous size is larger than our minimum
                        if (self.previous_size.width() >= min_size.width() and 
                            self.previous_size.height() >= min_size.height()):
                            self.view.resize(self.previous_size)
                elif self.current_tab == "Startup" or self.view.size() == self.startup_size:
                    self.view.resize(self.shared_tab_size)
            
            # Update the current tab
            self.current_tab = title
            
            # Center the window after resizing only if requested
            if center_window:
                self._center_window()
            
        except Exception as e:
            self.log_error("_resize_for_tab", e)

    def _close_tab(self, index):
        """Handle tab closure using CloseableTabWidget"""
        try:
            if 0 <= index < self.view.tab_widget.count():
                title = self.view.tab_widget.tabText(index)
                
                # Safely deactivate controller if it exists
                if title in self.tab_widgets['controllers']:
                    try:
                        controller = self.tab_widgets['controllers'][title]
                        if hasattr(controller, 'deactivate'):
                            controller.deactivate()
                    except Exception as e:
                        self.logger.error(f"Error deactivating controller for {title}: {str(e)}")
                
                # Clean up references
                if title in self.tab_widgets['widgets']:
                    del self.tab_widgets['widgets'][title]
                if title in self.tab_widgets['controllers']:
                    del self.tab_widgets['controllers'][title]

                self.logger.debug(f"Closed tab '{title}' at index {index}")
        except Exception as e:
            self.logger.error(f"Error in _close_tab: {str(e)}")

    def _toggle_theme(self):
        try:
            self.settings.set_theme("dark" if self.settings.get_theme() == "light" else "light")
            self.view.update_theme_icon()
            self.view.apply_theme()
        except Exception as e:
            show_error(self.settings, "Error toggling theme", str(e))

    def show(self):
        try:
            saved_position = self.settings.load_window_position("main_window")
            if saved_position:
                self.view.move(saved_position)
            self.view.show()
            self.view.apply_theme()
            # Center the window after showing it
            self._center_window()
            print("Window initialized")
        except Exception as e:
            self.log_error("show", e)
            show_error(self.settings, "Error showing main window", e)

    def open_new_genome_tab(self):
        # Check if the New Genome tab already exists
        existing_tab = self.find_tab_by_title("New Genome")
        if existing_tab:
            # If it exists, just switch to it
            self.view.tab_widget.setCurrentWidget(existing_tab)
        else:
            # If it doesn't exist, create a new one
            new_genome_controller = self.settings.get_new_genome_window()
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
            self.log_error("close_new_genome_and_switch_to_home", e)
            show_error(self.settings, "Error switching to Home tab", str(e))

    def _on_current_tab_changed(self, index):
        """Handle tab change events"""
        try:
            if index >= 0:
                new_tab_title = self.view.tab_widget.tabText(index)
                old_tab_title = self.current_tab
                
                # Store current size if coming from a non-Startup tab
                if old_tab_title and old_tab_title != "Startup":
                    self.previous_size = self.view.size()
                
                self._resize_for_tab(new_tab_title, center_window=False)
                
        except Exception as e:
            self.log_error("_on_current_tab_changed", e)

    def open_new_endonuclease_window(self):
        """Opens the new endonuclease as a separate window"""
        try:
            # Create the controller
            new_endonuclease_controller = self.settings.get_new_endonuclease_window()
            
            # Get the window from the controller
            window = new_endonuclease_controller.view
            
            # Set window properties
            window.setWindowModality(Qt.WindowModality.ApplicationModal)  # Make it modal
            window.setMinimumSize(QSize(800, 600))  # Set minimum size
            
            # Center the window relative to the main window
            main_window_center = self.view.geometry().center()
            window_geometry = window.frameGeometry()
            window_geometry.moveCenter(main_window_center)
            window.move(window_geometry.topLeft())
            
            # Show the window
            window.show()
            window.raise_()
            window.activateWindow()
            
            # Store reference to prevent garbage collection
            self._current_new_endonuclease_window = new_endonuclease_controller
            
            self.log_info("New Endonuclease window opened successfully")
        except Exception as e:
            self.log_error("open_new_endonuclease_window", e)
            show_error(self.settings, "Error opening new endonuclease window", str(e))

    def _on_db_validation_changed(self, is_valid, message):
        """Handle database validation state changes"""
        if not is_valid:
            self.logger.debug(f"Database validation failed: {message}")
            
            # Only show dialog if we're not in startup tab
            startup_tab = self.find_tab_by_title("Startup")
            if startup_tab and self.view.tab_widget.currentWidget() == startup_tab:
                return
            
            # Switch to startup tab for invalid paths (including when files are deleted)
            # but keep the current path
            if not startup_tab:
                self._open_startup_tab(keep_db_path=True)
                return
            
        else:
            self._cspr_deletion_dialog_shown = False  # Reset flag
            if "Successfully changed database directory to:" in message:
                show_message(
                    "Database Directory Change Complete",
                    message,
                    QtWidgets.QMessageBox.Icon.Information
                )
                # Refresh home tab if it exists
                home_tab = self.find_tab_by_title("Home")
                if home_tab and home_tab in self.tab_widgets['controllers']:
                    home_controller = self.tab_widgets['controllers']["Home"]
                    home_controller.refresh_data()

    def _on_db_state_changed(self, is_valid, message, changes):
        """Handle database state changes"""
        try:
            self.logger.debug(f"Database state changed - Valid: {is_valid}, Message: {message}, Changes: {changes}")
            
            # If path becomes invalid (including when files are deleted), switch to startup tab
            if not is_valid:
                # Skip only if we're already in startup tab
                startup_tab = self.find_tab_by_title("Startup")
                if startup_tab and self.view.tab_widget.currentWidget() == startup_tab:
                    return
                
                # Skip only if we're in the process of changing directory for new genome
                if self.settings.db_manager.is_changing_directory:
                    return
                
                # Otherwise, switch to startup tab with current path
                self._open_startup_tab(keep_db_path=True)
                return
            
            # Handle other state changes
            if changes:
                # Refresh home tab if it exists
                home_tab = self.find_tab_by_title("Home")
                if home_tab and home_tab in self.tab_widgets['controllers']:
                    home_controller = self.tab_widgets['controllers']["Home"]
                    home_controller.refresh_data()
                
        except Exception as e:
            self.logger.error(f"Error handling database state change: {str(e)}")

