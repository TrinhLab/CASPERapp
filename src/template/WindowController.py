from PyQt6.QtWidgets import QMainWindow
from .WindowView import WindowView
from .WindowModel import WindowModel
from utils.ui import show_error, show_message, position_window

class WindowController(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        try:
            self.model = WindowModel(global_settings)
            self.view = WindowView(global_settings)
            self.setCentralWidget(self.view)
            self._setup_connections()
            self._initialize_ui()
        except Exception as e:
            show_error(global_settings, "Error initializing WindowController", str(e))
            raise

    def _setup_connections(self):
        """Set up signal-slot connections between view elements and controller methods."""
        try:
            if self.view.push_button_action:
                self.view.push_button_action.clicked.connect(self.handle_action)
            if self.view.combo_box_selection:
                self.view.combo_box_selection.currentIndexChanged.connect(self.handle_selection_change)
            if self.view.custom_button:
                self.view.custom_button.clicked.connect(self.handle_custom_action)
            # Add more connections as needed
        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in WindowController", str(e))

    def _initialize_ui(self):
        """Initialize the UI with data from the model."""
        try:
            self.model.load_data()
            items = self.model.get_items()
            self.view.update_combo_box(items)
            # Add more initialization steps as needed
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in WindowController", str(e))

    def handle_action(self):
        """Handle the main action button click."""
        try:
            input_text = self.view.get_input_text()
            result = self.model.process_input(input_text)
            self.view.show_message(f"Processed result: {result}")
        except Exception as e:
            show_error(self.global_settings, "Error handling action", str(e))

    def handle_selection_change(self, index):
        """Handle changes in the combo box selection."""
        try:
            selected_item = self.view.combo_box_selection.currentText()
            # Process the selection change
            self.view.show_message(f"Selected: {selected_item}")
        except Exception as e:
            show_error(self.global_settings, "Error handling selection change", str(e))

    def handle_custom_action(self):
        """Handle the custom action button click."""
        try:
            # Implement custom action logic
            self.view.show_message("Custom action performed")
        except Exception as e:
            show_error(self.global_settings, "Error handling custom action", str(e))

    def save_data(self):
        """Save data using the model."""
        try:
            data_to_save = {"key": "value"}  # Collect data to save
            if self.model.save_data(data_to_save):
                self.view.show_message("Data saved successfully")
            else:
                self.view.show_message("Failed to save data")
        except Exception as e:
            show_error(self.global_settings, "Error saving data", str(e))

    def load_data(self):
        """Load data using the model and update the view."""
        try:
            self.model.load_data()
            # Update view with loaded data
            self.view.update_combo_box(self.model.get_items())
        except Exception as e:
            show_error(self.global_settings, "Error loading data", str(e))

    def show(self):
        """Show the window and position it on the screen."""
        try:
            position_window(self.view)
            self.view.show()
        except Exception as e:
            show_error(self.global_settings, "Error showing window", str(e))

    # Add any additional methods specific to this controller
