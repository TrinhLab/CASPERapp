from PyQt6.QtWidgets import QMainWindow, QPushButton, QRadioButton, QComboBox, QPlainTextEdit, QProgressBar
from PyQt6.QtGui import QIcon
from PyQt6 import uic, QtWidgets, QtCore
from utils.ui import scale_ui, show_error
import os
from typing import Optional

class WindowView(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._init_ui()

    def _init_ui(self) -> None:
        """Initialize the user interface."""
        try:
            self._load_ui_file()
            self._init_window_properties()
            self._init_ui_elements()
            self._init_custom_elements()
            self._scale_ui()
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in WindowView", str(e))

    def _load_ui_file(self) -> None:
        """Load the UI file for this window."""
        ui_file = os.path.join(self.global_settings.get_ui_dir(), "name_window.ui")
        uic.loadUi(ui_file, self)

    def _init_window_properties(self) -> None:
        """Set window properties like title and icon."""
        self.setWindowTitle("Window Title")
        self.setWindowIcon(QIcon(os.path.join(self.global_settings.get_assets_dir(), "icon_name.ico")))

    def _init_ui_elements(self) -> None:
        """Initialize UI elements from the loaded UI file."""
        # Example element initializations
        self.push_button_action = self._find_widget("pbtnAction", QPushButton)
        self.combo_box_selection = self._find_widget("cbmSelection", QComboBox)
        self.text_edit_input = self._find_widget("txtedInput", QPlainTextEdit)
        self.progress_bar = self._find_widget("progBar", QProgressBar)

    def _init_custom_elements(self) -> None:
        """Initialize any custom UI elements not in the UI file."""
        # Example custom element
        self.custom_button = QPushButton("Custom Action", self)
        self.custom_button.setGeometry(10, 10, 100, 30)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        """Find a widget by name and type."""
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.global_settings.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def _scale_ui(self) -> None:
        """Scale the UI based on the screen size."""
        scale_ui(self, custom_scale_width=1000, custom_scale_height=600)

    def update_combo_box(self, items: list) -> None:
        """Update items in a combo box."""
        if self.combo_box_selection:
            self.combo_box_selection.clear()
            self.combo_box_selection.addItems(items)

    def set_progress(self, value: int) -> None:
        """Set the value of the progress bar."""
        if self.progress_bar:
            self.progress_bar.setValue(value)

    def reset_progress(self) -> None:
        """Reset the progress bar to 0."""
        self.set_progress(0)

    def show_message(self, message: str) -> None:
        """Display a message to the user."""
        QtWidgets.QMessageBox.information(self, "Information", message)

    def get_input_text(self) -> str:
        """Get the text from the input text edit."""
        return self.text_edit_input.toPlainText() if self.text_edit_input else ""

    def clear_input(self) -> None:
        """Clear the input text edit."""
        if self.text_edit_input:
            self.text_edit_input.clear()

    def enable_action_button(self, enable: bool) -> None:
        """Enable or disable the action button."""
        if self.push_button_action:
            self.push_button_action.setEnabled(enable)

    def resizeEvent(self, event: QtCore.QEvent) -> None:
        """Handle resize events for the window."""
        super().resizeEvent(event)
        # Add any custom resize logic here

    # Add any additional methods specific to this view
