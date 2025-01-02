from typing import Optional
from PyQt6 import QtWidgets, QtGui, QtCore, uic
from utils.ui import show_error
import os
import qdarktheme

class NewEndonucleaseView(QtWidgets.QMainWindow):
    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        self.logger = self.settings.get_logger()
        
        # Set window properties
        self.setWindowTitle("New Endonuclease")
        self.setMinimumSize(400, 500)  # Set minimum window size
        
        # Center the window on screen
        screen = QtGui.QGuiApplication.primaryScreen()
        screen_geometry = screen.geometry()
        centerPoint = screen_geometry.center()
        
        self.init_ui()
        
        # Calculate and set position to center
        frame_geometry = self.frameGeometry()
        frame_geometry.moveCenter(centerPoint)
        self.move(frame_geometry.topLeft())

        # Apply theme
        self.apply_theme()

    def apply_theme(self):
        current_theme = self.settings.get_theme()
        themes = {
            "dark": {
                "bg_color": "#2b2b2b",
                "fg_color": "#ffffff",
                "button_bg_color": "#3a3a3a",
                "button_border_color": "#5a5a5a",
                "button_hover_bg_color": "#4a4a4a",
                "input_bg_color": "#3a3a3a",
                "input_border_color": "#5a5a5a",
                "menu_bg_color": "#2b2b2b",
                "menu_item_hover_bg_color": "#3a3a3a",
                "divider_color": "#444444"
            },
            "light": {
                "bg_color": "#f0f0f0",
                "fg_color": "#000000",
                "button_bg_color": "#e0e0e0",
                "button_border_color": "#c0c0c0",
                "button_hover_bg_color": "#d0d0d0",
                "input_bg_color": "#ffffff",
                "input_border_color": "#c0c0c0",
                "menu_bg_color": "#f0f0f0",
                "menu_item_hover_bg_color": "#e0e0e0",
                "divider_color": "#c0c0c0"
            }
        }

        theme = themes["dark"] if current_theme == "dark" else themes["light"]
        qdarktheme.setup_theme(current_theme)

        self.setStyleSheet(f"""
            QMainWindow {{
                background-color: {theme['bg_color']};
                color: {theme['fg_color']};
            }}
            QWidget {{
                background-color: {theme['bg_color']};
                color: {theme['fg_color']};
            }}
            QPushButton {{
                background-color: {theme['button_bg_color']};
                border: 1px solid {theme['button_border_color']};
                padding: 5px;
                min-width: 80px;
            }}
            QPushButton:hover {{
                background-color: {theme['button_hover_bg_color']};
            }}
            QLineEdit {{
                background-color: {theme['input_bg_color']};
                border: 1px solid {theme['input_border_color']};
                padding: 5px;
            }}
            QComboBox {{
                background-color: {theme['input_bg_color']};
                border: 1px solid {theme['input_border_color']};
                padding: 5px;
            }}
            QComboBox:hover {{
                background-color: {theme['button_hover_bg_color']};
            }}
            QComboBox::drop-down {{
                border: none;
            }}
            QComboBox::down-arrow {{
                image: none;
                border: none;
            }}
            QRadioButton {{
                color: {theme['fg_color']};
            }}
            QRadioButton::indicator {{
                width: 13px;
                height: 13px;
            }}
            QRadioButton::indicator:checked {{
                background-color: {theme['button_hover_bg_color']};
                border: 2px solid {theme['button_border_color']};
                border-radius: 7px;
            }}
            QRadioButton::indicator:unchecked {{
                background-color: {theme['bg_color']};
                border: 2px solid {theme['button_border_color']};
                border-radius: 7px;
            }}
            QLabel {{
                color: {theme['fg_color']};
            }}
        """)

    def showEvent(self, event):
        """Override showEvent to apply theme and styles when window is shown"""
        super().showEvent(event)
        self._set_styles()

    def init_ui(self):
        try:
            uic.loadUi(os.path.join(self.settings.get_ui_dir_path(), 'new_endonuclease_window.ui'), self)
            self._init_ui_components()
            self._set_styles()  # Add style initialization
            self.disable_form_elements()
        except Exception as e:
            show_error(self.settings, "Error initializing NewEndonucleaseView", str(e))

    def _set_styles(self):
        """Apply the global groupbox style and theme"""
        try:
            # Apply global groupbox style
            style = self.settings.get_groupbox_style()
            for groupbox in self.findChildren(QtWidgets.QGroupBox):
                groupbox.setStyleSheet(style)
            
            # Apply theme
            self.apply_theme()
        except Exception as e:
            self.logger.error(f"Error setting styles: {str(e)}")

    def _init_ui_components(self):
        self.combo_box_select_endonuclease = self._find_widget('cmbSelectEndonuclease', QtWidgets.QComboBox)

        self._init_grpCasDetails()
        self._init_grpgRNADetails()
        self._init_grpScoringAlgorithms()
        self._init_boxlayhbotButtons()

    def _init_grpCasDetails(self):
        self.line_edit_organism = self._find_widget('ledOrganism', QtWidgets.QLineEdit)
        self.line_edit_abbreviation = self._find_widget('ledAbbreviation', QtWidgets.QLineEdit)
        self.line_edit_CRISPR_type = self._find_widget('ledCRISPRType', QtWidgets.QLineEdit)

        reg_exp_organism = QtCore.QRegularExpression(r"^[a-zA-Z\s]+$")
        reg_exp_abbreviation = QtCore.QRegularExpression(r"^[a-zA-Z]+$")
        reg_exp_CRISPR_type = QtCore.QRegularExpression(r"^[a-zA-Z0-9\s-]+$")

        self.line_edit_organism.setValidator(QtGui.QRegularExpressionValidator(reg_exp_organism, self))
        self.line_edit_abbreviation.setValidator(QtGui.QRegularExpressionValidator(reg_exp_abbreviation, self))
        self.line_edit_CRISPR_type.setValidator(QtGui.QRegularExpressionValidator(reg_exp_CRISPR_type, self))

    def _init_grpgRNADetails(self):
        self.line_edit_seed_length = self._find_widget('ledSeedLength', QtWidgets.QLineEdit)
        self.line_edit_five_prime_length = self._find_widget('ledFivePrimeLength', QtWidgets.QLineEdit)
        self.line_edit_three_prime_length = self._find_widget('ledThreePrimeLength', QtWidgets.QLineEdit)
        self.line_edit_PAM_sequence = self._find_widget('ledPAMSequence', QtWidgets.QLineEdit)
        self.radio_button_three_prime_pam = self._find_widget('rbtnThreePrimePAM', QtWidgets.QRadioButton)

    def _init_grpScoringAlgorithms(self):
        self.combo_box_on_target_matrix = self._find_widget('cmbOnTargetMatrix', QtWidgets.QComboBox)
        self.combo_box_off_target_matrix = self._find_widget('cmbOffTargetMatrix', QtWidgets.QComboBox)

    def _init_boxlayhbotButtons(self):
        self.push_button_reset_form = self._find_widget('pbtnResetForm', QtWidgets.QPushButton)
        self.push_button_delete = self._find_widget('pbtnDelete', QtWidgets.QPushButton)
        self.push_button_submit = self._find_widget('pbtnSubmit', QtWidgets.QPushButton)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget
    
    def disable_form_elements(self):
        for widget in self.findChildren((QtWidgets.QLineEdit, QtWidgets.QComboBox, QtWidgets.QPushButton)):
            if widget != self.combo_box_select_endonuclease:
                widget.setEnabled(False)
                if isinstance(widget, QtWidgets.QLineEdit):
                    widget.setStyleSheet("background-color: #f0f0f0;")

    def enable_form_elements(self):
        for widget in self.findChildren((QtWidgets.QLineEdit, QtWidgets.QComboBox, QtWidgets.QPushButton)):
            widget.setEnabled(True)
            if isinstance(widget, QtWidgets.QLineEdit):
                widget.setStyleSheet("")
    
    def enable_delete_button(self):
        self.push_button_delete.setEnabled(True)
        self.push_button_delete.setVisible(True)

    def disable_delete_button(self):
        self.push_button_delete.setEnabled(False)
        self.push_button_delete.setVisible(False)

    def update_combo_box_on_target_matrix(self, on_list):
        self.combo_box_on_target_matrix.clear()
        self.combo_box_on_target_matrix.addItems(on_list)

    def update_combo_box_off_target_matrix(self, off_list):
        self.combo_box_off_target_matrix.clear()
        self.combo_box_off_target_matrix.addItems(off_list)

    def get_form_data(self):
        return {
            'endonuclease_organism': self.line_edit_organism.text(),
            'endonuclease_abbreviation': self.line_edit_abbreviation.text(),
            'endonuclease_CRISPR_type': self.line_edit_CRISPR_type.text(),
            'endonuclease_seed_length': self.line_edit_seed_length.text(),
            'endonuclease_five_prime_length': self.line_edit_five_prime_length.text(),
            'endonuclease_three_prime_length': self.line_edit_three_prime_length.text(),
            'endonuclease_pam_sequence': self.line_edit_PAM_sequence.text(),
            'endonuclease_direction': '3' if self.radio_button_three_prime_pam.isChecked() else '5',
            'endonuclease_on_target_scoring': self.combo_box_on_target_matrix.currentText(),
            'endonuclease_off_target_scoring': self.combo_box_off_target_matrix.currentText()
        }

    def clear_form(self):
        self.line_edit_organism.clear()
        self.line_edit_abbreviation.clear()
        self.line_edit_CRISPR_type.clear()
        self.line_edit_seed_length.clear()
        self.line_edit_five_prime_length.clear()
        self.line_edit_three_prime_length.clear()
        self.line_edit_PAM_sequence.clear()
        self.combo_box_on_target_matrix.setCurrentIndex(0)
        self.combo_box_off_target_matrix.setCurrentIndex(0)

    def populate_form(self, data):
        self.line_edit_organism.setText(data['endonuclease_organism'])
        self.line_edit_abbreviation.setText(data['endonuclease_abbreviation'])
        self.line_edit_CRISPR_type.setText(data['endonuclease_CRISPR_type'])
        self.line_edit_seed_length.setText(data['endonuclease_seed_length'])
        self.line_edit_five_prime_length.setText(data['endonuclease_five_prime_length'])
        self.line_edit_three_prime_length.setText(data['endonuclease_three_prime_length'])
        self.line_edit_PAM_sequence.setText(data['endonuclease_pam_sequence'])
        self.combo_box_on_target_matrix.setCurrentText(data['endonuclease_on_target_scoring'])
        self.combo_box_off_target_matrix.setCurrentText(data['endonuclease_off_target_scoring'])

    def update_combo_box_select_endonuclease(self, endonucleases):
        self.combo_box_select_endonuclease.clear()
        self.combo_box_select_endonuclease.addItem("")
        self.combo_box_select_endonuclease.addItem("Define New Endonuclease")
        self.combo_box_select_endonuclease.addItems(endonucleases)

    def get_selected_endonuclease(self):
        return self.combo_box_select_endonuclease.currentText()
