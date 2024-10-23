from typing import Optional
from PyQt6 import QtWidgets, QtGui, QtCore, uic
from utils.ui import show_error
import os

class NewEndonucleaseView(QtWidgets.QMainWindow):
    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        self.logger = self.settings.get_logger()

        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(os.path.join(self.settings.get_ui_dir_path(), 'new_endonuclease_window.ui'), self)

            self._init_ui_components()
            self.disable_form_elements()
        except Exception as e:
            show_error(self.settings, "Error initializing NewEndonucleaseView", str(e))

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
