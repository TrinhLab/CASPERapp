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
            uic.loadUi(os.path.join(self.settings.get_ui_dir(), 'new_endonuclease_window.ui'), self)
            self._init_ui_components()
        except Exception as e:
            show_error(self.settings, "Error initializing NewEndonucleaseView", str(e))

    def _init_ui_components(self):
        self._init_grpCasDetails()
        self._init_grpgRNADetails()
        self._init_grpScoringAlgorithms()
        self._init_boxlayhbotButtons()

    def _init_grpCasDetails(self):
        self.line_edit_organism = self._find_widget('ledOrganism', QtWidgets.QLineEdit)
        self.line_edit_abbreviation = self._find_widget('ledAbbreviation', QtWidgets.QLineEdit)
        self.line_edit_crispr_type = self._find_widget('ledCRISPRType', QtWidgets.QLineEdit)

        reg_exp_organism = QtCore.QRegularExpression(r"^[a-zA-Z\s]+$")
        reg_exp_abbreviation = QtCore.QRegularExpression(r"^[a-zA-Z]+$")
        reg_exp_crispr_type = QtCore.QRegularExpression(r"^[a-zA-Z0-9\s]+$")

        self.line_edit_organism.setValidator(QtGui.QRegularExpressionValidator(reg_exp_organism, self))
        self.line_edit_abbreviation.setValidator(QtGui.QRegularExpressionValidator(reg_exp_abbreviation, self))
        self.line_edit_crispr_type.setValidator(QtGui.QRegularExpressionValidator(reg_exp_crispr_type, self))

    def _init_grpgRNADetails(self):
        self.line_edit_seed_length = self._find_widget('ledSeedLength', QtWidgets.QLineEdit)
        self.line_edit_five_prime_length = self._find_widget('ledFivePrimeLength', QtWidgets.QLineEdit)
        self.line_edit_three_prime_length = self._find_widget('ledThreePrimeLength', QtWidgets.QLineEdit)
        self.line_edit_PAM_sequence = self._find_widget('ledPAMSequence', QtWidgets.QLineEdit)

    def _init_grpScoringAlgorithms(self):
        self.combo_box_on_target_matrix = self._find_widget('cmbOnTargetMatrix', QtWidgets.QComboBox)
        self.combo_box_off_target_matrix = self._find_widget('cmbOffTargetMatrix', QtWidgets.QComboBox)

    def _init_boxlayhbotButtons(self):
        self.push_button_reset_form = self._find_widget('pbtnResetForm', QtWidgets.QPushButton)
        self.push_button_submit = self._find_widget('pbtnSubmit', QtWidgets.QPushButton)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget