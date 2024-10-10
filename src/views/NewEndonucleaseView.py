from PyQt6 import QtWidgets, QtGui, QtCore, uic
import os

class NewEndonucleaseView(QtWidgets.QWidget):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(os.path.join(self.global_settings.get_ui_dir(), 'new_endonuclease.ui'), self)
            self._init_ui_elements()
        except Exception as e:
            self.logger.error(f"Error initializing NewEndonucleaseView: {str(e)}", exc_info=True)
            raise

    def _init_ui_elements(self):
        self._init_grpCasDetails()
        # Initialize other UI elements...

    def _init_grpCasDetails(self):
        self.line_edit_organism = self.findChild(QtWidgets.QLineEdit, 'ledOrganism')
        self.line_edit_abbreviation = self.findChild(QtWidgets.QLineEdit, 'ledAbbreviation')
        self.line_edit_crispr_type = self.findChild(QtWidgets.QLineEdit, 'ledCRISPRType')

        if self.line_edit_organism is None:
            self.logger.error("ledOrganism not found in UI file")
            return
        if self.line_edit_abbreviation is None:
            self.logger.error("ledAbbreviation not found in UI file")
            return
        if self.line_edit_crispr_type is None:
            self.logger.error("ledCRISPRType not found in UI file")
            return

        reg_exp_organism = QtCore.QRegularExpression(r"^[a-zA-Z\s]+$")
        reg_exp_abbreviation = QtCore.QRegularExpression(r"^[a-zA-Z]+$")
        reg_exp_crispr_type = QtCore.QRegularExpression(r"^[a-zA-Z0-9\s]+$")

        self.line_edit_organism.setValidator(QtGui.QRegularExpressionValidator(reg_exp_organism, self))
        self.line_edit_abbreviation.setValidator(QtGui.QRegularExpressionValidator(reg_exp_abbreviation, self))
        self.line_edit_crispr_type.setValidator(QtGui.QRegularExpressionValidator(reg_exp_crispr_type, self))

    # Add other initialization methods as needed