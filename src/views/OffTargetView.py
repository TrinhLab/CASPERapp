from PyQt6 import QtWidgets, uic, QtCore
from utils.ui import show_error

class OffTargetView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        try:
            super().__init__()
            self.global_settings = global_settings
            self.logger = global_settings.get_logger()
            
            uic.loadUi(global_settings.get_ui_dir_path() + '/off_target.ui', self)
            self.setWindowTitle("Off-Target Analysis")
            
            self._init_ui_components()
            
            self.apply_theme()
            
        except Exception as e:
            show_error(self.global_settings, "Error initializing OffTargetView", str(e))

    def _init_ui_components(self):
        """Initialize UI components and set default values"""
        try:
            self._init_grpStep1()
            self._init_grpStep2()
            self._init_grpStep3()
            
            self.push_button_cancel = self.findChild(QtWidgets.QPushButton, 'pbtnCancel')
            self.push_button_submit = self.findChild(QtWidgets.QPushButton, 'pbtnSubmit')
            
            self.setStyleSheet(self.global_settings.get_stylesheet())
            
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI components", str(e))

    def _init_grpStep1(self):
        self.combo_box_organism = self.findChild(QtWidgets.QComboBox, 'cmbOrganism')
        self.combo_box_endonuclease = self.findChild(QtWidgets.QComboBox, 'cmbEndonuclease')

    def _init_grpStep2(self):
        self.double_spin_box_tolerance = self.findChild(QtWidgets.QDoubleSpinBox, 'dspnTolerance')
        self.combo_box_max_mismatches = self.findChild(QtWidgets.QComboBox, 'cmbMaxNoMismatches')
        
        self.radio_button_average_output_yes = self.findChild(QtWidgets.QRadioButton, 'rbtnAverageOutputYes')
        self.radio_button_average_output_no = self.findChild(QtWidgets.QRadioButton, 'rbtnAverageOutputNo')
        
        self.radio_button_save_output_yes = self.findChild(QtWidgets.QRadioButton, 'rbtnSaveOutputFileYes')
        self.radio_button_save_output_no = self.findChild(QtWidgets.QRadioButton, 'rbtnSaveOutputFileNo')
        
        self.average_output_group = QtWidgets.QButtonGroup(self)
        self.average_output_group.addButton(self.radio_button_average_output_yes)
        self.average_output_group.addButton(self.radio_button_average_output_no)
        
        self.save_output_group = QtWidgets.QButtonGroup(self)
        self.save_output_group.addButton(self.radio_button_save_output_yes)
        self.save_output_group.addButton(self.radio_button_save_output_no)
        
        self.radio_button_average_output_no.setChecked(True)
        self.radio_button_save_output_no.setChecked(True)
        
        self.line_edit_output_file = self.findChild(QtWidgets.QLineEdit, 'ledSaveOutputFile')
        self.radio_button_save_output_yes.toggled.connect(self._on_save_output_toggled)
        self.line_edit_output_file.setEnabled(False)  # Initially disabled

    def _init_grpStep3(self):
        """Initialize progress bar only"""
        self.prog_bar = self.findChild(QtWidgets.QProgressBar, 'progBar')

        self.prog_bar.setMinimum(0)
        self.prog_bar.setMaximum(100)
        self.prog_bar.setValue(0)

    def _on_save_output_toggled(self, checked):
        self.line_edit_output_file.setEnabled(checked)
        if not checked:
            self.line_edit_output_file.clear()

    def apply_theme(self):
        """Apply the current theme"""
        if self.global_settings.get_theme() == "dark":
            self.setStyleSheet(self.global_settings.get_dark_stylesheet())
        else:
            self.setStyleSheet(self.global_settings.get_light_stylesheet())

    def get_analysis_parameters(self):
        """Get all parameters needed for off-target analysis"""
        output_filename = self.line_edit_output_file.text().strip()
        save_output = self.radio_button_save_output_yes.isChecked()
        
        return {
            'organism': self.combo_box_organism.currentText(),
            'endonuclease': self.combo_box_endonuclease.currentText(),
            'max_mismatches': int(self.combo_box_max_mismatches.currentText()),
            'tolerance': self.double_spin_box_tolerance.value(),
            'average_output': self.radio_button_average_output_yes.isChecked(),
            'save_output': save_output,
            'output_filename': output_filename if save_output else ''
        }
    
    def set_combo_box_organism(self, organisms):
        self.combo_box_organism.clear()
        self.combo_box_organism.addItems(organisms)

    def set_combo_box_endonuclease(self, endonucleases):
        self.combo_box_endonuclease.clear()
        self.combo_box_endonuclease.addItems(endonucleases)

    def set_combo_box_max_mismatches(self, max_mismatches=10):
        self.combo_box_max_mismatches.clear()
        self.combo_box_max_mismatches.addItems([str(mismatch) for mismatch in range(max_mismatches)])
        self.combo_box_max_mismatches.setCurrentIndex(3) 

    def update_progress_bar(self, value, status=""):
        self.prog_bar.setValue(value)
        if status:
            self.logger.debug(f"Progress: {status}")

    def show_error(self, title, message):
        QtWidgets.QMessageBox.critical(self, title, message)

    def show_warning(self, title, message):
        QtWidgets.QMessageBox.warning(self, title, message)

    def closeEvent(self, event):
        try:
            self.push_button_cancel.clicked.emit()
            event.accept()
        except Exception as e:
            self.logger.error(f"Error in closeEvent: {str(e)}")
            event.accept()

    def set_endonucleases(self, endonucleases):
        """Set available endonucleases in combo box"""
        try:
            self.combo_box_endonuclease.clear()
            self.combo_box_endonuclease.addItems(endonucleases)
        except Exception as e:
            self.logger.error(f"Error setting endonucleases: {str(e)}")
