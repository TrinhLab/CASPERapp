from PyQt6 import QtWidgets, QtGui, QtCore, uic
import os
from typing import Optional
from PyQt6.QtCore import Qt

class NewGenomeWindowView(QtWidgets.QMainWindow): 
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.check_mark = "✓" 
        self.spinner_movie = QtGui.QMovie(os.path.join(self.global_settings.get_ui_dir_path(), "throbber.gif"))
        self.spinner_movie.setScaledSize(QtCore.QSize(20, 20)) 

        self._init_ui()

    def _init_ui(self):
        uic.loadUi(os.path.join(self.global_settings.get_ui_dir_path(), 'new_genome_window.ui'), self)
        # self.set_styles()

        self._init_ui_components()

    def set_styles(self):
        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        font: bold 14pt 'Arial';
                        margin-top: 10px;}"""

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1","Step2"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1","Step3"))

    def _init_ui_components(self):
        self._init_grpStep1()
        self._init_grpStep2()
        self._init_grpStep3()

        self.push_button_reset_form = self._find_widget('pbtnResetForm', QtWidgets.QPushButton)

    def _init_grpStep1(self):
        # left column 
        self.line_edit_organism_name = self._find_widget('ledOrganismName', QtWidgets.QLineEdit)
        self.line_edit_strain = self._find_widget('ledStrain', QtWidgets.QLineEdit)
        self.combo_box_endonuclease = self._find_widget('cmbEndonuclease', QtWidgets.QComboBox)
        self.line_edit_organism_code = self._find_widget('ledOrganismCode', QtWidgets.QLineEdit)

        reg_exp_organism_name_strain = QtCore.QRegularExpression("[^/\\\\_]+") 
        reg_exp_organism_code = QtCore.QRegularExpression("\\S+")

        self.line_edit_organism_name.setValidator(QtGui.QRegularExpressionValidator(reg_exp_organism_name_strain, self))
        self.line_edit_strain.setValidator(QtGui.QRegularExpressionValidator(reg_exp_organism_name_strain, self))
        self.line_edit_organism_code.setValidator(QtGui.QRegularExpressionValidator(reg_exp_organism_code, self))

        self.line_edit_organism_name.setFocus()

        # right column
        self.line_edit_seed_length = self._find_widget('ledSeedLength', QtWidgets.QLineEdit)
        self.line_edit_five_length = self._find_widget('ledFiveLength', QtWidgets.QLineEdit)
        self.line_edit_three_length = self._find_widget('ledThreeLength', QtWidgets.QLineEdit)
        self.check_box_generate_repeats = self._find_widget('chckGenerateRepeats', QtWidgets.QCheckBox)
        self.check_box_multithreading = self._find_widget('chckMultithreading', QtWidgets.QCheckBox)

        self.line_edit_seed_length.setEnabled(False)
        self.line_edit_five_length.setEnabled(False)
        self.line_edit_three_length.setEnabled(False)
        self.check_box_generate_repeats.setEnabled(False)

    def _init_grpStep2(self):
        self.push_button_ncbi_search = self._find_widget('pbtnNCBISearch', QtWidgets.QPushButton)
        self.push_button_browse_file = self._find_widget('pbtnBrowseFile', QtWidgets.QPushButton)
        self.line_edit_selected_file = self._find_widget('ledSelectedFile', QtWidgets.QLineEdit)
        self.push_button_add_job = self._find_widget('pbtnAddJob', QtWidgets.QPushButton)

        self.line_edit_selected_file.setPlaceholderText("Selected FASTA/FNA File")

    def _init_grpStep3(self):
        self.push_button_remove_job = self._find_widget('pbtnRemoveJob', QtWidgets.QPushButton)
        self.push_button_remove_all_jobs = self._find_widget('pbtnRemoveAllJobs', QtWidgets.QPushButton)
        self.push_button_run_all_jobs = self._find_widget('pbtnRunAllJobs', QtWidgets.QPushButton)
        self.table_widget_jobs = self._find_widget('tableJobs', QtWidgets.QTableWidget)
        self.progress_bar_jobs = self._find_widget('progbarJobs', QtWidgets.QProgressBar)

        # init jobs table
        self.table_widget_jobs.setColumnCount(2)
        self.table_widget_jobs.setHorizontalHeaderLabels(['Job Name', 'Job Status'])
        self.table_widget_jobs.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table_widget_jobs.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_widget_jobs.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.global_settings.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def set_progress_bar_jobs(self, value):
        self.progress_bar_jobs.setValue(round(value))
        if value >= 100:
            self.progress_bar_jobs.setFormat("All Jobs Completed")
        else:
            self.progress_bar_jobs.setFormat(f"{value:.1f}%")

    def reset_progress_bar_jobs(self):
        self.progress_bar_jobs.setValue(0)
        self.progress_bar_jobs.setFormat("%p%")

    def reset_table_widget_jobs(self):
        self.table_widget_jobs.clearContents()
        self.table_widget_jobs.setRowCount(0)

    def add_job_to_table(self, name):
        row_position = self.table_widget_jobs.rowCount()
        self.table_widget_jobs.insertRow(row_position)
        
        # Job Name
        name_item = QtWidgets.QTableWidgetItem(name)
        name_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter) 
        self.table_widget_jobs.setItem(row_position, 0, name_item)
        
        # Job Status
        status_item = QtWidgets.QTableWidgetItem("Queued")
        status_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.table_widget_jobs.setItem(row_position, 1, status_item)

    def start_spinner(self, row):
        spinner_label = QtWidgets.QLabel()
        spinner_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        spinner_label.setMovie(self.spinner_movie)
        self.spinner_movie.start()
        
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(widget)
        layout.addWidget(QtWidgets.QLabel("Processing"))
        layout.addWidget(spinner_label)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.table_widget_jobs.setCellWidget(row, 1, widget)

    def stop_spinner(self, row):
        self.spinner_movie.stop()
        self.table_widget_jobs.removeCellWidget(row, 1)
        queued_item = QtWidgets.QTableWidgetItem("Queued")
        queued_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.table_widget_jobs.setItem(row, 1, queued_item)

    def set_job_completed(self, row):
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(widget)
        layout.addWidget(QtWidgets.QLabel("Completed"))
        checkmark_label = QtWidgets.QLabel(self.check_mark)
        checkmark_label.setStyleSheet("color: green; font-weight: bold;")
        layout.addWidget(checkmark_label)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.table_widget_jobs.setCellWidget(row, 1, widget)

    def update_endonuclease_dropdown(self, endos):
        self.combo_box_endonuclease.clear()
        self.combo_box_endonuclease.addItems(endos)

    def get_selected_endonuclease(self):
        return self.combo_box_endonuclease.currentText()

    def get_organism_name(self):
        return self.line_edit_organism_name.text()

    def get_strain(self):
        return self.line_edit_strain.text()

    def get_organism_code(self):
        return self.line_edit_organism_code.text()

    def get_selected_file(self):
        return self.line_edit_selected_file.text()

    def set_selected_file(self, file_path):
        # file_name = os.path.basename(file_path)
        self.line_edit_selected_file.setText(file_path)

    def is_multithreading_checked(self):
        return self.check_box_multithreading.isChecked()

    def is_generate_repeats_checked(self):
        return self.check_box_generate_repeats.isChecked()

    def set_endonuclease_lengths(self, seed_length, five_prime_length, three_prime_length):
        self.line_edit_seed_length.setText(seed_length)
        self.line_edit_five_length.setText(five_prime_length)
        self.line_edit_three_length.setText(three_prime_length)

    def isVisible(self):
        return super().isVisible() and not self.isHidden()

    def get_selected_job_identifier(self):
        current_row = self.table_widget_jobs.currentRow()
        if current_row >= 0:
            return self.table_widget_jobs.item(current_row, 0).text()
        return None

    def remove_selected_job_from_table(self):
        current_row = self.table_widget_jobs.currentRow()
        if current_row >= 0:
            self.table_widget_jobs.removeRow(current_row)

    def get_job_name_from_row(self, row):
        item = self.table_widget_jobs.item(row, 0)
        if item:
            return item.text()
        return None

    def sizeHint(self):
        return QtCore.QSize(575, 700)
    