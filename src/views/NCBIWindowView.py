from PyQt6 import QtWidgets, QtGui, QtCore, uic
from PyQt6.QtWidgets import QWidget, QVBoxLayout
import os
from typing import Optional

class NCBIWindowView(QtWidgets.QMainWindow):
    def __init__(self, settings):
        super(NCBIWindowView, self).__init__()
        self.settings = settings
        self.logger = settings.get_logger()
        self.progress_bars = {}
        self.progress_labels = {}

        self.setup_ui()

    def setup_ui(self):
        uic.loadUi(os.path.join(self.settings.get_ui_dir_path(), "ncbi_window_v2.ui"), self)

        self._init_ui_components()

    def _init_ui_components(self) -> None:
        self._init_grpStep1()
        self._init_grpStep2()
        self._init_grpStep3()

    def _init_grpStep1(self) -> None:
        self.line_edit_organism = self._find_widget("ledOrganism", QtWidgets.QLineEdit)
        self.line_edit_strain = self._find_widget("ledStrain", QtWidgets.QLineEdit)
        self.line_edit_max_results = self._find_widget("ledMaxResults", QtWidgets.QLineEdit)
        self.check_box_complete_genomes_only = self._find_widget("chkCompleteGenomesOnly", QtWidgets.QCheckBox)

    def _init_grpStep2(self) -> None:
        self.push_button_search = self._find_widget("pbtnSearch", QtWidgets.QPushButton)
        self.check_box_select_all_rows = self._find_widget("chkSelectAllRows", QtWidgets.QCheckBox)
        self.table_ncbi_results = self._find_widget("tblNCBIResults", QtWidgets.QTableView)

    def _init_grpStep3(self) -> None:
        self.radio_button_collections_refseq = self._find_widget("rbtnCollectionsRefSeq", QtWidgets.QRadioButton)
        self.radio_button_collections_genbank = self._find_widget("rbtnCollectionsGenBank", QtWidgets.QRadioButton)

        self.check_box_file_types_fna = self._find_widget("chkFileTypesFNA", QtWidgets.QCheckBox)
        self.check_box_file_types_gbff = self._find_widget("chkFileTypesGBFF", QtWidgets.QCheckBox)

        self.push_button_download_files = self._find_widget("pbtnDownloadFiles", QtWidgets.QPushButton)

        self.progress_bar_download_files = self._find_widget("pbDownloadFiles", QtWidgets.QProgressBar)
        self.label_download_files_status = self._find_widget("lblDownloadFilesStatus", QtWidgets.QLabel)

        self.progress_bar_download_files.setValue(0)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.settings.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

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

    def populate_ncbi_table(self, model):
        self.table_ncbi_results.setModel(model)
        
        # Set selection behavior to select entire rows
        self.table_ncbi_results.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_ncbi_results.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
        
        # Set the horizontal header to resize mode
        header = self.table_ncbi_results.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        
        # Enable horizontal scrolling
        self.table_ncbi_results.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table_ncbi_results.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Adjust row height to fit content
        self.table_ncbi_results.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        
        # Center-align the header text
        self.table_ncbi_results.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        
        # Set word wrap for all cells
        self.table_ncbi_results.setWordWrap(True)
        
        # Adjust the size of the table to fit its contents
        self.table_ncbi_results.resizeColumnsToContents()
        self.table_ncbi_results.resizeRowsToContents()
        
        # Set a minimum width for each column based on header text
        font_metrics = QtGui.QFontMetrics(self.table_ncbi_results.font())
        for i in range(self.table_ncbi_results.model().columnCount()):
            header_text = self.table_ncbi_results.model().headerData(i, QtCore.Qt.Orientation.Horizontal)
            width = font_metrics.horizontalAdvance(header_text) + 20  # Add some padding
            self.table_ncbi_results.setColumnWidth(i, max(width, self.table_ncbi_results.columnWidth(i)))
        
        # Ensure the last column doesn't stretch
        header.setStretchLastSection(False)

    def get_search_parameters(self):
        return {
            'organism': self.line_edit_organism.text(),
            'strain': self.line_edit_strain.text(),
            'max_results': self.line_edit_max_results.text(),
            'complete_genomes_only': self.check_box_complete_genomes_only.isChecked(),
            'refseq': self.radio_button_collections_refseq.isChecked(),
            'genbank': self.radio_button_collections_genbank.isChecked(),
            'fna': self.check_box_file_types_fna.isChecked(),
            'gbff': self.check_box_file_types_gbff.isChecked()
        }

    def get_selected_rows(self):
        return self.table_ncbi_results.selectionModel().selectedRows()

    def reset_progress(self):
        self.progress_bar_download_files.setValue(0)
        self.set_download_files_status_label("Ready to download")

    def set_progress(self, value):
        self.progress_bar_download_files.setValue(value)

    def set_download_files_status_label(self, text):
        self.label_download_files_status.setText(text)

    # def clear_form(self):
    #     self.organism_line_edit.clear()
    #     self.infra_name_line_edit.clear()
    #     self.ret_max_line_edit.setText("100")
    #     self.yes_box.setChecked(False)
    #     self.refseq_checkbox.setChecked(True)
    #     self.genbank_checkbox.setChecked(False)
    #     self.fna_checkbox.setChecked(True)
    #     self.gbff_checkbox.setChecked(True)

    def clear_table(self):
        source_model = self.table_ncbi_results.model().sourceModel()
        if source_model:
            source_model.clear()
        self.table_ncbi_results.reset()
