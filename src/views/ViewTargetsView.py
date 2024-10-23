from typing import Optional
from PyQt6 import QtWidgets, uic
from PyQt6.QtWidgets import QTableWidgetItem, QAbstractItemView
from PyQt6.QtGui import QTextDocument
from utils.ui import show_error

class ViewTargetsView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()

        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/view_targets.ui', self)
            self._init_ui_components()
        except Exception as e:
            show_error(self.settings, "Error initializing ViewTargetsView", str(e))

    def _init_ui_components(self):
        self._init_grpGuideViewer()
        self._init_grpGuideAnalysis()
        self._init_grpGeneViewer()

        self.push_button_export_grna = self._find_widget('pbtnExportgRNA', QtWidgets.QPushButton)

    def _init_grpGuideViewer(self):
        self.combo_box_gene = self._find_widget('cmbGene', QtWidgets.QComboBox)
        self.combo_box_endonuclease = self._find_widget('cmbEndonuclease', QtWidgets.QComboBox)
        self.check_box_select_all = self._find_widget('chkSelectAll', QtWidgets.QCheckBox)
        self.push_button_filter_options = self._find_widget('pbtnFilterOptions', QtWidgets.QPushButton)
        self.push_button_scoring_options = self._find_widget('pbtnScoringOptions', QtWidgets.QPushButton)
        self.table_targets = self._find_widget('tblTargets', QtWidgets.QTableWidget)

        self.table_targets.setColumnCount(8)
        self.table_targets.setHorizontalHeaderLabels(["Location", "Endonuclease", "Sequence", "Strand", "PAM", "Score", "Off-Target", "Details"])
        self.table_targets.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_targets.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_targets.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.table_targets.horizontalHeader().setSectionResizeMode(7, QtWidgets.QHeaderView.ResizeMode.Stretch)

    def _init_grpGuideAnalysis(self):
        self.push_button_off_target = self._find_widget('pbtnOffTarget', QtWidgets.QPushButton)
        self.push_button_cotargeting = self._find_widget('pbtnCoTargeting', QtWidgets.QPushButton)

    def _init_grpGeneViewer(self):
        self.push_button_highlight_guides = self._find_widget('pbtnHighlightGuides', QtWidgets.QPushButton)
        self.push_button_clear_guides = self._find_widget('pbtnClearGuides', QtWidgets.QPushButton)
        self.line_edit_start_location = self._find_widget('ledStartLocation', QtWidgets.QLineEdit)
        self.line_edit_stop_location = self._find_widget('ledStopLocation', QtWidgets.QLineEdit)
        self.push_button_change_location = self._find_widget('pbtnChangeLocation', QtWidgets.QPushButton)
        self.text_edit_gene_viewer = self._find_widget('txtedGeneViewer', QtWidgets.QTextEdit)
        self.push_button_reset_location = self._find_widget('pbtnResetLocation', QtWidgets.QPushButton)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget 

    def display_targets_in_table(self, targets):
        self.table_targets.setRowCount(len(targets))
        for row, target in enumerate(targets):
            self.table_targets.setItem(row, 0, QTableWidgetItem(str(target[0])))  # Location
            self.table_targets.setItem(row, 1, QTableWidgetItem(target[5]))  # Endonuclease
            self.table_targets.setItem(row, 2, QTableWidgetItem(target[1]))  # Sequence
            self.table_targets.setItem(row, 3, QTableWidgetItem(target[4]))  # Strand
            self.table_targets.setItem(row, 4, QTableWidgetItem(target[2]))  # PAM
            self.table_targets.setItem(row, 5, QTableWidgetItem(str(target[3])))  # Score
            self.table_targets.setItem(row, 6, QTableWidgetItem("N/A"))  # Off-Target (placeholder)
            
            details_button = QtWidgets.QPushButton("Details")
            self.table_targets.setCellWidget(row, 7, details_button)

        self.table_targets.resizeColumnsToContents()

    def get_selected_targets(self):
        selected_rows = set(index.row() for index in self.table_targets.selectedIndexes())
        return [self.get_row_data(row) for row in selected_rows]

    def get_row_data(self, row):
        return {
            'location': self.table_targets.item(row, 0).text(),
            'endonuclease': self.table_targets.item(row, 1).text(),
            'sequence': self.table_targets.item(row, 2).text(),
            'strand': self.table_targets.item(row, 3).text(),
            'pam': self.table_targets.item(row, 4).text(),
            'score': self.table_targets.item(row, 5).text(),
            'off_target': self.table_targets.item(row, 6).text()
        }
    
    def set_combo_box_endonuclease(self, endonucleases):
        self.combo_box_endonuclease.addItems(endonucleases)

    def set_combo_box_gene(self, genes):
        self.combo_box_gene.addItems(genes)

    def set_text_edit_gene_viewer(self, sequence):
        self.text_edit_gene_viewer.setText(sequence)

    def update_gene_info(self, info):
        # Implement this method if you have a widget to display gene info
        pass

    def update_gene_viewer(self, sequence):
        self.text_edit_gene_viewer.clear()
        doc = QTextDocument()
        doc.setHtml(sequence)
        self.text_edit_gene_viewer.setDocument(doc)

    def select_all_targets(self, select):
        for row in range(self.table_targets.rowCount()):
            self.table_targets.selectRow(row) if select else self.table_targets.clearSelection()

    def show_filter_options_dialog(self, options):
        # Implement this method to show filter options dialog
        pass

    def filter_options_accepted(self):
        # Implement this method to check if filter options were accepted
        return True

    def get_filter_options(self):
        # Implement this method to return new filter options
        return {}

    def show_scoring_options_dialog(self, options):
        # Implement this method to show scoring options dialog
        pass

    def scoring_options_accepted(self):
        # Implement this method to check if scoring options were accepted
        return True

    def get_scoring_options(self):
        # Implement this method to return new scoring options
        return {}

    def get_export_file_path(self):
        # Implement this method to get the export file path from the user
        return QtWidgets.QFileDialog.getSaveFileName(self, 'Save File')[0]
