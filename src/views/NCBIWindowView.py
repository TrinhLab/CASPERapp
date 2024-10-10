from PyQt6 import QtWidgets, QtGui, QtCore, uic
from PyQt6.QtWidgets import QWidget, QVBoxLayout
import os
from typing import Optional

class NCBIWindowView(QWidget):
    def __init__(self, settings):
        super(NCBIWindowView, self).__init__()
        self.settings = settings
        self.setup_ui()
        self.progress_bars = {}
        self.progress_labels = {}

    def setup_ui(self):
        uic.loadUi(self.settings.get_ui_dir() + "/ncbi_window_v2.ui", self)
        print("NCBIWindowView initialized")

        self._init_ui_elements()
        self.progressBar.setValue(0)

    def _init_ui_elements(self) -> None:
        # Create a main layout to hold everything
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Create a widget to hold the original content
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        for child in self.children():
            if isinstance(child, QWidget):
                content_layout.addWidget(child)

        main_layout.addWidget(content_widget)

        self._init_grpStep1()
        self._init_grpStep2()
        self._init_grpStep3()

    def _init_grpStep1(self) -> None:
        self.organism_line_edit = self._find_widget("organism_line_edit", QtWidgets.QLineEdit)
        self.infra_name_line_edit = self._find_widget("infra_name_line_edit", QtWidgets.QLineEdit)
        self.ret_max_line_edit = self._find_widget("ret_max_line_edit", QtWidgets.QLineEdit)
        self.yes_box = self._find_widget("yes_box", QtWidgets.QCheckBox)

    def _init_grpStep2(self) -> None:
        self.search_button = self._find_widget("search_button", QtWidgets.QPushButton)
        self.all_rows = self._find_widget("all_rows", QtWidgets.QCheckBox)
        self.ncbi_table = self._find_widget("ncbi_table", QtWidgets.QTableView)

    def _init_grpStep3(self) -> None:
        self.genbank_checkbox = self._find_widget("genbank_checkbox", QtWidgets.QRadioButton)
        self.refseq_checkbox = self._find_widget("refseq_checkbox", QtWidgets.QRadioButton)
        self.fna_checkbox = self._find_widget("fna_checkbox", QtWidgets.QCheckBox)
        self.gbff_checkbox = self._find_widget("gbff_checkbox", QtWidgets.QCheckBox)
        self.download_button = self._find_widget("download_button", QtWidgets.QPushButton)
        self.progressBar = self._find_widget("progressBar", QtWidgets.QProgressBar)
        self.progressLabel = self._find_widget("progressLabel", QtWidgets.QLabel)

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

    def populate_table(self, model):
        self.ncbi_table.setModel(model)
        
        # Set the horizontal header to resize mode
        header = self.ncbi_table.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        
        # Enable horizontal scrolling
        self.ncbi_table.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.ncbi_table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Adjust row height to fit content
        self.ncbi_table.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        
        # Center-align the header text
        self.ncbi_table.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        
        # Set word wrap for all cells
        self.ncbi_table.setWordWrap(True)
        
        # Adjust the size of the table to fit its contents
        self.ncbi_table.resizeColumnsToContents()
        self.ncbi_table.resizeRowsToContents()
        
        # Set a minimum width for each column based on header text
        font_metrics = QtGui.QFontMetrics(self.ncbi_table.font())
        for i in range(self.ncbi_table.model().columnCount()):
            header_text = self.ncbi_table.model().headerData(i, QtCore.Qt.Orientation.Horizontal)
            width = font_metrics.horizontalAdvance(header_text) + 20  # Add some padding
            self.ncbi_table.setColumnWidth(i, max(width, self.ncbi_table.columnWidth(i)))
        
        # Ensure the last column doesn't stretch
        header.setStretchLastSection(False)
        
        # Connect the selectionChanged signal after populating the table
        if self.ncbi_table.selectionModel():
            self.ncbi_table.selectionModel().selectionChanged.connect(self.on_selection_changed)

    def get_search_parameters(self):
        return {
            'organism': self.organism_line_edit.text(),
            'strain': self.infra_name_line_edit.text(),
            'ret_max': self.ret_max_line_edit.text(),
            'complete_genomes_only': self.yes_box.isChecked(),
            'refseq': self.refseq_checkbox.isChecked(),
            'genbank': self.genbank_checkbox.isChecked(),
            'fna': self.fna_checkbox.isChecked(),
            'gbff': self.gbff_checkbox.isChecked()
        }

    def get_selected_rows(self):
        return self.ncbi_table.selectionModel().selectedRows()

    def reset_progress(self):
        self.progressBar.setValue(0)
        self.statusLabel.setText("Ready to download")

    def set_progress(self, value):
        self.progressBar.setValue(value)

    def update_status_label(self, text):
        self.statusLabel.setText(text)

    # Remove or comment out these methods as they're no longer needed
    # def add_progress_bar(self, id):
    # def update_progress_bar(self, id, value, max_value):
    # def update_progress_label(self, id, text):
    # def clear_progress_bars(self):

    def clear_form(self):
        self.organism_line_edit.clear()
        self.infra_name_line_edit.clear()
        self.ret_max_line_edit.setText("100")
        self.yes_box.setChecked(False)
        self.refseq_checkbox.setChecked(True)
        self.genbank_checkbox.setChecked(False)
        self.fna_checkbox.setChecked(True)
        self.gbff_checkbox.setChecked(True)

    def clear_table(self):
        source_model = self.ncbi_table.model().sourceModel()
        if source_model:
            source_model.clear()
        self.ncbi_table.reset()

    def on_selection_changed(self, selected, deselected):
        # Update the "Select All" checkbox based on the selection state
        if self.ncbi_table.model():
            all_selected = len(self.ncbi_table.selectionModel().selectedRows()) == self.ncbi_table.model().rowCount()
            self.all_rows.setChecked(all_selected)

    def show_window(self):
        print("NCBIWindowView: show_window called")
        self.show()
        if self.search_button:
            self.search_button.setFocus()
            self.search_button.setVisible(True)
            self.search_button.raise_()
            self.search_button.setEnabled(True)
            print(f"NCBIWindowView: Search button has focus: {self.search_button.hasFocus()}")
            print(f"NCBIWindowView: Search button is visible: {self.search_button.isVisible()}")
            print(f"NCBIWindowView: Search button is enabled: {self.search_button.isEnabled()}")
            print(f"NCBIWindowView: Search button parent: {self.search_button.parent()}")
            print(f"NCBIWindowView: Search button geometry: {self.search_button.geometry()}")
        else:
            print("NCBIWindowView: Search button not found")
        
        # Ensure all widgets are properly laid out
        self.adjustSize()
