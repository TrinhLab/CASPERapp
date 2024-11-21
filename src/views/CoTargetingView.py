from typing import Optional
from PyQt6 import QtWidgets, uic
from PyQt6.QtWidgets import QTableWidgetItem, QAbstractItemView, QMessageBox
from utils.ui import show_error

class CoTargetingView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/cotargeting.ui', self)
            self.setWindowTitle("Co-targeting")


            self.line_edit_organism = self._find_widget('ledOrganism', QtWidgets.QLineEdit)
            self.table_endonucleases = self._find_widget('tblEndonucleases', QtWidgets.QTableWidget)

            self.push_button_cancel = self._find_widget('pbtnCancel', QtWidgets.QPushButton)
            self.push_button_submit = self._find_widget('pbtnSubmit', QtWidgets.QPushButton)
            
            # Initialize table
            self.table_endonucleases.setColumnCount(1)
            self.table_endonucleases.setShowGrid(True)
            self.table_endonucleases.setHorizontalHeaderLabels(["Endonuclease"])
            self.table_endonucleases.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
            self.table_endonucleases.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            self.table_endonucleases.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
            self.table_endonucleases.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
            
        except Exception as e:
            show_error(self.settings, "Error initializing CoTargeting UI", str(e))

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget 

    def populate_table(self, endo_choices):
        """Populate table with endonuclease choices"""
        try:
            # Filter original endonucleases (no co-targeted ones)
            filtered_endos = [item for item in endo_choices 
                            if len(item.split(",")) == 1 and "|" not in item]
            
            self.table_endonucleases.setRowCount(len(filtered_endos))
            for i, endo in enumerate(filtered_endos):
                self.table_endonucleases.setItem(i, 0, QTableWidgetItem(endo))
                
            self.table_endonucleases.resizeColumnsToContents()
            
        except Exception as e:
            show_error(self.settings, "Error populating table", str(e))

    def get_selected_endonucleases(self):
        """Get list of selected endonucleases"""
        try:
            selected = []
            for item in self.table_endonucleases.selectedItems():
                selected.append(item.text())
            return selected
        except Exception as e:
            self.logger.error(f"Error getting selected endonucleases: {str(e)}")
            return []

    def show_error(self, title, message):
        """Show error message"""
        QMessageBox.critical(self, title, message)

