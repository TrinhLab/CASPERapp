from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, QPushButton, QHBoxLayout, QLabel
from PyQt6 import uic
from PyQt6.QtCore import Qt, QTimer

class FindTargetsView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._init_ui()

    def _init_ui(self):
        uic.loadUi(self.global_settings.get_ui_dir_path() + '/find_targets.ui', self)
        self.results_table = self.findChild(QTableWidget, 'tblTargets')
        self.results_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        
        # Optimize table performance
        self.results_table.setUpdatesEnabled(False)  # Disable updates during setup
        self.results_table.setSortingEnabled(False)  # Disable sorting during setup
        self.results_table.horizontalHeader().setStretchLastSection(True)
        
        # Set up the table columns
        self.results_table.setColumnCount(7)
        self.results_table.setHorizontalHeaderLabels([
            "Feature Type", "Chromosome/Scaffold #", "Feature ID/Locus Tag", 
            "Feature Name", "Feature Description", "Location", "Strand"
        ])

        self.push_button_view_targets = self.findChild(QPushButton, 'pbtnViewTargets')
        
        # Pre-allocate items for better performance
        self._cached_items = {}

    def _get_table_item(self, text):
        """Cache and reuse QTableWidgetItems for better performance"""
        if text not in self._cached_items:
            item = QTableWidgetItem(str(text))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)  # Make item read-only
            self._cached_items[text] = item
        return self._cached_items[text].clone()

    def display_results(self, results):
        # Disable updates for bulk operations
        self.results_table.setUpdatesEnabled(False)
        self.results_table.setSortingEnabled(False)
        
        # Set row count once
        self.results_table.setRowCount(len(results))

        # Batch insert items
        for row, result in enumerate(results):
            self.results_table.setItem(row, 0, self._get_table_item(result['feature_type']))
            self.results_table.setItem(row, 1, self._get_table_item(str(result['chromosome'])))
            self.results_table.setItem(row, 2, self._get_table_item(result['feature_id']))
            self.results_table.setItem(row, 3, self._get_table_item(result['feature_name']))
            self.results_table.setItem(row, 4, self._get_table_item(result['feature_description']))
            self.results_table.setItem(row, 5, self._get_table_item(result['location']))
            self.results_table.setItem(row, 6, self._get_table_item(result['strand']))

        # Re-enable updates and adjust columns
        QTimer.singleShot(0, self._finish_table_update)

    def _finish_table_update(self):
        """Complete table update in the next event loop iteration"""
        self.results_table.resizeColumnsToContents()
        self.results_table.setUpdatesEnabled(True)
        self.results_table.setSortingEnabled(True)

    def get_selected_targets(self):
        selected_rows = set(index.row() for index in self.results_table.selectedIndexes())
        selected_targets = []
        
        # Get column indices once
        columns = {
            'feature_type': 0,
            'chromosome': 1,
            'feature_id': 2,
            'feature_name': 3,
            'feature_description': 4,
            'location': 5,
            'strand': 6
        }
        
        for row in selected_rows:
            target = {
                'feature_type': self.results_table.item(row, columns['feature_type']).text(),
                'chromosome': self.results_table.item(row, columns['chromosome']).text(),
                'feature_id': self.results_table.item(row, columns['feature_id']).text(),
                'feature_name': self.results_table.item(row, columns['feature_name']).text(),
                'feature_description': self.results_table.item(row, columns['feature_description']).text(),
                'location': self.results_table.item(row, columns['location']).text(),
                'strand': self.results_table.item(row, columns['strand']).text()
            }
            selected_targets.append(target)
        return selected_targets
