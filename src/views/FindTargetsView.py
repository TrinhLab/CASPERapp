from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, QPushButton, QHBoxLayout, QLabel
from PyQt6 import uic
from PyQt6.QtCore import Qt

class FindTargetsView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._init_ui()

    def _init_ui(self):
        uic.loadUi(self.global_settings.get_ui_dir_path() + '/find_targets.ui', self)
        self.results_table = self.findChild(QTableWidget, 'tblTargets')
        self.results_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        
        # Set up the table columns
        self.results_table.setColumnCount(7)
        self.results_table.setHorizontalHeaderLabels([
            "Feature Type", "Chromosome/Scaffold #", "Feature ID/Locus Tag", 
            "Feature Name", "Feature Description", "Location", "Strand"
        ])

        self.push_button_view_targets = self.findChild(QPushButton, 'pbtnViewTargets')
        
        # Add buttons
        # button_layout = QHBoxLayout()
        # self.export_button = QPushButton("Export")
        # self.view_targets_button = QPushButton("View Targets")
        # self.generate_library_button = QPushButton("Generate Library")
        # button_layout.addWidget(self.export_button)
        # button_layout.addWidget(self.view_targets_button)
        # button_layout.addWidget(self.generate_library_button)
        
        # Add result count label
        # self.result_count_label = QLabel("Results: 0")
        
        # main_layout = QVBoxLayout(self)
        # main_layout.addWidget(self.result_count_label)
        # main_layout.addWidget(self.results_table)
        # main_layout.addLayout(button_layout)

    def display_results(self, results):
        self.results_table.setRowCount(len(results))

        for row, result in enumerate(results):
            self.results_table.setItem(row, 0, QTableWidgetItem(result['feature_type']))
            self.results_table.setItem(row, 1, QTableWidgetItem(str(result['chromosome'])))
            self.results_table.setItem(row, 2, QTableWidgetItem(result['feature_id']))
            self.results_table.setItem(row, 3, QTableWidgetItem(result['feature_name']))
            self.results_table.setItem(row, 4, QTableWidgetItem(result['feature_description']))
            self.results_table.setItem(row, 5, QTableWidgetItem(result['location']))
            self.results_table.setItem(row, 6, QTableWidgetItem(result['strand']))

        self.results_table.resizeColumnsToContents()
        # self.result_count_label.setText(f"Results: {len(results)}")

    def get_selected_targets(self):
        selected_rows = set(index.row() for index in self.results_table.selectedIndexes())
        selected_targets = []
        for row in selected_rows:
            target = {
                'feature_type': self.results_table.item(row, 0).text(),
                'chromosome': self.results_table.item(row, 1).text(),
                'feature_id': self.results_table.item(row, 2).text(),
                'feature_name': self.results_table.item(row, 3).text(),
                'feature_description': self.results_table.item(row, 4).text(),
                'location': self.results_table.item(row, 5).text(),
                'strand': self.results_table.item(row, 6).text()
            }
            selected_targets.append(target)
        return selected_targets
