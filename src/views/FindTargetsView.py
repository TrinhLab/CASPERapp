from PyQt6 import QtWidgets
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, 
                            QPushButton, QHBoxLayout, QLabel, QAbstractItemView)
from PyQt6 import uic
from PyQt6.QtCore import Qt, QTimer
import time

class FindTargetsView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._init_ui()
        self.batch_size = 100  # Number of rows to load at once
        self._all_results = []  # Store all results
        self._loaded_rows = 0   # Track number of loaded rows

    def _init_ui(self):
        uic.loadUi(self.global_settings.get_ui_dir_path() + '/find_targets.ui', self)
        self.results_table = self.findChild(QTableWidget, 'tblTargets')
        
        # Optimize table settings for large datasets
        self.results_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.results_table.setShowGrid(False)
        self.results_table.setAlternatingRowColors(True)
        
        # Enable virtual scrolling mode
        self.results_table.setVerticalScrollMode(QTableWidget.ScrollMode.ScrollPerPixel)
        self.results_table.setHorizontalScrollMode(QTableWidget.ScrollMode.ScrollPerPixel)
        
        # Optimize viewport updates
        self.results_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.results_table.viewport().setProperty("cursor", Qt.CursorShape.ArrowCursor)
        
        # Set table properties for better performance
        self.results_table.setColumnCount(5)  # Reduced from 7 to 5 columns
        headers = [
            "Feature Type", "Chromosome/Scaffold #", "Feature ID/Locus Tag", 
            "Feature Name", "Feature Description"
        ]
        self.results_table.setHorizontalHeaderLabels(headers)
        
        # Set optimized column widths
        column_widths = [100, 150, 150, 150, 300]  # Adjusted widths
        for i, width in enumerate(column_widths):
            self.results_table.setColumnWidth(i, width)
        
        self.results_table.horizontalHeader().setStretchLastSection(True)
        
        # Connect scroll events for virtual scrolling
        self.results_table.verticalScrollBar().valueChanged.connect(self._handle_scroll)
        
        self.push_button_view_targets = self.findChild(QPushButton, 'pbtnViewTargets')

    def _create_table_item(self, text):
        """Optimized item creation"""
        item = QTableWidgetItem(str(text))
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        return item

    def _create_row_items(self, result):
        """Create all items for a row at once"""
        return [
            self._create_table_item(result['feature_type']),
            self._create_table_item(str(result['chromosome'])),
            self._create_table_item(result['feature_id']),
            self._create_table_item(result['feature_name']),
            self._create_table_item(result['feature_description'])
        ]

    def display_results(self, results):
        start_time = time.time()
        
        # Store all results and reset loaded count
        self._all_results = results
        self._loaded_rows = 0
        
        # Disable visual updates
        self.results_table.setUpdatesEnabled(False)
        self.results_table.setSortingEnabled(False)
        self.results_table.setVisible(False)
        
        # Set total row count
        total_rows = len(results)
        self.results_table.setRowCount(total_rows)
        
        # Load initial batch
        self._load_batch(0, min(self.batch_size, total_rows))
        
        # Re-enable table and updates
        self.results_table.setVisible(True)
        self.results_table.setUpdatesEnabled(True)
        self.results_table.setSortingEnabled(True)
        
        total_time = time.time() - start_time
        self.global_settings.logger.debug(f"Initial display time: {total_time:.2f} seconds")

    def _load_batch(self, start_idx, end_idx):
        """Load a batch of rows efficiently"""
        if start_idx >= len(self._all_results) or start_idx >= end_idx:
            return
        
        batch_items = []
        for row in range(start_idx, end_idx):
            if row >= len(self._all_results):
                break
            row_items = self._create_row_items(self._all_results[row])
            batch_items.append((row, row_items))
        
        # Batch set items
        for row, items in batch_items:
            for col, item in enumerate(items):
                self.results_table.setItem(row, col, item)
        
        self._loaded_rows = end_idx

    def _handle_scroll(self, value):
        """Handle scroll events for virtual scrolling"""
        viewport_height = self.results_table.viewport().height()
        row_height = self.results_table.rowHeight(0)
        visible_rows = viewport_height // row_height
        
        # Calculate which rows should be visible
        scroll_position = value
        start_row = max(0, scroll_position - visible_rows)
        end_row = min(len(self._all_results), scroll_position + visible_rows * 2)
        
        # Load more rows if needed
        if end_row > self._loaded_rows:
            self._load_batch(self._loaded_rows, end_row)

    def get_selected_targets(self):
        selected_rows = set(index.row() for index in self.results_table.selectedIndexes())
        selected_targets = []
        
        for row in selected_rows:
            if row < len(self._all_results):
                selected_targets.append(self._all_results[row])
        
        return selected_targets

    def clear_results(self):
        """Clear all results from the table"""
        self.results_table.setUpdatesEnabled(False)
        self.results_table.clearContents()
        self.results_table.setRowCount(0)
        self._all_results = []
        self._loaded_rows = 0
        self.results_table.setUpdatesEnabled(True)
