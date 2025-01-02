from PyQt6 import QtWidgets
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, 
                            QPushButton, QHBoxLayout, QLabel, QAbstractItemView, QCheckBox)
from PyQt6 import uic
from PyQt6.QtCore import Qt, QTimer
import time

class FindTargetsView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.logger
        self._init_ui()
        self.batch_size = 100  # Number of rows to load at once
        self._all_results = []  # Store all results
        self._loaded_rows = 0   # Track number of loaded rows

    def _init_ui(self):
        uic.loadUi(self.global_settings.get_ui_dir_path() + '/find_targets.ui', self)
        self.checkbox_select_all = self.findChild(QCheckBox, 'chkSelectAll')
        self.results_table = self.findChild(QTableWidget, 'tblTargets')
        
        # Connect select all checkbox signal
        self.checkbox_select_all.stateChanged.connect(self._on_select_all_changed)
        
        # Optimize table settings for large datasets
        self.results_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.results_table.setSelectionMode(QTableWidget.SelectionMode.MultiSelection)
        self.results_table.setShowGrid(False)
        self.results_table.setAlternatingRowColors(True)
        
        # Enable virtual scrolling mode
        self.results_table.setVerticalScrollMode(QTableWidget.ScrollMode.ScrollPerPixel)
        self.results_table.setHorizontalScrollMode(QTableWidget.ScrollMode.ScrollPerPixel)
        
        self.results_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.results_table.viewport().setProperty("cursor", Qt.CursorShape.ArrowCursor)
        
        self.results_table.setColumnCount(5) 
        headers = [
            "Feature Type", "Chromosome/Scaffold #", "Feature ID/Locus Tag", 
            "Feature Name", "Feature Description"
        ]
        self.results_table.setHorizontalHeaderLabels(headers)
        
        column_widths = [100, 150, 150, 150, 300]
        for i, width in enumerate(column_widths):
            self.results_table.setColumnWidth(i, width)
        
        self.results_table.horizontalHeader().setStretchLastSection(True)
        
        self.results_table.verticalScrollBar().valueChanged.connect(self._handle_scroll)

        self.push_button_generate_library = self.findChild(QPushButton, 'pbtnGenerateLibrary')
        self.push_button_view_targets = self.findChild(QPushButton, 'pbtnViewTargets')

        self.push_button_generate_library.clicked.connect(self._on_generate_library_clicked)

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
        """Display results with filtering support"""
        try:
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
            
        except Exception as e:
            self.logger.error(f"Error displaying results: {str(e)}")

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
        end_row = min(len(self._all_results), scroll_position + visible_rows * 2)
        
        # Load more rows if needed
        if end_row > self._loaded_rows:
            self._load_batch(self._loaded_rows, end_row)

    def get_selected_targets(self):
        """Get selected targets from the currently displayed (filtered) results"""
        try:
            # Get indices of selected rows in the current view
            selected_rows = set(index.row() for index in self.results_table.selectedIndexes())
            selected_targets = []
            
            if not selected_rows:
                self.logger.debug("No rows selected")
                return []
            
            # Get the currently visible rows from the table
            for row in selected_rows:
                try:
                    # Check if all required cells have valid data
                    cells = [self.results_table.item(row, col) for col in range(5)]
                    if any(cell is None for cell in cells):
                        self.logger.warning(f"Row {row} has missing data, skipping")
                        continue
                    
                    # Get data from visible row
                    target_data = {
                        'feature_type': cells[0].text(),
                        'chromosome': cells[1].text(),
                        'feature_id': cells[2].text(),
                        'feature_name': cells[3].text(),
                        'feature_description': cells[4].text()
                    }
                    
                    # Find corresponding full target data from _all_results
                    for full_target in self._all_results:
                        if (full_target['feature_id'] == target_data['feature_id'] and 
                            full_target['feature_type'] == target_data['feature_type']):
                            selected_targets.append(full_target)
                            break

                except Exception as row_error:
                    self.logger.warning(f"Error processing row {row}: {str(row_error)}")
                    continue

            self.logger.debug(f"Selected {len(selected_targets)} targets from filtered view")
            return selected_targets
            
        except Exception as e:
            self.logger.error(f"Error getting selected targets: {str(e)}")
            return []
    
    def clear_results(self):
        self.results_table.clearContents()
        self.results_table.setRowCount(0)   

    def _on_generate_library_clicked(self):
        """Handle generate library button click"""
        try:
            selected_targets = self.get_selected_targets()

            if not selected_targets:
                QtWidgets.QMessageBox.warning(
                    self,
                    "No Selection",
                    "Please select targets to generate library."
                )
                return
            
            # Store selected targets in global settings for persistence
            self.global_settings._current_selected_targets = selected_targets
            
            # Create and show generate library window
            from controllers.GenerateLibraryController import GenerateLibraryController
            generate_library_controller = GenerateLibraryController(
                self.global_settings,
                selected_targets
            )
            generate_library_controller.show()
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in generate library: {str(e)}")
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                f"An error occurred while opening the generate library window: {str(e)}"
            )

    def _on_select_all_changed(self, state):
        """Handle select all checkbox state changes"""
        try:
            self.results_table.setUpdatesEnabled(False)  # Disable updates for performance
            if state == Qt.CheckState.Checked.value:
                self.results_table.selectAll()
            else:
                self.results_table.clearSelection()
            self.results_table.setUpdatesEnabled(True)  # Re-enable updates
        except Exception as e:
            self.logger.error(f"Error in select all handler: {str(e)}")
