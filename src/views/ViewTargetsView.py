from typing import Optional
from PyQt6 import QtWidgets, uic, QtCore
from PyQt6.QtWidgets import QTableWidgetItem, QAbstractItemView
from PyQt6.QtGui import QTextDocument
from PyQt6.QtCore import Qt, pyqtSignal
from utils.ui import show_error
import time
import traceback

class ViewTargetsView(QtWidgets.QMainWindow):
    # Define the signal
    gene_selected = pyqtSignal(str)  # Signal to emit when gene is selected
    
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

        # Connect gene selection change with direct signal
        self.combo_box_gene.currentTextChanged.connect(self._on_gene_changed)

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
        """Ultra-fast target display using virtual table and minimal UI updates"""
        try:
            start_time = time.time()
            
            # Store complete set of targets if not already stored
            if not hasattr(self, '_complete_targets'):
                self._complete_targets = targets
            
            # Filter targets for currently selected gene
            selected_text = self.combo_box_gene.currentText()
            # Extract locus tag from "locus_tag: gene_name" format
            selected_locus = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
            
            if selected_locus:
                # Filter targets with more robust comparison
                filtered_targets = []
                for target in self._complete_targets:
                    target_locus = str(target.get('feature_id', '')).strip()
                    if target_locus.lower() == selected_locus.lower():
                        filtered_targets.append(target)
                
                # Store filtered results
                self._all_results = filtered_targets
            else:
                filtered_targets = self._complete_targets
                self._all_results = filtered_targets
                
            total_rows = len(filtered_targets)
            
            # Completely freeze UI
            self.setUpdatesEnabled(False)
            self.table_targets.setUpdatesEnabled(False)
            self.table_targets.setSortingEnabled(False)
            self.table_targets.setVisible(False)
            
            try:
                # Pre-allocate table
                self.table_targets.clearContents()
                self.table_targets.setRowCount(total_rows)
                
                # Get current headers to check for Azimuth column
                headers = self.get_table_headers()
                azimuth_index = headers.index("Azimuth 2.0") if "Azimuth 2.0" in headers else None
                
                # Pre-create flags once
                flags = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                
                # Load ALL rows at once
                for row in range(total_rows):
                    target = filtered_targets[row]
                    
                    # Create and set basic items
                    for col, value in enumerate([
                        target['location'], target['endonuclease'],
                        target['sequence'], target['strand'], target['pam']
                    ]):
                        item = QTableWidgetItem(str(value))
                        item.setFlags(flags)
                        self.table_targets.setItem(row, col, item)
                    
                    # Handle score separately for numeric sorting
                    score_item = QTableWidgetItem()
                    score_item.setData(QtCore.Qt.ItemDataRole.EditRole, float(target['score']))
                    self.table_targets.setItem(row, 5, score_item)
                    
                    # Add off-target placeholder
                    ot_item = QTableWidgetItem("--.--")
                    self.table_targets.setItem(row, 6, ot_item)
                    
                    # Create details button
                    details_button = QtWidgets.QPushButton("Details")
                    self.table_targets.setCellWidget(row, 7, details_button)
                    
                    # Add Azimuth score if column exists
                    if azimuth_index is not None and 'azimuth_score' in target:
                        azimuth_item = QTableWidgetItem()
                        azimuth_item.setData(QtCore.Qt.ItemDataRole.EditRole, float(target['azimuth_score']))
                        self.table_targets.setItem(row, azimuth_index, azimuth_item)
                
                # Set column widths
                column_widths = [100, 100, 200, 80, 80, 80, 80, 100]
                for col, width in enumerate(column_widths):
                    self.table_targets.setColumnWidth(col, width)
                
            finally:
                # Re-enable UI
                self.table_targets.setVisible(True)
                self.table_targets.setUpdatesEnabled(True)
                self.setUpdatesEnabled(True)
                self.table_targets.setSortingEnabled(True)
                
                total_time = time.time() - start_time
                self.logger.debug(f"Display time: {total_time:.2f} seconds for {total_rows} rows")
                
        except Exception as e:
            self.logger.error(f"Error in display_results: {str(e)}")
            show_error(self.settings, "Error displaying targets", str(e))

    def _handle_scroll_virtual(self, value, total_rows, row_height, buffer_rows):
        """Handle virtual scrolling with minimal updates"""
        try:
            if not hasattr(self, '_all_results') or not self._all_results:
                return
                
            # Calculate visible range with safety checks
            viewport_height = max(1, self.table_targets.viewport().height())
            row_height = max(1, row_height)  # Ensure non-zero
            visible_rows = viewport_height // row_height
            
            # Calculate which rows should be visible
            current_row = value // row_height if row_height > 0 else 0
            start_row = max(0, current_row - buffer_rows)
            end_row = min(total_rows, current_row + visible_rows + buffer_rows)
            
            # Only update rows that aren't already loaded
            for row in range(start_row, end_row):
                if row < len(self._all_results) and not self.table_targets.item(row, 0):
                    target = self._all_results[row]
                    
                    # Create and set items efficiently
                    for col, value in enumerate([
                        target['location'], target['endonuclease'],
                        target['sequence'], target['strand'], target['pam'],
                        target['score'], "--.--"
                    ]):
                        item = QTableWidgetItem(str(value))
                        item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                        self.table_targets.setItem(row, col, item)
                    
                    if not self.table_targets.cellWidget(row, 7):
                        details_button = QtWidgets.QPushButton("Details")
                        self.table_targets.setCellWidget(row, 7, details_button)
                        
        except Exception as e:
            self.logger.error(f"Error in _handle_scroll_virtual: {str(e)}")

    def get_selected_targets(self):
        """Get selected targets with all necessary data"""
        try:
            selected_rows = set(index.row() for index in self.table_targets.selectedIndexes())
            selected_targets = []
            
            # Get column indices once
            columns = {
                'location': 0,
                'endonuclease': 1,
                'sequence': 2,
                'strand': 3,
                'pam': 4,
                'score': 5
            }
            
            for row in sorted(selected_rows):
                # Verify all required cells have data
                if all(self.table_targets.item(row, col) is not None 
                      for col in columns.values()):
                    
                    target = {
                        'location': self.table_targets.item(row, columns['location']).text(),
                        'endonuclease': self.table_targets.item(row, columns['endonuclease']).text(),
                        'sequence': self.table_targets.item(row, columns['sequence']).text(),
                        'strand': self.table_targets.item(row, columns['strand']).text(),
                        'pam': self.table_targets.item(row, columns['pam']).text(),
                        'score': self.table_targets.item(row, columns['score']).text()
                    }
                    selected_targets.append(target)
                else:
                    self.logger.warning(f"Skipping row {row} due to missing data")
                    
            if not selected_targets:
                self.logger.warning("No valid targets selected")
                
            return selected_targets
            
        except Exception as e:
            self.logger.error(f"Error getting selected targets: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return []

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
        """Set genes in combo box with optimized performance"""
        try:
            start_time = time.time()
            
            # Disable UI updates
            self.combo_box_gene.blockSignals(True)
            self.combo_box_gene.setUpdatesEnabled(False)
            
            # Clear existing items efficiently
            self.combo_box_gene.clear()
            
            # Debug logging
            self.logger.debug(f"Received {len(genes)} genes")
            
            # Add items in a single batch
            if genes:
                # Pre-allocate size
                self.combo_box_gene.insertItems(0, genes)
                
                # Set first item without triggering updates
                if self.combo_box_gene.count() > 0:
                    self.combo_box_gene.setCurrentIndex(0)
                    
                self.logger.debug(f"Added {len(genes)} genes to combo box")
                
            # Re-enable UI updates
            self.combo_box_gene.setUpdatesEnabled(True)
            self.combo_box_gene.blockSignals(False)
            
            total_time = time.time() - start_time
            self.logger.debug(f"Combo box update time: {total_time:.2f} seconds")
            
        except Exception as e:
            self.logger.error(f"Error setting genes in combo box: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def set_text_edit_gene_viewer(self, sequence):
        """Update gene viewer with new sequence"""
        try:
            if sequence:
                self.text_edit_gene_viewer.setText(sequence)
                self.logger.debug(f"Updated gene viewer with sequence of length: {len(sequence)}")
            else:
                self.text_edit_gene_viewer.clear()
                self.logger.debug("Cleared gene viewer - no sequence provided")
        except Exception as e:
            self.logger.error(f"Error setting gene viewer text: {str(e)}")

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

    def _on_gene_changed(self, selected_text):
        """Handle gene selection change and emit signal"""
        try:
            self.logger.debug(f"Gene selection changed to: {selected_text}")
            
            # Reset scroll position
            self.table_targets.verticalScrollBar().setValue(0)
            
            # Filter and display targets
            if hasattr(self, '_complete_targets'):
                self.display_targets_in_table(self._complete_targets)
            
            # Emit signal for controller to update gene sequence
            self.gene_selected.emit(selected_text)
            
        except Exception as e:
            self.logger.error(f"Error in _on_gene_changed: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def get_table_headers(self):
        """Get current table headers"""
        headers = []
        for i in range(self.table_targets.columnCount()):
            headers.append(self.table_targets.horizontalHeaderItem(i).text())
        return headers

    def add_scoring_column(self, algorithm_name, position=None):
        """Add a new column for alternative scoring method at specified position"""
        if position is None:
            # Add to end if no position specified
            position = self.table_targets.columnCount()
        
        self.table_targets.insertColumn(position)
        self.table_targets.setHorizontalHeaderItem(
            position, 
            QtWidgets.QTableWidgetItem(algorithm_name)
        )
        
        # Shift any existing columns after the insertion point
        for i in range(self.table_targets.columnCount() - 1, position, -1):
            for row in range(self.table_targets.rowCount()):
                self.table_targets.setItem(row, i, self.table_targets.takeItem(row, i-1))
                
            # Move column header
            header_item = self.table_targets.takeHorizontalHeaderItem(i-1)
            if header_item:
                self.table_targets.setHorizontalHeaderItem(i, header_item)
        
        return position
