from typing import Optional
from PyQt6 import QtWidgets, uic
from PyQt6.QtWidgets import QTableWidgetItem, QAbstractItemView
from PyQt6.QtGui import QTextDocument
from PyQt6.QtCore import Qt, pyqtSignal
from utils.ui import show_error
import traceback
from views.DNAFeatureViewer import DNAFeatureViewer

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

        self.push_button_export_selected_grnas = self._find_widget('pbtnExportSelectedgRNAs', QtWidgets.QPushButton)

    def _init_grpGuideViewer(self):
        self.combo_box_gene = self._find_widget('cmbGene', QtWidgets.QComboBox)
        self.combo_box_endonuclease = self._find_widget('cmbEndonuclease', QtWidgets.QComboBox)
        self.check_box_filter_5_prime_g_sequences = self._find_widget('chkFilter5PrimeG', QtWidgets.QCheckBox)
        self.spin_box_minimum_on_target_score = self._find_widget('spnMinOTScore', QtWidgets.QSpinBox)
        self.check_box_select_all = self._find_widget('chkSelectAll', QtWidgets.QCheckBox)
        self.push_button_scoring_options = self._find_widget('pbtnScoringOptions', QtWidgets.QPushButton)
        self.table_guides = self._find_widget('tblGuides', QtWidgets.QTableWidget)

        self.combo_box_gene.currentTextChanged.connect(self._on_gene_changed)

        self.table_guides.setColumnCount(8)
        self.table_guides.setHorizontalHeaderLabels(["Location", "Endonuclease", "Sequence", "Strand", "PAM", "Score", "Off-Target", "Details"])
        self.table_guides.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_guides.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_guides.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        
        # Enable sorting
        self.table_guides.setSortingEnabled(True)
        
        # Enable horizontal scrolling
        self.table_guides.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table_guides.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Set size policy to allow table to shrink and expand
        self.table_guides.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding
        )
        
        # Set resize mode for header
        header = self.table_guides.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(False) 
        
        # Set minimum section size to prevent columns from becoming too narrow
        header.setMinimumSectionSize(80)

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
        self.check_box_view_exons_only = self._find_widget('chkViewExonsOnly', QtWidgets.QCheckBox)

        self.text_edit_gene_viewer.setReadOnly(True)

        # Create DNA feature viewer
        self.dna_feature_viewer = DNAFeatureViewer()
        
        # Get the layout of the gene viewer group
        gene_viewer_group = self.findChild(QtWidgets.QGroupBox, 'grpGeneViewer')
        gene_viewer_layout = gene_viewer_group.layout()
        
        # Find the row index of the text editor
        text_editor_row = -1
        for i in range(gene_viewer_layout.rowCount()):
            item = gene_viewer_layout.itemAtPosition(i, 0)
            if item and item.widget() == self.text_edit_gene_viewer:
                text_editor_row = i
                break
        
        if text_editor_row != -1:
            # Insert DNA feature viewer above the text editor
            gene_viewer_layout.addWidget(self.dna_feature_viewer, text_editor_row, 0, 1, -1)
        
        # Connect signals
        self.dna_feature_viewer.sequence_selected.connect(self._on_sequence_selected)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget 

    def display_guides_in_table(self, guides):
        try:
            self._all_guides = guides
            
            selected_text = self.combo_box_gene.currentText()
            
            # First filter by position/feature
            if selected_text and "chrom" in selected_text and "start:" in selected_text:
                filtered_guides = []
                for guide in self._all_guides:
                    if guide.get('feature_id') == selected_text:
                        filtered_guides.append(guide)
                self.logger.debug(f"Filtered to {len(filtered_guides)} guides for position {selected_text}")
            else:
                selected_locus = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
                
                if selected_locus:
                    filtered_guides = []
                    for guide in self._all_guides:
                        guide_locus = str(guide.get('feature_id', '')).strip()
                        if guide_locus.lower() == selected_locus.lower():
                            filtered_guides.append(guide)
                else:
                    filtered_guides = self._all_guides
            
            # Apply additional filters
            final_guides = []
            for guide in filtered_guides:
                # Filter by minimum score
                min_score = self.spin_box_minimum_on_target_score.value()
                if float(guide.get('score', 0)) < min_score:
                    continue
                    
                # Filter by 5' G sequences
                if self.check_box_filter_5_prime_g_sequences.isChecked():
                    sequence = guide.get('sequence', '')
                    if not sequence or not sequence.startswith('G'):
                        continue
                        
                final_guides.append(guide)
                
            # Update table with new guides
            total_rows = len(final_guides)
            self.logger.debug(f"Processing {total_rows} rows for display after filtering")
            
            # Completely freeze UI
            self.setUpdatesEnabled(False)
            self.table_guides.setUpdatesEnabled(False)
            self.table_guides.setSortingEnabled(False)
            self.table_guides.setVisible(False)
            
            try:
                # Clear and resize table
                self.table_guides.clearContents()
                self.table_guides.setRowCount(total_rows)
                
                # Get current headers to check for Azimuth column
                headers = self.get_table_headers()
                azimuth_index = headers.index("Azimuth 2.0") if "Azimuth 2.0" in headers else None
                
                # Pre-create flags once
                flags = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                
                # Load ALL rows at once
                for row in range(total_rows):
                    guide = final_guides[row]
                    
                    # Extract start position from location (format: "start-end")
                    location = guide['location']
                    start_pos = location.split('-')[0] if '-' in location else location
                    
                    # Create items with proper data roles for sorting
                    items = [
                        (0, self._create_sortable_item(start_pos, int(start_pos))),  # Location as number
                        (1, QTableWidgetItem(guide['endonuclease'])),
                        (2, QTableWidgetItem(guide['sequence'])),
                        (3, QTableWidgetItem(guide['strand'])),
                        (4, QTableWidgetItem(guide['pam'])),
                        (5, self._create_sortable_item(str(guide['score']), float(guide['score']))),  # Score as number
                        (6, QTableWidgetItem("--.--"))  # Off-target placeholder
                    ]
                    
                    for col, item in items:
                        item.setFlags(flags)
                        self.table_guides.setItem(row, col, item)
                    
                    # Only add details button if sequence has off-target details
                    sequence = guide['sequence']
                    if hasattr(self, '_off_target_details') and sequence in self._off_target_details:
                        details_button = QtWidgets.QPushButton("Details")
                        details_button.clicked.connect(self._show_details)
                        self.table_guides.setCellWidget(row, 7, details_button)
                    
                    # Add Azimuth score if column exists
                    if azimuth_index is not None and 'azimuth_score' in guide:
                        azimuth_score = float(guide['azimuth_score'])
                        azimuth_item = self._create_sortable_item(str(azimuth_score), azimuth_score)
                        self.table_guides.setItem(row, azimuth_index, azimuth_item)
                
                # Updated column widths
                column_widths = [
                    80,  # Location
                    100,  # Endonuclease
                    200,  # Sequence
                    10,   # Strand
                    80,  # PAM
                    10,   # Score
                    30,   # Off-Target
                    80   # Details
                ]
                
                # Set the column widths
                for col, width in enumerate(column_widths):
                    self.table_guides.setColumnWidth(col, width)
                
                essential_columns_width = sum(column_widths[:8])  # First 6 columns
                self.table_guides.setMinimumWidth(essential_columns_width)
                
                # Update the group box to properly handle scrolling
                guide_viewer_group = self.findChild(QtWidgets.QGroupBox, 'grpGuideViewer')
                guide_viewer_group.setMinimumWidth(essential_columns_width + 50)  # Add some padding for scrollbar
                
            finally:
                # Re-enable UI
                self.table_guides.setVisible(True)
                self.table_guides.setUpdatesEnabled(True)
                self.setUpdatesEnabled(True)
                self.table_guides.setSortingEnabled(True)
                
        except Exception as e:
            self.logger.error(f"Error in display_guides: {str(e)}")
            show_error(self.settings, "Error displaying guides", str(e))

    def _create_sortable_item(self, display_text, sort_value):
        """Create a table item that displays text but sorts by numeric value"""
        item = QTableWidgetItem()
        item.setData(Qt.ItemDataRole.DisplayRole, sort_value)  # Use raw value for display
        item.setData(Qt.ItemDataRole.EditRole, sort_value)    # Used for sorting
        
        # Format display text based on value type
        if isinstance(sort_value, (int, float)):
            if isinstance(sort_value, int):
                # For integers (like positions), show full number
                item.setText(f"{sort_value:d}")
            else:
                # For floats (like scores), show with 2 decimal places
                item.setText(f"{sort_value:.2f}")
        else:
            item.setText(str(sort_value))
        
        return item

    def _handle_scroll_virtual(self, value, total_rows, row_height, buffer_rows):
        try:
            if not hasattr(self, '_all_guides') or not self._all_guides:
                return
                
            # Calculate visible range with safety checks
            viewport_height = max(1, self.table_guides.viewport().height())
            row_height = max(1, row_height)  # Ensure non-zero
            visible_rows = viewport_height // row_height
            
            # Calculate which rows should be visible
            current_row = value // row_height if row_height > 0 else 0
            start_row = max(0, current_row - buffer_rows)
            end_row = min(total_rows, current_row + visible_rows + buffer_rows)
            
            # Only update rows that aren't already loaded
            for row in range(start_row, end_row):
                if row < len(self._all_guides) and not self.table_guides.item(row, 0):
                    guide = self._all_guides[row]
                    
                    # Extract start position from location
                    location = guide['location']
                    start_pos = location.split('-')[0] if '-' in location else location
                    
                    # Create items with proper data roles for sorting
                    items = [
                        (0, self._create_sortable_item(start_pos, int(start_pos))),  # Location as number
                        (1, QTableWidgetItem(guide['endonuclease'])),
                        (2, QTableWidgetItem(guide['sequence'])),
                        (3, QTableWidgetItem(guide['strand'])),
                        (4, QTableWidgetItem(guide['pam'])),
                        (5, self._create_sortable_item(str(guide['score']), float(guide['score']))),  # Score as number
                        (6, QTableWidgetItem("--.--"))  # Off-target placeholder
                    ]
                    
                    # Set items with flags
                    for col, item in items:
                        item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                        self.table_guides.setItem(row, col, item)
                    
                    # Add details button
                    if not self.table_guides.cellWidget(row, 7):
                        details_button = QtWidgets.QPushButton("Details")
                        self.table_guides.setCellWidget(row, 7, details_button)
                        
        except Exception as e:
            self.logger.error(f"Error in _handle_scroll_virtual: {str(e)}")

    def get_selected_guides(self):
        """Get selected guides with all necessary data"""
        try:
            selected_rows = sorted(set(item.row() for item in self.table_guides.selectedItems()))
            selected_guides = []
            
            # Get column indices once
            columns = {
                'location': 0,
                'endonuclease': 1,
                'sequence': 2,
                'strand': 3,
                'pam': 4,
                'score': 5,
                'off_target': 6
            }
            
            # Get current gene information from combo box
            current_gene = self.combo_box_gene.currentText()
            if ': ' in current_gene:  # Format is "locus_tag: gene_name"
                locus_tag, gene_name = current_gene.split(': ', 1)
            else:
                locus_tag = current_gene
                gene_name = current_gene
            
            for row in selected_rows:
                # Create guide dictionary directly from table items
                guide = {}
                valid_row = True
                
                for col_name, col_index in columns.items():
                    item = self.table_guides.item(row, col_index)
                    if item is None:
                        valid_row = False
                        self.logger.warning(f"Missing data in row {row}, column {col_name}")
                        break
                    guide[col_name] = item.text()
                
                if valid_row:
                    # Add gene information
                    guide['locus_tag'] = locus_tag.strip()
                    guide['gene_name'] = gene_name.strip()
                    selected_guides.append(guide)
                    
            if not selected_guides:
                self.logger.warning("No valid guides selected")
                
            return selected_guides
            
        except Exception as e:
            self.logger.error(f"Error getting selected guides: {str(e)}")
            return []

    def get_row_data(self, row):
        return {
            'location': self.table_guides.item(row, 0).text(),
            'endonuclease': self.table_guides.item(row, 1).text(),
            'sequence': self.table_guides.item(row, 2).text(),
            'strand': self.table_guides.item(row, 3).text(),
            'pam': self.table_guides.item(row, 4).text(),
            'score': self.table_guides.item(row, 5).text(),
            'off_target': self.table_guides.item(row, 6).text()
        }
    
    def set_combo_box_endonuclease(self, endonucleases):
        self.combo_box_endonuclease.addItems(endonucleases)

    def set_combo_box_gene(self, genes):
        try:
            # Disable UI updates
            self.combo_box_gene.blockSignals(True)
            self.combo_box_gene.setUpdatesEnabled(False)
            
            # Clear existing items efficiently
            self.combo_box_gene.clear()
            
            # Debug logging
            self.logger.debug(f"Received {len(genes)} genes")
            
            # Use a set to ensure uniqueness
            unique_genes = list(set(genes))
            
            # Add items in a single batch
            if unique_genes:
                # Pre-allocate size
                self.combo_box_gene.insertItems(0, unique_genes)
                
                # Set first item without triggering updates
                if self.combo_box_gene.count() > 0:
                    self.combo_box_gene.setCurrentIndex(0)
                    
                self.logger.debug(f"Added {len(unique_genes)} unique genes to combo box")
                
            # Re-enable UI updates
            self.combo_box_gene.setUpdatesEnabled(True)
            self.combo_box_gene.blockSignals(False)
            
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

    def update_gene_viewer(self, sequence, features=None):
        """Update both text editor and DNA feature viewer"""
        # Update text editor
        self.text_edit_gene_viewer.clear()
        doc = QTextDocument()
        doc.setHtml(sequence)
        self.text_edit_gene_viewer.setDocument(doc)
        
        # Get start position from line edit
        try:
            start_pos = int(self.line_edit_start_location.text())
        except (ValueError, TypeError):
            start_pos = 1
        
        # Update DNA feature viewer
        if features is None:
            features = []
        self.dna_feature_viewer.set_data(sequence, features, start_pos)

    def select_all_guides(self, select):
        for row in range(self.table_guides.rowCount()):
            self.table_guides.selectRow(row) if select else self.table_guides.clearSelection()

    def get_export_file_path(self):
        # Implement this method to get the export file path from the user
        return QtWidgets.QFileDialog.getSaveFileName(self, 'Save File')[0]

    def _on_gene_changed(self, selected_text):
        """Handle gene selection change and emit signal"""
        try:
            self.logger.debug(f"Gene selection changed to: {selected_text}")
            
            # Reset scroll position
            self.table_guides.verticalScrollBar().setValue(0)
            
            # Filter and display targets
            if hasattr(self, '_complete_targets'):
                self.display_guides_in_table(self._complete_targets)
            
            # Emit signal for controller to update gene sequence
            self.gene_selected.emit(selected_text)
            
        except Exception as e:
            self.logger.error(f"Error in _on_gene_changed: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def get_table_headers(self):
        """Get current table headers"""
        headers = []
        for i in range(self.table_guides.columnCount()):
            headers.append(self.table_guides.horizontalHeaderItem(i).text())
        return headers

    def add_scoring_column(self, algorithm_name, position=None):
        """Add a new column for alternative scoring method at specified position"""
        if position is None:
            # Add to end if no position specified
            position = self.table_guides.columnCount()
        
        self.table_guides.insertColumn(position)
        self.table_guides.setHorizontalHeaderItem(
            position, 
            QtWidgets.QTableWidgetItem(algorithm_name)
        )
        
        # Shift any existing columns after the insertion point
        for i in range(self.table_guides.columnCount() - 1, position, -1):
            for row in range(self.table_guides.rowCount()):
                self.table_guides.setItem(row, i, self.table_guides.takeItem(row, i-1))
                
            # Move column header
            header_item = self.table_guides.takeHorizontalHeaderItem(i-1)
            if header_item:
                self.table_guides.setHorizontalHeaderItem(i, header_item)
        
        return position

    def update_off_target_details(self, off_target_results, detailed_results=None):
        """Update off-target scores and details"""
        try:
            # Store detailed results if provided
            if detailed_results:
                self._off_target_details = detailed_results
            
            # Update off-target scores in table
            for row in range(self.table_guides.rowCount()):
                sequence = self.table_guides.item(row, 2).text()
                if sequence in off_target_results:
                    score = off_target_results[sequence]
                    score_item = QTableWidgetItem(str(score))
                    self.table_guides.setItem(row, 6, score_item)
                    
                    # Add details button if detailed results exist
                    if detailed_results and sequence in detailed_results:
                        details_button = QtWidgets.QPushButton("Details")
                        details_button.clicked.connect(self._show_details)
                        self.table_guides.setCellWidget(row, 7, details_button)
                    
            self.table_guides.resizeColumnsToContents()
            
        except Exception as e:
            self.logger.error(f"Error updating off-target details: {str(e)}")
            show_error(self.settings, "Error updating off-target details", str(e))

    def _show_details(self):
        """Show off-target details dialog"""
        try:
            button = self.sender()
            index = self.table_guides.indexAt(button.pos())
            sequence = self.table_guides.item(index.row(), 2).text()
            
            if sequence in self._off_target_details:
                details = self._off_target_details[sequence]
                
                msg = QtWidgets.QMessageBox()
                msg.setWindowTitle("Details")
                
                # Format details message
                chromo_str = "<html><b>Reference gRNA:</b><br>Location, Sequence, Strand, PAM, On Score<br></html>"
                input_str = (f"{self.table_guides.item(index.row(),0).text()}, {sequence}, "
                            f"{self.table_guides.item(index.row(),3).text()}, "
                            f"{self.table_guides.item(index.row(),4).text()}, "
                            f"{self.table_guides.item(index.row(),5).text()}<br><br>")
                detail_str = "<html><b>Off-Target Hits:</b><br>Off Score, Chromosome, Location, Sequence<br></html>"
                
                msg.setText(chromo_str + input_str + detail_str + "<br>".join(details))
                msg.exec()
                
        except Exception as e:
            self.logger.error(f"Error showing details: {str(e)}")
            show_error(self.settings, "Error showing details", str(e))

    def _on_sequence_selected(self, start, end):
        """Handle sequence selection in DNA feature viewer"""
        self.line_edit_start_location.setText(str(start))
        self.line_edit_stop_location.setText(str(end))

    def highlight_guides_in_viewer(self, guides_to_highlight, sequence):
        """Highlight guides in viewer"""
        try:
            for guide in guides_to_highlight:
                sequence_to_find = guide['sequence']
                strand = guide['strand']
                
                if strand == '-':
                    sequence_to_find = str(Seq(sequence_to_find).reverse_complement())
                
                sequence_upper = sequence.upper()
                target_upper = sequence_to_find.upper()
                
                pos = sequence_upper.find(target_upper)
                if pos != -1:
                    # Set color based on strand
                    color = QColor(255, 0, 0, 100) if strand == '-' else QColor(0, 255, 0, 100)
                    
                    # Highlight sequence in viewer
                    self.dna_feature_viewer.sequence_viewer.highlight_sequence(
                        pos, 
                        pos + len(sequence_to_find) - 1,
                        color
                    )
                    
        except Exception as e:
            self.logger.error(f"Error highlighting guides: {str(e)}")
            show_error(self.settings, "Error highlighting guides", str(e))

    def clear_highlights(self):
        """Clear highlights in viewer"""
        self.dna_feature_viewer.sequence_viewer.clear_highlights()