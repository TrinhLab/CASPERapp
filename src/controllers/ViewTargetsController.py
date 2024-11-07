import logging
from controllers.ScoringOptionsController import ScoringOptionsController
from models.ViewTargetsModel import ViewTargetsModel
from views.ViewTargetsView import ViewTargetsView
from PyQt6.QtWidgets import QMessageBox
from utils.ui import show_error
import time
from PyQt6 import QtWidgets, QtCore
import traceback
import threading
from Bio.Seq import Seq

class ViewTargetsController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        start_time = time.time()
        self.model = ViewTargetsModel(global_settings)
        self.view = ViewTargetsView(global_settings)
        init_time = time.time() - start_time
        self.global_settings.logger.debug(f"ViewTargets initialization took: {init_time:.2f} seconds")
        
        self.setup_connections()
        self.organism = ""
        self.endonuclease = ""

    def setup_connections(self):
        self.view.push_button_off_target.clicked.connect(self.perform_off_target_analysis)
        self.view.push_button_cotargeting.clicked.connect(self.perform_cotargeting)
        self.view.push_button_highlight_guides.clicked.connect(self.highlight_gene_viewer)
        self.view.push_button_export_grna.clicked.connect(self.export_targets)
        self.view.push_button_filter_options.clicked.connect(self.show_filter_options)
        self.view.push_button_scoring_options.clicked.connect(self.show_scoring_options)
        self.view.push_button_change_location.clicked.connect(self.change_indices)
        self.view.push_button_reset_location.clicked.connect(self.reset_location)
        self.view.check_box_select_all.stateChanged.connect(self.select_all)
        self.view.combo_box_gene.currentIndexChanged.connect(self.display_gene_data)
        self.view.gene_selected.connect(self.on_gene_selected)

    def load_targets(self, selected_targets, organism, endonuclease):
        try:
            total_start = time.time()
            
            self.organism = organism
            self.endonuclease = endonuclease
            
            # Time model loading
            model_start = time.time()
            self.model.load_targets(selected_targets, organism, endonuclease)
            model_time = time.time() - model_start
            self.global_settings.logger.debug(f"Model load_targets took: {model_time:.2f} seconds")
            
            # Time getting targets
            targets_start = time.time()
            targets = self.model.get_targets()
            targets_time = time.time() - targets_start
            self.global_settings.logger.debug(f"Getting targets took: {targets_time:.2f} seconds")
            
            # Get feature ID mapping from FindTargetsModel - Optimized with timing
            genes_start = time.time()
            
            # Time set creation
            set_start = time.time()
            seen_genes = set()
            formatted_genes = []
            set_time = time.time() - set_start
            self.global_settings.logger.debug(f"Set initialization took: {set_time:.2f} seconds")
            
            # Time target processing
            process_start = time.time()
            for target in selected_targets:
                gene_name = target.get('feature_name')
                feature_id = target.get('feature_id')
                
                if gene_name and feature_id and gene_name not in seen_genes:
                    seen_genes.add(gene_name)
                    formatted_genes.append(f"{feature_id}: {gene_name}")
            process_time = time.time() - process_start
            self.global_settings.logger.debug(f"Target processing took: {process_time:.2f} seconds")
            
            # Time sorting
            sort_start = time.time()
            formatted_genes.sort()
            sort_time = time.time() - sort_start
            self.global_settings.logger.debug(f"Sorting took: {sort_time:.2f} seconds")
            
            # Time view update
            view_start = time.time()
            self.view.set_combo_box_gene(formatted_genes)
            view_time = time.time() - view_start
            self.global_settings.logger.debug(f"View update took: {view_time:.2f} seconds")
            
            genes_time = time.time() - genes_start
            self.global_settings.logger.debug(f"Total setting genes took: {genes_time:.2f} seconds")
            
            # Time displaying targets
            display_start = time.time()
            self.view.display_targets_in_table(targets)
            display_time = time.time() - display_start
            self.global_settings.logger.debug(f"Displaying targets took: {display_time:.2f} seconds")
            
            # Time setting endonuclease
            endo_start = time.time()
            self.view.set_combo_box_endonuclease([endonuclease])
            endo_time = time.time() - endo_start
            self.global_settings.logger.debug(f"Setting endonuclease took: {endo_time:.2f} seconds")
            
            # Time loading gene viewer
            gene_start = time.time()
            self.load_gene_viewer()
            gene_time = time.time() - gene_start
            self.global_settings.logger.debug(f"Loading gene viewer took: {gene_time:.2f} seconds")
            
            total_time = time.time() - total_start
            self.global_settings.logger.debug(f"Total load_targets took: {total_time:.2f} seconds")
            
        except Exception as e:
            show_error(self.global_settings, "Error loading targets", str(e))

    def load_gene_viewer(self):
        """Load gene viewer with sequence and location information"""
        try:
            total_start = time.time()
            
            # Get selected gene from combo box
            combo_start = time.time()
            selected_text = self.view.combo_box_gene.currentText()
            if not selected_text:
                self.global_settings.logger.debug("No gene selected")
                return
            combo_time = time.time() - combo_start
            self.global_settings.logger.debug(f"Combo box access time: {combo_time:.2f} seconds")
            
            # Extract locus tag from "locus_tag: gene_name" format
            parse_start = time.time()
            locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
            self.global_settings.logger.debug(f"Loading sequence for locus tag: {locus_tag}")
            parse_time = time.time() - parse_start
            self.global_settings.logger.debug(f"Locus tag parsing time: {parse_time:.2f} seconds")
            
            # Get gene sequence with padding
            sequence_start = time.time()
            sequence_data = self.model.get_gene_sequence(locus_tag)
            sequence_time = time.time() - sequence_start
            self.global_settings.logger.debug(f"Sequence retrieval time: {sequence_time:.2f} seconds")
            
            if sequence_data:
                # Update gene viewer with sequence
                viewer_start = time.time()
                self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
                viewer_time = time.time() - viewer_start
                self.global_settings.logger.debug(f"Text viewer update time: {viewer_time:.2f} seconds")
                
                # Update location fields
                location_start = time.time()
                self.view.line_edit_start_location.setText(str(sequence_data['start']))
                self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                location_time = time.time() - location_start
                self.global_settings.logger.debug(f"Location fields update time: {location_time:.2f} seconds")
                
                total_time = time.time() - total_start
                self.global_settings.logger.debug(f"Total gene viewer loading took: {total_time:.2f} seconds")
            else:
                self.global_settings.logger.warning(f"No sequence data found for locus tag {locus_tag}")
                
        except Exception as e:
            self.global_settings.logger.error(f"Error in load_gene_viewer: {str(e)}")
            self.global_settings.logger.error(f"Stack trace: {traceback.format_exc()}")

    def perform_off_target_analysis(self):
        try:
            selected_targets = self.view.get_selected_targets()
            if not selected_targets:
                QMessageBox.warning(self.view, "No Selection", "Please select targets for off-target analysis.")
                return
            # Implement off-target analysis logic here
            # You might want to create a new controller for off-target analysis
            off_target_controller = self.global_settings.get_off_target_window()
            off_target_controller.analyze(selected_targets, self.organism, self.endonuclease)
        except Exception as e:
            show_error(self.global_settings, "Error in off-target analysis", str(e))

    def perform_cotargeting(self):
        try:
            selected_targets = self.view.get_selected_targets()
            if not selected_targets:
                QMessageBox.warning(self.view, "No Selection", "Please select targets for cotargeting.")
                return
            # Implement cotargeting logic here
            # You might want to create a new controller for cotargeting
            cotargeting_controller = self.global_settings.get_cotargeting_window()
            cotargeting_controller.analyze(selected_targets, self.organism, self.endonuclease)
        except Exception as e:
            show_error(self.global_settings, "Error in cotargeting", str(e))

    def highlight_gene_viewer(self):
        """Highlight selected targets in gene viewer"""
        try:
            self.global_settings.logger.debug("Starting highlight_gene_viewer")
            
            # Get selected targets
            selected_rows = self.view.get_selected_targets()
            self.global_settings.logger.debug(f"Selected targets: {selected_rows}")
            
            if not selected_rows:
                QMessageBox.warning(self.view, "No Selection", 
                                  "Please select targets to highlight in the gene viewer.")
                return

            # Convert table selections to the format expected by the model
            targets_to_highlight = []
            for target in selected_rows:
                target_info = {
                    'location': target['location'],
                    'sequence': target['sequence'],
                    'strand': target['strand']
                }
                targets_to_highlight.append(target_info)
                self.global_settings.logger.debug(f"Target to highlight: {target_info}")

            # Get current gene sequence
            current_gene = self.view.combo_box_gene.currentText()
            locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
            self.global_settings.logger.debug(f"Getting sequence for locus tag: {locus_tag}")
            
            # Get gene sequence with padding
            sequence_data = self.model.get_gene_sequence(locus_tag)
            if not sequence_data or 'sequence' not in sequence_data:
                self.global_settings.logger.error("No sequence data found")
                QMessageBox.warning(self.view, "No Gene Data", 
                                  "Could not get gene sequence for highlighting.")
                return

            self.global_settings.logger.debug(f"Gene sequence length: {len(sequence_data['sequence'])}")
            
            # Highlight the sequences
            if targets_to_highlight:
                self.global_settings.logger.debug("Attempting to highlight sequences")
                self.highlight_targets_in_gene_viewer(targets_to_highlight)
            else:
                self.global_settings.logger.error("No valid targets to highlight")
                QMessageBox.warning(self.view, "No Valid Targets", 
                                  "Could not get sequence information from the selected rows.")

        except Exception as e:
            self.global_settings.logger.error(f"Error in highlight_gene_viewer: {str(e)}\n{traceback.format_exc()}")
            show_error(self.global_settings, "Error highlighting gene viewer", str(e))

    def export_targets(self):
        try:
            selected_targets = self.view.get_selected_targets()
            if not selected_targets:
                QMessageBox.warning(self.view, "No Selection", "Please select targets to export.")
                return
            file_path = self.view.get_export_file_path()
            if file_path:
                self.model.export_targets(selected_targets, file_path)
                QMessageBox.information(self.view, "Export Successful", "Selected targets have been exported successfully.")
        except Exception as e:
            show_error(self.global_settings, "Error exporting targets", str(e))

    def show_filter_options(self):
        try:
            filter_options = self.model.get_filter_options()
            self.view.show_filter_options_dialog(filter_options)
            if self.view.filter_options_accepted():
                new_options = self.view.get_filter_options()
                self.model.set_filter_options(new_options)
                self.refresh_targets_display()
        except Exception as e:
            show_error(self.global_settings, "Error showing filter options", str(e))

    def show_scoring_options(self):
        """Show scoring options window"""
        try:
            # Create scoring options controller if not exists
            if not hasattr(self, '_scoring_options_controller'):
                # Create controller with self as view_targets_controller
                self._scoring_options_controller = ScoringOptionsController(
                    global_settings=self.global_settings,
                    view_targets_controller=self
                )
                
            # Show scoring options window
            self._scoring_options_controller.show()
            
        except Exception as e:
            self.global_settings.logger.error(f"Error showing scoring options: {str(e)}")
            self.global_settings.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.global_settings, "Error", f"Could not show scoring options: {str(e)}")

    def change_indices(self):
        try:
            start = int(self.view.line_edit_start_location.text())
            end = int(self.view.line_edit_stop_location.text())
            if self.model.update_gene_viewer_indices(start, end):
                self.view.set_text_edit_gene_viewer(self.model.gene_sequence)
            else:
                QMessageBox.warning(self.view, "Invalid Range", "Please enter valid start and end positions within the gene range.")
        except ValueError:
            QMessageBox.warning(self.view, "Invalid Input", "Please enter valid integer values for start and end positions.")
        except Exception as e:
            show_error(self.global_settings, "Error changing indices", str(e))

    def reset_location(self):
        try:
            self.model.reset_gene_viewer_indices()
            gene_data = self.model.get_gene_data(self.view.combo_box_gene.currentText())
            self.view.set_text_edit_gene_viewer(gene_data['sequence'])
            self.view.line_edit_start_location.setText(str(gene_data['start']))
            self.view.line_edit_stop_location.setText(str(gene_data['end']))
        except Exception as e:
            show_error(self.global_settings, "Error resetting location", str(e))

    def select_all(self, state):
        try:
            self.view.select_all_targets(state == 2)  # 2 corresponds to Qt.Checked
        except Exception as e:
            show_error(self.global_settings, "Error selecting all targets", str(e))

    def display_gene_data(self, gene_name):
        try:
            gene_data = self.model.get_gene_data(gene_name)
            if gene_data and gene_data['sequence']:
                self.view.set_text_edit_gene_viewer(gene_data['sequence'])
                self.view.line_edit_start_location.setText(str(gene_data['start']))
                self.view.line_edit_stop_location.setText(str(gene_data['end']))
            else:
                self.view.set_text_edit_gene_viewer("No sequence data available for this gene")
                self.view.line_edit_start_location.clear()
                self.view.line_edit_stop_location.clear()
        except Exception as e:
            show_error(self.global_settings, "Error displaying gene data", str(e))

    def refresh_targets_display(self):
        try:
            filtered_targets = self.model.get_filtered_targets()
            self.view.display_targets_in_table(filtered_targets)
        except Exception as e:
            show_error(self.global_settings, "Error refreshing targets display", str(e))

    def show(self):
        self.view.show()

    def on_gene_selected(self, selected_text):
        """Handle gene selection signal"""
        try:
            # Extract locus tag from "locus_tag: gene_name" format
            locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
            self.global_settings.logger.debug(f"Loading sequence for locus tag: {locus_tag}")
            
            # Get gene sequence with padding using locus tag
            sequence_data = self.model.get_gene_sequence(locus_tag)
            if sequence_data:
                # Update gene viewer with sequence
                self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
                
                # Update location fields
                self.view.line_edit_start_location.setText(str(sequence_data['start']))
                self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                
                self.global_settings.logger.debug(f"Updated gene viewer with sequence of length: {len(sequence_data['sequence'])}")
            else:
                self.global_settings.logger.warning(f"No sequence data found for locus tag {locus_tag}")
                self.view.set_text_edit_gene_viewer("No sequence data available for this gene")
                self.view.line_edit_start_location.clear()
                self.view.line_edit_stop_location.clear()
                    
        except Exception as e:
            self.global_settings.logger.error(f"Error handling gene selection: {str(e)}")
            self.global_settings.logger.error(f"Stack trace: {traceback.format_exc()}")

    def highlight_targets_in_gene_viewer(self, targets_to_highlight=None):
        """Highlight selected targets in gene viewer"""
        try:
            self.global_settings.logger.debug("Starting highlight_gene_viewer")
            
            # Get selected targets if none provided
            if targets_to_highlight is None:
                targets_to_highlight = self.view.get_selected_targets()
            
            self.global_settings.logger.debug(f"Selected targets: {targets_to_highlight}")
            
            if not targets_to_highlight:
                QMessageBox.warning(self.view, "No Selection", 
                                  "Please select targets to highlight in the gene viewer.")
                return

            # Get current gene sequence
            selected_text = self.view.combo_box_gene.currentText()
            locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
            
            sequence_data = self.model.get_gene_sequence(locus_tag)
            if not sequence_data or 'sequence' not in sequence_data:
                self.global_settings.logger.error("No sequence data available for highlighting")
                return
                
            sequence = sequence_data['sequence']
            
            # Sort targets by position for efficient highlighting
            highlights = []
            sequences_found = 0
            total_sequences = len(targets_to_highlight)
            
            for target in targets_to_highlight:
                self.global_settings.logger.debug(f"Processing target: {target}")
                sequence_to_find = target['sequence']
                strand = target['strand']
                
                # For negative strand, we need to use reverse complement
                if strand == '-':
                    sequence_to_find = str(Seq(sequence_to_find).reverse_complement())
                    self.global_settings.logger.debug(f"Reverse complemented sequence: {sequence_to_find}")
                
                # Search for the sequence in the gene viewer text
                sequence_upper = sequence.upper()
                target_upper = sequence_to_find.upper()
                
                self.global_settings.logger.debug(f"Searching for sequence: {target_upper}")
                
                # Find all occurrences
                pos = sequence_upper.find(target_upper)
                if pos != -1:
                    self.global_settings.logger.debug(f"Found sequence at position: {pos}")
                    color = 'red' if strand == '-' else 'green'
                    highlights.append((pos, len(sequence_to_find), color))
                    sequences_found += 1
                else:
                    self.global_settings.logger.debug(f"Sequence not found: {target_upper}")

            # Only show warning if NO sequences were found
            if sequences_found == 0:
                self.global_settings.logger.warning("No sequences could be highlighted")
                QMessageBox.warning(self.view, "Highlighting Failed", 
                                  "Could not highlight any of the selected sequences in the current gene view.")
                return

            self.global_settings.logger.debug(f"Found {sequences_found} out of {total_sequences} sequences to highlight")

            # Build highlighted sequence
            result = []
            last_pos = 0
            for pos, length, color in sorted(highlights):  # Sort highlights by position
                result.append(sequence[last_pos:pos])
                result.append(f"<span style='background-color: {color};'>")
                result.append(sequence[pos:pos+length])
                result.append("</span>")
                last_pos = pos + length
            
            result.append(sequence[last_pos:])
            highlighted_sequence = ''.join(result)
            
            # Update the view with highlighted sequence
            self.view.update_gene_viewer(highlighted_sequence)
            self.global_settings.logger.debug(f"Successfully highlighted {sequences_found} sequences")
            
        except Exception as e:
            self.global_settings.logger.error(f"Error highlighting targets: {str(e)}")
            self.global_settings.logger.error(f"Stack trace: {traceback.format_exc()}")

    def update_scores(self, scores, algorithm):
        """Update the table with new scores from alternative scoring methods"""
        try:
            # Get current table headers
            headers = self.view.get_table_headers()
            
            # Get selected rows
            selected_rows = sorted(set(index.row() for index in self.view.table_targets.selectedIndexes()))
            if not selected_rows:
                self.global_settings.logger.warning("No rows selected for scoring")
                return
                
            # Determine the position for the new column (after the "Score" column)
            score_index = headers.index("Score")
            desired_index = score_index + 1
            
            # Disable updates to prevent crashes
            self.view.table_targets.setUpdatesEnabled(False)
            
            try:
                # Add new column for algorithm if it doesn't exist
                if algorithm not in headers:
                    # Store current column count
                    current_cols = self.view.table_targets.columnCount()
                    
                    # Insert new column after Score
                    self.view.table_targets.insertColumn(desired_index)
                    
                    # Set header for new column
                    self.view.table_targets.setHorizontalHeaderItem(
                        desired_index,
                        QtWidgets.QTableWidgetItem(algorithm)
                    )
                    
                    # Move Off-Target and Details columns one position right
                    for row in range(self.view.table_targets.rowCount()):
                        # Move Off-Target
                        off_target_item = self.view.table_targets.takeItem(row, desired_index)
                        if off_target_item:
                            self.view.table_targets.setItem(row, desired_index + 1, off_target_item)
                        
                        # Move Details button
                        details_widget = self.view.table_targets.cellWidget(row, desired_index)
                        if details_widget:
                            self.view.table_targets.setCellWidget(row, desired_index + 1, details_widget)
                    
                    col_index = desired_index
                else:
                    col_index = headers.index(algorithm)
                
                # Update scores in the table for selected rows only
                for score_idx, row in enumerate(selected_rows):
                    if score_idx < len(scores) and scores[score_idx] != -1:
                        score_item = QtWidgets.QTableWidgetItem()
                        # Round to 2 decimal places
                        rounded_score = round(float(scores[score_idx]), 2)
                        score_item.setData(QtCore.Qt.ItemDataRole.EditRole, rounded_score)
                        self.view.table_targets.setItem(row, col_index, score_item)
                        
                        # Also update the target data to preserve score during filtering/sorting
                        if hasattr(self.view, '_all_results'):
                            self.view._all_results[row]['azimuth_score'] = rounded_score
                
                # Resize columns to fit new content
                self.view.table_targets.resizeColumnsToContents()
                
                self.global_settings.logger.debug(f"Updated scores for algorithm: {algorithm}")
                self.global_settings.logger.debug(f"Updated rows: {selected_rows}")
                
            finally:
                # Re-enable updates
                self.view.table_targets.setUpdatesEnabled(True)
                
        except Exception as e:
            self.global_settings.logger.error(f"Error updating scores: {str(e)}")
            raise
