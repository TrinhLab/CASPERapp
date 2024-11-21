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
import os

class ViewTargetsController:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.model = ViewTargetsModel(global_settings)
        self.view = ViewTargetsView(global_settings)
        self.logger = global_settings.get_logger()
        
        self.setup_connections()
        self.organism = ""
        self.endonuclease = ""
        self.selected_targets = None

    def setup_connections(self):
        self.view.push_button_off_target.clicked.connect(self.perform_off_target_analysis)
        self.view.push_button_cotargeting.clicked.connect(self.perform_cotargeting)
        self.view.push_button_highlight_guides.clicked.connect(self.highlight_gene_viewer)
        self.view.push_button_clear_guides.clicked.connect(self.clear_highlighted_guides)
        self.view.push_button_export_selected_grnas.clicked.connect(self.export_targets)
        self.view.push_button_scoring_options.clicked.connect(self.show_scoring_options)
        self.view.push_button_change_location.clicked.connect(self.change_indices)
        self.view.push_button_reset_location.clicked.connect(self.reset_location)
        self.view.check_box_select_all.stateChanged.connect(self.select_all)
        # self.view.combo_box_gene.currentIndexChanged.connect(self.display_gene_data)
        self.view.gene_selected.connect(self.on_gene_selected)
        
        self.view.check_box_filter_5_prime_g_sequences.stateChanged.connect(self.refresh_guides_display)
        self.view.spin_box_minimum_on_target_score.valueChanged.connect(self.refresh_guides_display)

    def load_guides(self, selected_targets, organism, endonuclease):
        try:
            self.organism = organism
            self.endonuclease = endonuclease
            self.selected_targets = selected_targets

            print(f"Loading guides for {organism} and {selected_targets} with {endonuclease}")
            
            self.model.load_guides(selected_targets, organism, endonuclease)
            
            # Get available endonucleases for this organism
            org_to_endo = self.settings.get_organism_to_endonuclease()
            if organism in org_to_endo:
                available_endos = org_to_endo[organism]
                self.view.combo_box_endonuclease.clear()
                self.view.combo_box_endonuclease.addItems(available_endos)
                
                # Set current endonuclease
                current_index = self.view.combo_box_endonuclease.findText(endonuclease)
                if current_index >= 0:
                    self.view.combo_box_endonuclease.setCurrentIndex(current_index)
                    
                self.view.combo_box_endonuclease.currentTextChanged.connect(self._on_endonuclease_changed)
            
            # Format gene names for display
            seen_positions = set()
            formatted_genes = []
            
            if selected_targets and selected_targets[0].get('feature_type') == 'Position':
                position_groups = {}
                for target in selected_targets:
                    position_name = target['feature_id']
                    if position_name not in position_groups:
                        position_groups[position_name] = target
                        formatted_genes.append(position_name)
                
                if formatted_genes:
                    first_guide = position_groups[formatted_genes[0]]
                    self.view.line_edit_start_location.setText(str(first_guide['start']))
                    self.view.line_edit_stop_location.setText(str(first_guide['end']))
                    
                    if 'gene_sequence' in first_guide:
                        self.view.set_text_edit_gene_viewer(first_guide['gene_sequence'])
            else:
                for target in selected_targets:
                    gene_name = target.get('feature_name')
                    feature_id = target.get('feature_id')
                    
                    if gene_name and feature_id and gene_name not in seen_positions:
                        seen_positions.add(gene_name)
                        formatted_genes.append(f"{feature_id}: {gene_name}")
            
            formatted_genes.sort()
            self.view.set_combo_box_gene(formatted_genes)
            
            guides = self.model.get_guides()
            self.view.display_guides_in_table(guides)
            
            # Trigger gene sequence retrieval for first entry
            if formatted_genes:
                first_gene = formatted_genes[0]
                self.on_gene_selected(first_gene)

        except Exception as e:
            self.logger.error(f"Error in load_guides: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error loading guides", str(e))

    def _on_endonuclease_changed(self, new_endonuclease):
        try:
            if new_endonuclease != self.endonuclease:
                self.logger.debug(f"Changing endonuclease from {self.endonuclease} to {new_endonuclease}")
                self.endonuclease = new_endonuclease
                
                # Check if this is a co-targeting endonuclease combination
                if '|' in new_endonuclease:
                    selected_endos = new_endonuclease.split('|')
                    # Rerun co-targeting logic with selected endonucleases
                    self.handle_cotargeting_result(selected_endos)
                else:
                    # Regular single endonuclease handling
                    updated_targets = []
                    for target in self.selected_targets:
                        new_target = target.copy()
                        new_target['endonuclease'] = new_endonuclease
                        updated_targets.append(new_target)
                    
                    self.logger.debug(f"Created {len(updated_targets)} updated targets for {new_endonuclease}")
                    
                    # Update model with new targets
                    self.model.load_guides(updated_targets, self.organism, new_endonuclease)
                    guides = self.model.get_guides()
                    
                    self.logger.debug(f"Got {len(guides)} guides from model")
                    
                    # Update display
                    self.view.display_guides_in_table(guides)
                    
        except Exception as e:
            self.logger.error(f"Error changing endonuclease: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error", f"Could not change endonuclease: {str(e)}")

    def load_gene_viewer(self):
        try:
            
            # Get selected gene from combo box
            selected_text = self.view.combo_box_gene.currentText()
            if not selected_text:
                self.logger.debug("No gene selected")
                return
            
            # Extract locus tag from "locus_tag: gene_name" format
            locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
            self.logger.debug(f"Loading sequence for locus tag: {locus_tag}")
            
            # Get gene sequence with padding
            sequence_data = self.model.get_gene_sequence(locus_tag)
            
            if sequence_data:
                # Update gene viewer with sequence
                self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
                
                # Update location fields
                self.view.line_edit_start_location.setText(str(sequence_data['start']))
                self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                
            else:
                self.logger.warning(f"No sequence data found for locus tag {locus_tag}")
                
        except Exception as e:
            self.logger.error(f"Error in load_gene_viewer: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def perform_off_target_analysis(self):
        """Launch off-target analysis for selected guides"""
        try:
            # Get selected guides
            selected_guides = self.view.get_selected_guides()
            if not selected_guides:
                QtWidgets.QMessageBox.warning(
                    self.view, 
                    "No Selection", 
                    "Please select guides for off-target analysis."
                )
                return

            # Create off-target controller if not exists
            if not hasattr(self, '_off_target_controller'):
                from controllers.OffTargetController import OffTargetController
                self._off_target_controller = OffTargetController(self.settings)
                
                # Connect to results signal
                self._off_target_controller.off_target_results_ready.connect(
                    self._handle_off_target_results
                )
            
            # Set initial parameters based on current organism/endonuclease
            parameters = {
                'organism': self.organism,
                'endonuclease': self.endonuclease,
                'guides': selected_guides  # Pass the selected guides
            }
            
            # Initialize analysis with parameters
            self._off_target_controller.initialize_analysis(parameters)
            
            # Show and bring window to front
            self._off_target_controller.show()
            
        except Exception as e:
            self.logger.error(f"Error launching off-target analysis: {str(e)}")
            show_error(self.settings, "Error", f"Could not launch off-target analysis: {str(e)}")

    def _handle_off_target_results(self, results):
        """Handle off-target analysis results"""
        try:
            scores, details = results  # Unpack the tuple of results
            
            # Get current table headers
            headers = self.view.get_table_headers()
            
            # Find Score column index
            score_index = headers.index("Score")
            
            # Add Off-Target column if it doesn't exist
            if "Off-Target" not in headers:
                self.view.table_guides.insertColumn(score_index + 1)
                self.view.table_guides.setHorizontalHeaderItem(
                    score_index + 1,
                    QtWidgets.QTableWidgetItem("Off-Target")
                )
            
            off_target_index = headers.index("Off-Target") if "Off-Target" in headers else score_index + 1
            
            # Update the view with results and details
            self.view.update_off_target_details(scores, details)
            
            self.logger.debug(f"Updated off-target scores in table with {len(details) if details else 0} detailed results")
            
        except Exception as e:
            self.logger.error(f"Error handling off-target results: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error", f"Could not update off-target scores: {str(e)}")

    def perform_cotargeting(self):
        """Launch co-targeting analysis"""
        try:
            # Get current endonuclease choices
            current_items = [self.view.combo_box_endonuclease.itemText(i) 
                            for i in range(self.view.combo_box_endonuclease.count())]
            
            if len(current_items) <= 1:
                QtWidgets.QMessageBox.warning(
                    self.view,
                    "Not Enough Endonucleases",
                    "There are not enough endonucleases with this organism. At least 2 endonucleases are required for this function."
                )
                return

            # Get cotargeting controller and launch, passing self as view_targets_controller
            cotargeting_controller = self.settings.get_cotargeting_window(self)
            cotargeting_controller.launch(
                endo_choices=current_items,
                org_name=self.organism
            )

        except Exception as e:
            self.logger.error(f"Error in perform_cotargeting: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error", f"Could not launch co-targeting: {str(e)}")

    def highlight_gene_viewer(self):
        """Highlight selected guides in gene viewer"""
        try:
            self.logger.debug("Starting highlight_gene_viewer")
            
            # Get selected guides
            selected_rows = self.view.get_selected_guides()
            self.logger.debug(f"Selected guides: {selected_rows}")
            
            if not selected_rows:
                QMessageBox.warning(self.view, "No Selection", 
                                  "Please select guides to highlight in the gene viewer.")
                return

            # Convert table selections to the format expected by the model
            guides_to_highlight = []
            for guide in selected_rows:
                guide_info = {
                    'location': guide['location'],
                    'sequence': guide['sequence'],
                    'strand': guide['strand']
                }
                guides_to_highlight.append(guide_info)
                self.logger.debug(f"Guide to highlight: {guide_info}")

            # Get current gene sequence
            current_gene = self.view.combo_box_gene.currentText()
            
            # Check if this is a position-based search
            if "chrom" in current_gene and "start:" in current_gene:
                # Parse position from the text (format: "chrom X, start: Y, end: Z")
                try:
                    parts = current_gene.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly from model's _get_sequence_for_position
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if not sequence:
                        raise ValueError("Could not get sequence for position")
                        
                    sequence_data = {
                        'sequence': sequence,
                        'start': start,
                        'end': end
                    }
                    self.logger.debug(f"Got position-based sequence of length: {len(sequence)}")
                except Exception as e:
                    self.logger.error(f"Error getting position sequence: {str(e)}")
                    QMessageBox.warning(self.view, "Error", 
                                      "Could not get sequence for the specified position.")
                    return
            else:
                # Regular gene-based search
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                self.logger.debug(f"Getting sequence for locus tag: {locus_tag}")
                
                # Get gene sequence with padding
                sequence_data = self.model.get_gene_sequence(locus_tag)
                if not sequence_data or 'sequence' not in sequence_data:
                    self.logger.error("No sequence data found")
                    QMessageBox.warning(self.view, "No Gene Data", 
                                      "Could not get gene sequence for highlighting.")
                    return

            self.logger.debug(f"Gene sequence length: {len(sequence_data['sequence'])}")
            
            # Highlight the sequences
            if guides_to_highlight:
                self.logger.debug("Attempting to highlight sequences")
                self.highlight_guides_in_gene_viewer(guides_to_highlight)
            else:
                self.logger.error("No valid guides to highlight")
                QMessageBox.warning(self.view, "No Valid Guides", 
                                  "Could not get sequence information from the selected rows.")

        except Exception as e:
            self.logger.error(f"Error in highlight_gene_viewer: {str(e)}\n{traceback.format_exc()}")
            show_error(self.settings, "Error highlighting gene viewer", str(e))

    def export_targets(self):
        try:
            # Get selected guides
            selected_guides = self.view.get_selected_guides()
            if not selected_guides:
                QtWidgets.QMessageBox.warning(
                    self.view,
                    "No Selection",
                    "Please select guides to export."
                )
                return

            export_controller = self.settings.get_export_selected_grnas_window()
            export_controller.show_dialog(selected_guides, "View Targets")
        except Exception as e:
            self.logger.error(f"Error in export_targets: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Export Error", str(e))

    def show_scoring_options(self):
        try:
            # Create scoring options controller if not exists
            if not hasattr(self, '_scoring_options_controller'):
                # Create controller with self as view_targets_controller
                self._scoring_options_controller = ScoringOptionsController(
                    global_settings=self.settings,
                    view_targets_controller=self
                )
                
            self._scoring_options_controller.show()
            
        except Exception as e:
            self.logger.error(f"Error showing scoring options: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error", f"Could not show scoring options: {str(e)}")

    def change_indices(self):
        """Change the start and end positions for gene viewer"""
        try:
            # Make sure gene viewer has content
            if not self.view.text_edit_gene_viewer.toPlainText():
                QMessageBox.warning(
                    self.view,
                    "Gene Viewer Error",
                    "Gene Viewer display is empty! Please ensure there is sequence data to view."
                )
                return

            # Get current gene/position info
            current_gene = self.view.combo_box_gene.currentText()
            
            try:
                new_start = int(self.view.line_edit_start_location.text())
                new_end = int(self.view.line_edit_stop_location.text())
            except ValueError:
                QMessageBox.warning(
                    self.view,
                    "Invalid Input",
                    "Please enter valid integer values for start and end positions."
                )
                return

            # Validate start and end positions
            if new_start <= 0 or new_end <= 0:
                QMessageBox.warning(
                    self.view,
                    "Invalid location indices",
                    "Location indices cannot be negative or zero! Please set values larger than 0."
                )
                return

            if new_start >= new_end:
                QMessageBox.warning(
                    self.view,
                    "Invalid location indices",
                    "Start location must be less than stop location."
                )
                return

            if abs(new_start - new_end) > 50000:
                QMessageBox.warning(
                    self.view,
                    "Sequence Too Long",
                    "The sequence is too long! Please choose indices that will make the sequence less than 50,000!"
                )
                return

            # Get sequence for new range
            if "chrom" in current_gene and "start:" in current_gene:
                # For position-based searches
                try:
                    parts = current_gene.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    
                    # Get sequence for new range
                    sequence = self.model._get_sequence_for_position(chrom, new_start, new_end)
                    
                    if sequence:
                        self.view.set_text_edit_gene_viewer(sequence)
                        # Update the line edits with new positions
                        self.view.line_edit_start_location.setText(str(new_start))
                        self.view.line_edit_stop_location.setText(str(new_end))
                    else:
                        raise ValueError("Could not get sequence for new position")
                        
                except Exception as e:
                    QMessageBox.warning(
                        self.view,
                        "Position Error",
                        f"Error changing position: {str(e)}"
                    )
                    return
            else:
                # For feature-based searches
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                gene_data = self.model.get_gene_data(locus_tag)
                
                if not gene_data or 'info' not in gene_data:
                    QMessageBox.warning(
                        self.view,
                        "Gene Data Error",
                        "Could not get gene data for the current selection."
                    )
                    return
                    
                # Get new sequence for the range
                sequence_data = self.model.get_gene_sequence_for_range(locus_tag, new_start, new_end)
                if sequence_data:
                    self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
                    # Update the line edits with new positions
                    self.view.line_edit_start_location.setText(str(new_start))
                    self.view.line_edit_stop_location.setText(str(new_end))
                else:
                    QMessageBox.warning(
                        self.view,
                        "Sequence Error",
                        "Could not get sequence for the specified range."
                    )
                    return

        except Exception as e:
            self.logger.error(f"Error in change_indices: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error changing indices", str(e))

    def reset_location(self):
        """Reset gene viewer to the original sequence and location"""
        try:
            # Get current gene/position
            current_gene = self.view.combo_box_gene.currentText()
            
            # Check if this is a position-based search
            if "chrom" in current_gene and "start:" in current_gene:
                try:
                    # Parse position from the text (format: "chrom X, start: Y, end: Z")
                    parts = current_gene.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly using model's method
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        # Update gene viewer with sequence
                        self.view.set_text_edit_gene_viewer(sequence)
                        
                        # Update location fields
                        self.view.line_edit_start_location.setText(str(start))
                        self.view.line_edit_stop_location.setText(str(end))
                    else:
                        raise ValueError("Could not get sequence for position")
                        
                except Exception as e:
                    self.logger.error(f"Error resetting position: {str(e)}")
                    QMessageBox.warning(
                        self.view,
                        "Position Error",
                        f"Error resetting position: {str(e)}"
                    )
                    return
            else:
                # For feature-based searches
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                sequence_data = self.model.get_gene_sequence(locus_tag)
                
                if sequence_data:
                    # Update gene viewer with sequence
                    self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
                    
                    # Update location fields
                    self.view.line_edit_start_location.setText(str(sequence_data['start']))
                    self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                else:
                    self.logger.warning(f"No sequence data found for locus tag {locus_tag}")
                    QMessageBox.warning(
                        self.view,
                        "Gene Data Error",
                        "Could not get gene sequence for resetting location."
                    )
                    return
                
        except Exception as e:
            self.logger.error(f"Error in reset_location: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error resetting location", str(e))

    def select_all(self, state):
        try:
            self.view.select_all_guides(state == 2)  # 2 corresponds to Qt.Checked
        except Exception as e:
            show_error(self.settings, "Error selecting all guides", str(e))

    def display_gene_data(self, gene_name):
        try:
            gene_data = self.model.get_gene_data(gene_name)
            print(f"Gene data: {gene_data}")
            if gene_data and gene_data['sequence']:
                self.view.set_text_edit_gene_viewer(gene_data['sequence'])
                self.view.line_edit_start_location.setText(str(gene_data['start']))
                self.view.line_edit_stop_location.setText(str(gene_data['end']))
            else:
                self.view.set_text_edit_gene_viewer("No sequence data available for this gene")
                self.view.line_edit_start_location.clear()
                self.view.line_edit_stop_location.clear()
        except Exception as e:
            show_error(self.settings, "Error displaying gene data", str(e))

    def refresh_guides_display(self):
        """Refresh the guides display when filters change"""
        try:
            if hasattr(self, 'selected_targets'):
                # Get current guides from model
                self.model.load_guides(self.selected_targets, self.organism, self.endonuclease)
                guides = self.model.get_guides()
                self.view.display_guides_in_table(guides)
        except Exception as e:
            show_error(self.settings, "Error refreshing guides display", str(e))

    def show(self):
        self.view.show()

    def on_gene_selected(self, selected_text):
        """Handle gene selection signal"""
        try:
            self.logger.debug(f"Gene selection changed to: {selected_text}")
            
            # Check if this is a position-based search
            if "chrom" in selected_text and "start:" in selected_text:
                try:
                    # Parse position from the text (format: "chrom X, start: Y, end: Z")
                    parts = selected_text.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly using _get_sequence_for_position
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        # Update gene viewer with sequence
                        self.view.set_text_edit_gene_viewer(sequence)
                        
                        # Update location fields
                        self.view.line_edit_start_location.setText(str(start))
                        self.view.line_edit_stop_location.setText(str(end))
                        
                        self.logger.debug(f"Updated position view with sequence of length: {len(sequence)}")
                        
                        # Filter guides for this position
                        position_guides = [g for g in self.model.guides 
                                         if g.get('feature_id') == selected_text]
                        self.view.display_guides_in_table(position_guides)
                    else:
                        self.logger.warning(f"No sequence found for position {chrom}:{start}-{end}")
                        self.view.set_text_edit_gene_viewer("No sequence data available for this position")
                        self.view.line_edit_start_location.clear()
                        self.view.line_edit_stop_location.clear()
                except Exception as e:
                    self.logger.error(f"Error handling position selection: {str(e)}")
                    self.logger.error(f"Stack trace: {traceback.format_exc()}")
            else:
                # Regular gene-based search
                locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
                self.logger.debug(f"Loading sequence for locus tag: {locus_tag}")
                
                # Get gene sequence with padding using locus tag
                sequence_data = self.model.get_gene_sequence(locus_tag)
                if sequence_data:
                    # Update gene viewer with sequence
                    self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
                    
                    # Update location fields
                    self.view.line_edit_start_location.setText(str(sequence_data['start']))
                    self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                    
                    self.logger.debug(f"Updated gene viewer with sequence of length: {len(sequence_data['sequence'])}")
                    
                    # Filter guides for this gene
                    gene_guides = [g for g in self.model.guides 
                                  if str(g.get('feature_id', '')).strip().lower() == locus_tag.lower()]
                    self.view.display_guides_in_table(gene_guides)
                else:
                    self.logger.warning(f"No sequence data found for locus tag {locus_tag}")
                    self.view.set_text_edit_gene_viewer("No sequence data available for this gene")
                    self.view.line_edit_start_location.clear()
                    self.view.line_edit_stop_location.clear()
                    
        except Exception as e:
            self.logger.error(f"Error handling gene selection: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def highlight_guides_in_gene_viewer(self, guides_to_highlight=None):
        """Highlight selected guides in gene viewer"""
        try:
            self.logger.debug("Starting highlight_gene_viewer")
            
            if guides_to_highlight is None:
                guides_to_highlight = self.view.get_selected_guides()
            
            self.logger.debug(f"Selected guides: {guides_to_highlight}")
            
            if not guides_to_highlight:
                QMessageBox.warning(self.view, "No Selection", 
                                  "Please select guides to highlight in the gene viewer.")
                return

            # Get current gene sequence
            selected_text = self.view.combo_box_gene.currentText()
            
            # For position-based searches, get sequence directly from model
            if "chrom" in selected_text and "start:" in selected_text:
                try:
                    # Parse position from the text (format: "chrom X, start: Y, end: Z")
                    parts = selected_text.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly from FindTargetsModel
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if not sequence:
                        raise ValueError("Could not get sequence for position")
                    
                    self.logger.debug(f"Got sequence of length {len(sequence)} for position-based search")
                    
                except Exception as e:
                    self.logger.error(f"Error parsing position or getting sequence: {str(e)}")
                    return
            else:
                # Regular gene-based search
                locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
                sequence_data = self.model.get_gene_sequence(locus_tag)
                if not sequence_data or 'sequence' not in sequence_data:
                    self.logger.error("No sequence data available for highlighting")
                    return
                sequence = sequence_data['sequence']
                
            # Process highlights
            highlights = []
            sequences_found = 0
            total_sequences = len(guides_to_highlight)
            
            for guide in guides_to_highlight:
                self.logger.debug(f"Processing guide: {guide}")
                sequence_to_find = guide['sequence']
                strand = guide['strand']
                
                if strand == '-':
                    sequence_to_find = str(Seq(sequence_to_find).reverse_complement())
                    self.logger.debug(f"Reverse complemented sequence: {sequence_to_find}")
                
                sequence_upper = sequence.upper()
                target_upper = sequence_to_find.upper()
                
                self.logger.debug(f"Searching for sequence: {target_upper}")
                
                pos = sequence_upper.find(target_upper)
                if pos != -1:
                    self.logger.debug(f"Found sequence at position: {pos}")
                    color = 'red' if strand == '-' else 'green'
                    highlights.append((pos, len(sequence_to_find), color))
                    sequences_found += 1
                else:
                    self.logger.debug(f"Sequence not found: {target_upper}")

            if sequences_found == 0:
                self.logger.warning("No sequences could be highlighted")
                QMessageBox.warning(self.view, "Highlighting Failed", 
                                  "Could not highlight any of the selected sequences in the current gene view.")
                return

            # Build highlighted sequence
            result = []
            last_pos = 0
            for pos, length, color in sorted(highlights):
                result.append(sequence[last_pos:pos])
                result.append(f"<span style='background-color: {color};'>")
                result.append(sequence[pos:pos+length])
                result.append("</span>")
                last_pos = pos + length
            
            result.append(sequence[last_pos:])
            highlighted_sequence = ''.join(result)
            
            # Update the view with highlighted sequence
            self.view.update_gene_viewer(highlighted_sequence)
            self.logger.debug(f"Successfully highlighted {sequences_found} sequences")
            
        except Exception as e:
            self.logger.error(f"Error highlighting guides: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def update_scores(self, scores, algorithm):
        """Update the table with new scores from alternative scoring methods"""
        try:
            # Get current table headers
            headers = self.view.get_table_headers()
            
            # Get selected rows
            selected_rows = sorted(set(index.row() for index in self.view.table_guides.selectedIndexes()))
            if not selected_rows:
                self.logger.warning("No rows selected for scoring")
                return
                
            # Determine the position for the new column (after the "Score" column)
            score_index = headers.index("Score")
            desired_index = score_index + 1
            
            # Disable updates to prevent crashes
            self.view.table_guides.setUpdatesEnabled(False)
            
            try:
                # Add new column for algorithm if it doesn't exist
                if algorithm not in headers:
                    
                    # Insert new column after Score
                    self.view.table_guides.insertColumn(desired_index)
                    
                    # Set header for new column
                    self.view.table_guides.setHorizontalHeaderItem(
                        desired_index,
                        QtWidgets.QTableWidgetItem(algorithm)
                    )
                    
                    # Move Off-Target and Details columns one position right
                    for row in range(self.view.table_guides.rowCount()):
                        # Move Off-Target
                        off_target_item = self.view.table_guides.takeItem(row, desired_index)
                        if off_target_item:
                            self.view.table_guides.setItem(row, desired_index + 1, off_target_item)
                        
                        # Move Details button
                        details_widget = self.view.table_guides.cellWidget(row, desired_index)
                        if details_widget:
                            self.view.table_guides.setCellWidget(row, desired_index + 1, details_widget)
                    
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
                        self.view.table_guides.setItem(row, col_index, score_item)
                        
                        # Also update the guide data to preserve score during filtering/sorting
                        if hasattr(self.view, '_all_results'):
                            self.view._all_results[row]['azimuth_score'] = rounded_score
                
                # Resize columns to fit new content
                self.view.table_guides.resizeColumnsToContents()
                
                self.logger.debug(f"Updated scores for algorithm: {algorithm}")
                self.logger.debug(f"Updated rows: {selected_rows}")
                
            finally:
                # Re-enable updates
                self.view.table_guides.setUpdatesEnabled(True)
                
        except Exception as e:
            self.logger.error(f"Error updating scores: {str(e)}")
            raise

    def clear_highlighted_guides(self):
        """Clear highlights in gene viewer and unselect rows in target table"""
        try:
            # Get current gene/position
            current_gene = self.view.combo_box_gene.currentText()
            
            # Reset gene viewer to original sequence
            if "chrom" in current_gene and "start:" in current_gene:
                # For position-based searches
                try:
                    parts = current_gene.split(',')
                    chrom = int(parts[0].split('chrom')[1].strip())
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        self.view.set_text_edit_gene_viewer(sequence)
                except Exception as e:
                    self.logger.error(f"Error resetting position sequence: {str(e)}")
            else:
                # For feature-based searches
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                sequence_data = self.model.get_gene_sequence(locus_tag)
                if sequence_data and 'sequence' in sequence_data:
                    self.view.set_text_edit_gene_viewer(sequence_data['sequence'])
            
            # Clear selected rows in table
            self.view.table_guides.clearSelection()
            
            # Uncheck select all checkbox
            self.view.check_box_select_all.setChecked(False)
            
            self.logger.debug("Cleared highlighted guides and selections")
            
        except Exception as e:
            self.logger.error(f"Error clearing highlighted guides: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error clearing guides", str(e))

    def handle_cotargeting_result(self, selected_endos):
        try:
            # Format combined endonuclease string
            endo_string = "|".join(selected_endos)
            
            # Add combined choice to endonuclease combo box
            current_items = [self.view.combo_box_endonuclease.itemText(i) 
                            for i in range(self.view.combo_box_endonuclease.count())]
            
            if endo_string not in current_items:
                self.view.combo_box_endonuclease.insertItem(0, endo_string)
                
            # Set as current selection
            self.view.combo_box_endonuclease.setCurrentText(endo_string)

            # Define PAM compatibility rules
            pam_rules = {
                'SpCas9': 'NGG',
                'SaCas9': 'NNGRRT',
                'NmCas9': 'NNNNGATT',
                'St1Cas9': 'NNAGAAW',
                'St3Cas9': 'NGGNG',
                'TdCas9': 'NAAAAC',
                'CjCas9': 'NNNNRYAC'
            }
            
            # Function to check if a PAM sequence matches the required pattern
            def is_pam_compatible(pam_seq, pattern):
                if len(pam_seq) != len(pattern):
                    return False
                
                for p, t in zip(pam_seq.upper(), pattern.upper()):
                    if t == 'N':
                        continue
                    elif t == 'R' and p not in 'AG':
                        return False
                    elif t == 'Y' and p not in 'CT':
                        return False
                    elif t == 'W' and p not in 'AT':
                        return False
                    elif t != p and t != 'N':
                        return False
                return True
            
            # Get guides for each individual endonuclease
            all_guides = []
            for endo in selected_endos:
                # Update guides with current endonuclease
                updated_guides = []
                for guide in self.selected_targets:
                    new_guide = guide.copy()
                    new_guide['endonuclease'] = endo
                    
                    # Extract chromosome number from full identifier
                    if 'chromosome' in new_guide:
                        chrom_id = new_guide['chromosome']
                        if isinstance(chrom_id, str) and '.' in chrom_id:
                            chrom_num = chrom_id.split('.')[-1]
                            new_guide['chromosome'] = chrom_num
                    
                    updated_guides.append(new_guide)
                
                # Load guides for this endonuclease
                self.model.load_guides(updated_guides, self.organism, endo)
                guides = self.model.get_guides()
                all_guides.extend(guides)
            
            # Group guides by sequence and collect their PAMs
            sequence_groups = {}
            for guide in all_guides:
                seq = guide['sequence']
                if seq not in sequence_groups:
                    sequence_groups[seq] = {'guides': [], 'pams': set(), 'endos': set()}
                sequence_groups[seq]['guides'].append(guide)
                sequence_groups[seq]['pams'].add(guide['pam'])
                sequence_groups[seq]['endos'].add(guide['endonuclease'])
            
            # Filter for guides that have compatible PAMs across all endonucleases
            cotargeted_guides = []
            for seq_info in sequence_groups.values():
                if len(seq_info['endos']) == len(selected_endos):
                    # Get the most stringent PAM from the guides
                    guide_pams = list(seq_info['pams'])
                    
                    # Sort PAMs by length (longer PAMs are typically more stringent)
                    guide_pams.sort(key=len, reverse=True)
                    stringent_pam = guide_pams[0]
                    
                    # Check if this PAM is compatible with all selected endonucleases
                    is_compatible = True
                    for endo in selected_endos:
                        for cas9_type, pam_pattern in pam_rules.items():
                            if cas9_type in endo:
                                if not is_pam_compatible(stringent_pam, pam_pattern):
                                    is_compatible = False
                                    break
                        if not is_compatible:
                            break
                    
                    if is_compatible:
                        # Use the guide with the most stringent PAM
                        for guide in seq_info['guides']:
                            if guide['pam'] == stringent_pam:
                                combined_guide = guide.copy()
                                combined_guide['endonuclease'] = endo_string
                                cotargeted_guides.append(combined_guide)
                                break
            
            # Update display with co-targeted guides
            self.view.display_guides_in_table(cotargeted_guides)
            
            self.logger.debug(f"Found {len(cotargeted_guides)} co-targeted guides with compatible PAMs")
            
        except Exception as e:
            self.logger.error(f"Error handling co-targeting result: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Co-targeting Error", str(e))
