from controllers.ScoringOptionsController import ScoringOptionsController
from models.ViewTargetsModel import ViewTargetsModel
from views.ViewTargetsView import ViewTargetsView
from PyQt6.QtWidgets import QMessageBox
from utils.ui import show_error
from PyQt6 import QtWidgets
import traceback
from views.LoadingDialog import LoadingDialog
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QColor
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
        self.view.gene_selected.connect(self.on_gene_selected)
        
        self.view.check_box_filter_5_prime_g_sequences.stateChanged.connect(self.refresh_guides_display)
        self.view.spin_box_minimum_on_target_score.valueChanged.connect(self.refresh_guides_display)
        self.view.check_box_view_exons_only.stateChanged.connect(self._on_view_exons_changed)

    def _on_view_exons_changed(self, state):
        """Handle view exons only checkbox state change"""
        try:
            is_checked = self.view.check_box_view_exons_only.isChecked()
            self.model.set_view_exons_only(is_checked)
            self.refresh_gene_viewer()
        except Exception as e:
            self.logger.error(f"Error handling view exons change: {str(e)}")

    def _clear_viewer_state(self):
        """Clear all highlights and cursor states from the gene viewer"""
        try:
            if hasattr(self.view, 'dna_feature_viewer'):
                # Clear sequence viewer highlights and cursor
                self.view.dna_feature_viewer.sequence_viewer.clear_highlights()
                for nuc in self.view.dna_feature_viewer.sequence_viewer.nucleotides:
                    nuc.show_cursor = False
                    nuc.update()
                self.view.dna_feature_viewer.sequence_viewer.selection_active = False
                self.view.dna_feature_viewer.sequence_viewer.selection_start = None
                self.view.dna_feature_viewer.sequence_viewer.selection_end = None
                
                # Clear insertion zone cursor
                if hasattr(self.view.dna_feature_viewer, 'insertion_zone'):
                    if hasattr(self.view.dna_feature_viewer.insertion_zone, 'sequence_cursor'):
                        self.view.dna_feature_viewer.insertion_zone.sequence_cursor.hide()
                    self.view.dna_feature_viewer.insertion_zone.current_cursor_pos = None
                    self.view.dna_feature_viewer.insertion_zone.selection_start = None
                    self.view.dna_feature_viewer.insertion_zone.selection_end = None
        except Exception as e:
            self.logger.error(f"Error clearing viewer state: {str(e)}")

    def load_guides(self, selected_targets, organism, endonuclease, loading_dialog=None):
        try:
            self.organism = organism
            self.endonuclease = endonuclease
            self.selected_targets = selected_targets

            # Use existing loading dialog if provided, otherwise create new one
            using_existing_dialog = loading_dialog is not None
            if not loading_dialog:
                loading_dialog = LoadingDialog(self.view, "Loading guides...")
                loading_dialog.show()
                QApplication.processEvents()
            
            try:
                # Clear any existing highlights and cursor
                self._clear_viewer_state()
                
                loading_dialog.set_message("Loading guides...", 60)
                QApplication.processEvents()
                
                self.model.load_guides(selected_targets, organism, endonuclease)
                
                # Initialize endonuclease combo box
                loading_dialog.set_message("Setting up endonucleases...", 65)
                QApplication.processEvents()
                
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
                
                loading_dialog.set_message("Processing guides...", 70)
                QApplication.processEvents()
                
                guides = self.model.get_guides()
                
                loading_dialog.set_message("Updating display...", 80)
                QApplication.processEvents()
                
                # Get unique position names or gene IDs
                unique_entries = set()
                for target in selected_targets:
                    if 'feature_id' in target:
                        if "chromosome" in str(target['feature_id']):
                            unique_entries.add(target['feature_id'])
                        else:
                            locus_tag = target['feature_id']
                            gene_data = self.model.get_gene_data(locus_tag)
                            if gene_data and 'info' in gene_data:
                                gene_name = gene_data['info'].get('gene_name', '')
                                display_text = f"{locus_tag}: {gene_name}" if gene_name else locus_tag
                                unique_entries.add(display_text)
                
                # Convert set to list for combo box
                entries = list(unique_entries)
                self.logger.debug(f"Found {len(entries)} unique entries")
                
                if entries:
                    first_entry = entries[0]
                    
                    # Block signals during setup
                    self.view.combo_box_gene.blockSignals(True)
                    
                    # Update combo box
                    self.view.set_combo_box_gene(entries)
                    self.view.combo_box_gene.setCurrentIndex(0)
                    
                    # Display guides for the first entry
                    if "chromosome" in first_entry and "start:" in first_entry:
                        position_guides = [g for g in guides if g.get('feature_id') == first_entry]
                        self.view.display_guides_in_table(position_guides)
                        
                        # Parse position from the text
                        parts = first_entry.split(',')
                        chrom = parts[0].split('chromosome')[1].strip()
                        start = int(parts[1].split('start:')[1].strip())
                        end = int(parts[2].split('end:')[1].strip())
                        
                        # Update location fields - add 1 to start for display
                        self.view.line_edit_start_location.setText(str(start + 1))
                        self.view.line_edit_stop_location.setText(str(end))
                        
                        # Get sequence directly for position-based search
                        sequence = self.model._get_sequence_for_position(chrom, start, end)
                        if sequence:
                            # Single call to set_data
                            self.view.dna_feature_viewer.set_data(sequence, [], start)
                    else:
                        # Regular gene-based search
                        locus_tag = first_entry.split(': ')[0] if ': ' in first_entry else first_entry
                        gene_guides = [g for g in guides if str(g.get('feature_id', '')).strip().lower() == locus_tag.lower()]
                        self.view.display_guides_in_table(gene_guides)
                        
                        # Get sequence data and features
                        sequence_data = self.model.get_gene_sequence(locus_tag)
                        if sequence_data:
                            # Update location fields - add 1 to start for display
                            self.view.line_edit_start_location.setText(str(sequence_data['start'] + 1))
                            self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                            
                            # Get features and update viewer in a single call
                            features = self.model.get_features_for_gene(locus_tag)
                            self.view.dna_feature_viewer.set_data(sequence_data['sequence'], features, sequence_data['start'])
                    
                    # Now unblock signals after everything is set up
                    self.view.combo_box_gene.blockSignals(False)
                
                loading_dialog.set_progress(100)
                QApplication.processEvents()
                
            finally:
                # Only close the dialog if we created it
                if not using_existing_dialog:
                    loading_dialog.close()
                    QApplication.processEvents()
                
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
            
            # Get current annotation file
            current_annotation_file = self.settings.get_current_annotation_file()
            if not current_annotation_file:
                QtWidgets.QMessageBox.warning(
                    self.view,
                    "No Annotation File",
                    "Please select an annotation file before performing off-target analysis."
                )
                return
                
            # Verify annotation file exists
            annotation_path = os.path.join(self.settings.get_db_path(), 'GBFF', current_annotation_file)
            if not os.path.isfile(annotation_path):
                # Try without GBFF subdirectory
                annotation_path = os.path.join(self.settings.get_db_path(), current_annotation_file)
                if not os.path.isfile(annotation_path):
                    QtWidgets.QMessageBox.warning(
                        self.view,
                        "Invalid Annotation File",
                        f"Could not find annotation file at {annotation_path}"
                    )
                    return
            
            # Set initial parameters based on current organism/endonuclease
            parameters = {
                'organism': self.organism,
                'endonuclease': self.endonuclease,
                'guides': selected_guides,  # Pass the selected guides
                'annotation_file': current_annotation_file  # Add annotation file
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
            scores, details = results  
            
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

            # Get current gene/position
            current_gene = self.view.combo_box_gene.currentText()
            sequence = None
            
            # Check if this is a position-based search
            if "chromosome" in current_gene and "start:" in current_gene:
                # Parse position from the text
                try:
                    parts = current_gene.split(',')
                    chrom = parts[0].split('chromosome')[1].strip()
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly for position-based search
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        self.logger.debug(f"Got sequence of length {len(sequence)} for position-based search")
                    else:
                        raise ValueError(f"Could not get sequence for position {start}-{end} in chromosome {chrom}")
                    
                except Exception as e:
                    self.logger.error(f"Error parsing position or getting sequence: {str(e)}")
                    QMessageBox.warning(
                        self.view,
                        "Sequence Error",
                        f"Could not get sequence for the selected position: {str(e)}"
                    )
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
                sequence = sequence_data['sequence']

            self.logger.debug(f"Got sequence of length: {len(sequence)}")
            
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

            # Highlight the sequences
            if guides_to_highlight:
                self.logger.debug("Attempting to highlight sequences")
                self.highlight_guides_in_gene_viewer(guides_to_highlight)
            else:
                self.logger.error("No valid guides to highlight")
                QMessageBox.warning(self.view, "No Valid Guides", 
                                  "Could not get sequence information from the selected rows.")

        except Exception as e:
            self.logger.error(f"Error in highlight_gene_viewer: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
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
        """Change the displayed sequence range"""
        try:
            # Get new start and end positions
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

            # Get current gene/position
            current_gene = self.view.combo_box_gene.currentText()
            
            # Handle position-based search
            if "chromosome" in current_gene and "start:" in current_gene:
                # For position-based searches
                try:
                    parts = current_gene.split(',')
                    # Get full chromosome identifier instead of just the number
                    chrom = parts[0].split('chromosome')[1].strip()  # This will now keep the full identifier

                    # Get sequence for new range - subtract 1 from start for 0-based indexing
                    sequence = self.model._get_sequence_for_position(chrom, new_start - 1, new_end)
                    
                    if sequence:
                        # Update DNA viewer with sequence
                        self.view.dna_feature_viewer.set_data(sequence, [], new_start - 1)
                        
                        # Update the line edits with new positions (keep display as 1-based)
                        self.view.line_edit_start_location.setText(str(new_start))
                        self.view.line_edit_stop_location.setText(str(new_end))
                    else:
                        raise ValueError("Could not get sequence for new position")
                        
                except Exception as e:
                    QMessageBox.warning(
                        self.view,
                        "Sequence Error",
                        f"Could not get sequence for range {new_start}-{new_end} in chromosome {chrom}"
                    )
            else:
                # For gene-based searches
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                
                # Get gene data to get chromosome information
                gene_data = self.model.get_gene_data(locus_tag)
                
                if gene_data and 'info' in gene_data:
                    try:
                        # Get chromosome from gene data
                        chrom = gene_data['info']['chromosome']
                        
                        # Get sequence directly using _get_sequence_for_position
                        sequence = self.model._get_sequence_for_position(chrom, new_start - 1, new_end)
                        
                        if sequence:
                            # Update line edits (keep display as 1-based)
                            self.view.line_edit_start_location.setText(str(new_start))
                            self.view.line_edit_stop_location.setText(str(new_end))
                            
                            # Get features for this gene
                            features = self.model.get_features_for_gene(locus_tag)
                            
                            # Update DNA viewer with new sequence
                            self.view.dna_feature_viewer.set_data(sequence, features, new_start - 1)
                            
                            # Update guides display for new range
                            gene_guides = [g for g in self.model.guides 
                                         if str(g.get('feature_id', '')).strip().lower() == locus_tag.lower() and
                                         new_start <= int(g['location'].split('-')[0]) <= new_end]
                            self.view.display_guides_in_table(gene_guides)
                        else:
                            raise ValueError("Could not get sequence for the specified range")
                            
                    except ValueError as ve:
                        QMessageBox.warning(
                            self.view,
                            "Range Error",
                            f"Invalid range: {str(ve)}"
                        )
                else:
                    QMessageBox.warning(
                        self.view,
                        "Gene Error",
                        f"Could not get sequence data for gene {locus_tag}"
                    )
                
        except Exception as e:
            self.logger.error(f"Error changing indices: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Error", f"Could not change sequence range: {str(e)}")

    def reset_location(self):
        """Reset gene viewer to the original sequence and location"""
        try:
            # Get current gene/position
            current_gene = self.view.combo_box_gene.currentText()
            
            # Check if this is a position-based search
            if "chromosome" in current_gene and "start:" in current_gene:
                try:
                    # Parse position from the text
                    parts = current_gene.split(',')
                    chrom = parts[0].split('chromosome')[1].strip()  # Keep full chromosome ID
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Get sequence directly using model's method
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        # Update DNA viewer with sequence
                        self.view.dna_feature_viewer.set_data(sequence, [], start)
                        
                        # Update location fields - add 1 to start for display
                        self.view.line_edit_start_location.setText(str(start + 1))
                        self.view.line_edit_stop_location.setText(str(end))
                        
                        # Update guides display
                        position_guides = [g for g in self.model.guides 
                                         if start <= int(g['location'].split('-')[0]) <= end]
                        self.view.display_guides_in_table(position_guides)
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
                # For gene-based searches
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                
                # Get original gene sequence
                sequence_data = self.model.get_gene_sequence(locus_tag)
                if sequence_data:
                    # Get features for this gene
                    features = self.model.get_features_for_gene(locus_tag)
                    
                    # Update DNA viewer with sequence and features
                    self.view.dna_feature_viewer.set_data(
                        sequence_data['sequence'], 
                        features, 
                        sequence_data['start']
                    )
                    
                    # Update location fields - add 1 to start for display
                    self.view.line_edit_start_location.setText(str(sequence_data['start'] + 1))
                    self.view.line_edit_stop_location.setText(str(sequence_data['end']))
                    
                    # Update guides display
                    gene_guides = [g for g in self.model.guides 
                                 if str(g.get('feature_id', '')).strip().lower() == locus_tag.lower()]
                    self.view.display_guides_in_table(gene_guides)
                else:
                    QMessageBox.warning(
                        self.view,
                        "Gene Error",
                        f"Could not get sequence data for gene {locus_tag}"
                    )
                    
        except Exception as e:
            self.logger.error(f"Error in reset_location: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            show_error(self.settings, "Reset Error", str(e))

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
                
                # Clear any existing highlights and cursor
                self._clear_viewer_state()
                
                self.view.display_guides_in_table(guides)
        except Exception as e:
            show_error(self.settings, "Error refreshing guides display", str(e))

    def show(self):
        self.view.show()

    def on_gene_selected(self, selected_text):
        """Handle gene selection signal"""
        try:
            # Create loading dialog
            loading_dialog = LoadingDialog(self.view, "Loading gene data...")
            loading_dialog.show()
            QApplication.processEvents()
            
            try:
                # Clear any existing highlights and cursor
                self._clear_viewer_state()
                
                # Load data in chunks
                loading_dialog.set_message("Loading sequence data...", 30)
                QApplication.processEvents()
                
                if "chromosome" in selected_text and "start:" in selected_text:
                    # Handle position-based search
                    parts = selected_text.split(',')
                    chrom = parts[0].split('chromosome')[1].strip()
                    start = int(parts[1].split('start:')[1].strip())
                    end = int(parts[2].split('end:')[1].strip())
                    
                    # Update location fields
                    self.view.line_edit_start_location.setText(str(start))
                    self.view.line_edit_stop_location.setText(str(end))
                    
                    # Get sequence directly for position-based search
                    sequence = self.model._get_sequence_for_position(chrom, start, end)
                    if sequence:
                        # Update DNA viewer with sequence
                        self.view.dna_feature_viewer.set_data(sequence, [], start)
                        self.logger.debug(f"Updated gene viewer with sequence of length {len(sequence)}")
                    else:
                        self.logger.error(f"Could not get sequence for position {start}-{end} in chromosome {chrom}")
                    
                    # Filter guides for this position
                    position_guides = [g for g in self.model.guides if g.get('feature_id') == selected_text]
                    self.view.display_guides_in_table(position_guides)
                    
                else:
                    # Regular gene-based search
                    locus_tag = selected_text.split(': ')[0] if ': ' in selected_text else selected_text
                    gene_data = self.model.get_gene_data(locus_tag)
                    
                    if gene_data and 'sequence' in gene_data and 'info' in gene_data:
                        self.view.line_edit_start_location.setText(str(gene_data['info']['start'] + 1))
                        self.view.line_edit_stop_location.setText(str(gene_data['info']['end']))
                        
                        loading_dialog.set_message("Updating display...", 80)
                        QApplication.processEvents()
                        
                        # Get features and update viewer in a single call
                        features = self.model.get_features_for_gene(locus_tag)
                        self.view.dna_feature_viewer.set_data(gene_data['sequence'], features, gene_data['info']['start'])
                        
                        # Filter guides for this gene
                        gene_guides = [g for g in self.model.guides 
                                     if str(g.get('feature_id', '')).strip().lower() == locus_tag.lower()]
                        self.view.display_guides_in_table(gene_guides)
                    else:
                        self.logger.warning(f"No valid gene data found for locus tag: {locus_tag}")
                
            finally:
                loading_dialog.close()
                
        except Exception as e:
            self.logger.error(f"Error handling gene selection: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def highlight_guides_in_gene_viewer(self, guides_to_highlight):
        """Highlight selected guides in gene viewer"""
        try:
            # Get current sequence
            current_gene = self.view.combo_box_gene.currentText()
            sequence_data = None
            
            # Get sequence based on view type
            if "chromosome" in current_gene and "start:" in current_gene:
                parts = current_gene.split(',')
                chrom = parts[0].split('chromosome')[1].strip()
                start = int(parts[1].split('start:')[1].strip())
                end = int(parts[2].split('end:')[1].strip())
                sequence = self.model._get_sequence_for_position(chrom, start, end)
                if sequence:
                    sequence_data = {'sequence': sequence}
            else:
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                sequence_data = self.model.get_gene_sequence(locus_tag)
                
            if not sequence_data or 'sequence' not in sequence_data:
                self.logger.error("No sequence data available")
                return
            
            sequence = sequence_data['sequence']
            sequence_upper = sequence.upper()
            
            # Clear existing highlights
            self.view.dna_feature_viewer.sequence_viewer.clear_highlights()
            
            for guide in guides_to_highlight:
                try:
                    print(f"Guide: {guide}")
                    guide_sequence = guide['sequence']
                    strand = guide['strand']
                    print(f"Strand: {strand}")
                    
                    # For negative strand guides
                    if strand == '-':
                        print("Negative strand")
                        # Convert sequence to complement for negative strand search
                        print(f"Sequence: {sequence_upper}")
                        complement_sequence = ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                                                     'K': 'M', 'Y': 'R', 'R': 'Y', 'M': 'K', 
                                                     'S': 'S'}[base] for base in sequence_upper)
                        print(f"Complement sequence: {complement_sequence}")
                        target_sequence = guide_sequence.upper()
                        print(f"Target sequence: {target_sequence}")
                        target_sequence = target_sequence[::-1]  # Reverse the sequence
                        print(f"Reversed target sequence: {target_sequence}")
                        pos = complement_sequence.find(target_sequence)
                        print(f"Position: {pos}")
                        
                        if pos != -1:
                            color = QColor(255, 0, 0, 100)  # Red for negative strand
                            self.logger.debug(f"Found negative strand sequence at position {pos}")
                            
                            # For negative strand, use the position directly but indicate strand
                            self.view.dna_feature_viewer.sequence_viewer.highlight_sequence(
                                pos,
                                pos + len(guide_sequence) - 1,
                                color,
                                strand='-'
                            )
                        else:
                            self.logger.warning(f"Negative strand sequence {target_sequence} not found")
                    else:
                        # For positive strand guides
                        target_sequence = guide_sequence.upper()
                        pos = sequence_upper.find(target_sequence)
                        
                        if pos != -1:
                            color = QColor(0, 255, 0, 100)  # Green for positive strand
                            self.logger.debug(f"Found positive strand sequence at position {pos}")
                            
                            self.view.dna_feature_viewer.sequence_viewer.highlight_sequence(
                                pos,
                                pos + len(guide_sequence) - 1,
                                color,
                                strand='+'
                            )
                        else:
                            self.logger.warning(f"Positive strand sequence {target_sequence} not found")
                    
                except Exception as e:
                    self.logger.error(f"Error highlighting guide: {str(e)}")
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error in highlight_guides_in_gene_viewer: {str(e)}")
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
            if "chromosome" in current_gene and "start:" in current_gene:
                # For position-based searches
                try:
                    parts = current_gene.split(',')
                    chrom = parts[0].split('chromosome')[1].strip()  # Keep full chromosome ID
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

    def refresh_gene_viewer(self):
        """Refresh gene viewer with sequence and features"""
        try:
            current_gene = self.view.combo_box_gene.currentText()
            if not current_gene:
                return

            self.logger.debug("Refreshing gene viewer")
            is_exons_only = self.view.check_box_view_exons_only.isChecked()
            self.logger.debug(f"View exons only is: {is_exons_only}")

            # Get gene data
            if "chromosome" in current_gene and "start:" in current_gene:
                # Handle position-based search
                parts = current_gene.split(',')
                chrom = parts[0].split('chromosome')[1].strip()
                start = int(parts[1].split('start:')[1].strip())
                end = int(parts[2].split('end:')[1].strip())
                
                self.logger.debug(f"Getting sequence for position: {chrom}:{start}-{end}")
                sequence = self.model._get_sequence_for_position(chrom, start, end)
                
                if sequence:
                    # Get features for this region
                    features = self.model.get_features_for_region(chrom, start, end)
                    self.logger.debug(f"Got sequence of length {len(sequence)} and {len(features)} features")
                    
                    # Verify DNA viewer exists
                    if not hasattr(self.view, 'dna_feature_viewer'):
                        self.logger.error("DNA viewer not initialized!")
                        return
                    
                    # Update DNA viewer directly
                    self.view.dna_feature_viewer.set_data(sequence, features, start)
                    self.logger.debug("Updated DNA viewer with sequence data")
                else:
                    self.logger.error("Failed to get sequence for position")
            else:
                # Regular gene-based search
                locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
                self.logger.debug(f"Getting sequence for locus tag: {locus_tag}")
                sequence_data = self.model.get_gene_sequence(locus_tag)
                
                if sequence_data and 'sequence' in sequence_data:
                    # Get features for this gene
                    features = self.model.get_features_for_gene(locus_tag)
                    sequence = sequence_data['sequence']
                    start_pos = sequence_data['start']
                    
                    self.logger.debug(f"Got sequence of length {len(sequence)}")
                    self.logger.debug(f"First 50 chars: {sequence[:50]}")
                    self.logger.debug(f"Start position: {start_pos}")
                    self.logger.debug(f"Number of features: {len(features)}")
                    
                    # Verify DNA viewer exists
                    if not hasattr(self.view, 'dna_feature_viewer'):
                        self.logger.error("DNA viewer not initialized!")
                        return
                    
                    # Update DNA viewer directly
                    self.view.dna_feature_viewer.set_data(sequence, features, start_pos)
                    self.logger.debug("Updated DNA viewer with sequence data")
                else:
                    self.logger.warning("No sequence data available")
                    
        except Exception as e:
            self.logger.error(f"Error refreshing gene viewer: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
