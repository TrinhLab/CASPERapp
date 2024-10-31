import logging
from models.ViewTargetsModel import ViewTargetsModel
from views.ViewTargetsView import ViewTargetsView
from PyQt6.QtWidgets import QMessageBox
from utils.ui import show_error
import time
import traceback
import threading

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
        try:
            start_time = time.time()
            
            # Get available genes from the model
            genes = self.model.get_available_genes()
            if genes:
                # Update the gene combo box
                self.view.combo_box_gene.clear()
                self.view.combo_box_gene.addItems(genes)
                
                # Fetch first gene immediately
                first_gene = genes[0]
                gene_data = self.model.get_gene_data(first_gene)
                
                if gene_data:
                    # Update the gene viewer with sequence
                    self.view.set_text_edit_gene_viewer(gene_data['sequence'])
                    
                    # Update location fields if available
                    if 'info' in gene_data and 'feature_location' in gene_data['info']:
                        location = gene_data['info']['feature_location']
                        if ':' in location:
                            start, end = location.split(':')[0], location.split(':')[1].split('(')[0]
                            self.view.line_edit_start_location.setText(start)
                            self.view.line_edit_stop_location.setText(end)
                    
                    # Pre-fetch next few genes in background thread
                    def prefetch_genes():
                        for gene in genes[1:5]:  # Pre-fetch next 4 genes
                            self.model.get_gene_data(gene)
                        
                    threading.Thread(target=prefetch_genes, daemon=True).start()
                
            execution_time = time.time() - start_time
            self.global_settings.logger.debug(f"Loading gene viewer took: {execution_time:.2f} seconds")
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in load_gene_viewer: {str(e)}\n{traceback.format_exc()}")

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
            self.global_settings.logger.debug(f"Current gene: {current_gene}")
            
            gene_data = self.model.get_gene_data(current_gene)
            if not gene_data:
                self.global_settings.logger.error("No gene data found")
                QMessageBox.warning(self.view, "No Gene Data", 
                                  "Could not get gene data for highlighting.")
                return

            self.global_settings.logger.debug(f"Gene sequence length: {len(gene_data['sequence'])}")
            
            # Highlight the sequences
            if targets_to_highlight:
                self.global_settings.logger.debug("Attempting to highlight sequences")
                highlighted_sequence = self.model.highlight_targets_in_gene_viewer(targets_to_highlight)
                
                if highlighted_sequence:
                    self.global_settings.logger.debug("Successfully highlighted sequences")
                    self.global_settings.logger.debug(f"Highlighted sequence length: {len(highlighted_sequence)}")
                    self.view.update_gene_viewer(highlighted_sequence)
                else:
                    self.global_settings.logger.error("Failed to highlight sequences - returned None")
                    QMessageBox.warning(self.view, "Highlighting Failed", 
                                      "Could not highlight the selected sequences. They may not be found in the current gene view.")
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
        try:
            scoring_options = self.model.get_scoring_options()
            self.view.show_scoring_options_dialog(scoring_options)
            if self.view.scoring_options_accepted():
                new_options = self.view.get_scoring_options()
                self.model.set_scoring_options(new_options)
                self.refresh_targets_display()
        except Exception as e:
            show_error(self.global_settings, "Error showing scoring options", str(e))

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
