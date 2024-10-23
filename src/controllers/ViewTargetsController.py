from models.ViewTargetsModel import ViewTargetsModel
from views.ViewTargetsView import ViewTargetsView
from PyQt6.QtWidgets import QMessageBox
from utils.ui import show_error

class ViewTargetsController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = ViewTargetsModel(global_settings)
        self.view = ViewTargetsView(global_settings)
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
            self.organism = organism
            self.endonuclease = endonuclease
            self.model.load_targets(selected_targets, organism, endonuclease)
            targets = self.model.get_targets()
            self.view.display_targets_in_table(targets)
            self.view.set_combo_box_endonuclease([endonuclease])
            self.load_gene_viewer(selected_targets)
        except Exception as e:
            show_error(self.global_settings, "Error loading targets", str(e))

    def load_gene_viewer(self, selected_targets):
        try:
            genes = self.model.available_genes
            self.view.set_combo_box_gene(genes)
            if genes:
                self.display_gene_data(genes[0])
            else:
                self.view.set_text_edit_gene_viewer("No gene data available")
                self.view.line_edit_start_location.clear()
                self.view.line_edit_stop_location.clear()
        except Exception as e:
            show_error(self.global_settings, "Error loading gene viewer", str(e))

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
            selected_targets = self.view.get_selected_targets()
            if not selected_targets:
                QMessageBox.warning(self.view, "No Selection", "Please select targets to highlight in the gene viewer.")
                return
            highlighted_sequence = self.model.highlight_targets_in_gene_viewer(selected_targets)
            self.view.update_gene_viewer(highlighted_sequence)
        except Exception as e:
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
