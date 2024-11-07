from models.ScoringOptionsModel import ScoringOptionsModel
from views.ScoringOptionsView import ScoringOptionsView

class ScoringOptionsController:
    def __init__(self, global_settings, view_targets_controller):
        self.global_settings = global_settings
        self.view_targets_controller = view_targets_controller
        self.model = ScoringOptionsModel(global_settings)
        self.view = ScoringOptionsView(global_settings)
        
        # Connect signals
        self.view.fasta_selected.connect(self._on_fasta_selected)
        self.view.submit_clicked.connect(self._on_submit)

    def show(self):
        """Show the scoring options window"""
        self.view.show()

    def _on_fasta_selected(self, fasta_path):
        """Handle FASTA file selection"""
        # Get current chromosome from view targets
        current_gene = self.view_targets_controller.view.combo_box_gene.currentText()
        locus_tag = current_gene.split(': ')[0] if ': ' in current_gene else current_gene
        print(f"Getting gene data for locus tag: {locus_tag}")
        gene_data = self.view_targets_controller.model.get_gene_data(locus_tag)
        
        if not gene_data or 'info' not in gene_data:
            self.view.show_error("Error", "Could not get chromosome information for current gene")
            return
            
        # Load FASTA file
        success = self.model.load_fasta(fasta_path, gene_data['info']['chromosome'])
        if not success:
            self.view.show_error("Error", "Failed to load FASTA file")
            return

    def _on_submit(self):
        """Handle submit button click"""
        try:
            # Validate inputs
            fasta_path = self.view.get_fasta_path()
            if not fasta_path:
                self.view.show_error("Error", "Please select a FASTA file")
                return
                
            algorithm = self.view.get_selected_algorithm()
            if not algorithm:
                self.view.show_error("Error", "Please select a scoring algorithm")
                return
                
            # Get selected targets from view targets
            selected_targets = self.view_targets_controller.view.get_selected_targets()
            if not selected_targets:
                self.view.show_error("Error", "No targets selected")
                return
                
            # Score sequences
            scores, reject_list, guide_list = self.model.score_sequences(selected_targets, algorithm)
            
            if scores is None:
                self.view.show_error("Error", "Failed to score sequences")
                return
                
            # Report rejected sequences
            if reject_list:
                rejected_seqs = "\n".join([guide_list[i] for i in reject_list])
                self.view.show_info(
                    "Sequences Not Found",
                    f"The following sequences were not found and scored as -1:\n{rejected_seqs}"
                )
                
            # Update scores in view targets
            self.view_targets_controller.update_scores(scores, algorithm)
            
            self.view.close()
            
        except Exception as e:
            self.global_settings.logger.error(f"Error in scoring submission: {str(e)}")
            self.view.show_error("Error", f"Error processing scores: {str(e)}")
