from models.CSPRparser import CSPRparser
from models.HomeWindowModel import HomeWindowModel
from models.AnnotationParser import AnnotationParser
import os
from Bio import SeqIO
from Bio.Seq import Seq

class ViewTargetsModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.targets = []
        self.cspr_parser = None
        self.annotation_parser = AnnotationParser(global_settings)
        self.gene_sequence = ""
        self.highlighted_sequence = ""
        self.gene_info = {}
        self.available_genes = []
        self.filter_options = {}
        self.scoring_options = {}
        self.annotation_path = ""
        self.current_gene_start = 0
        self.current_gene_end = 0
        self.extended_sequence = ""
        self.chromosome = ""

    def load_targets(self, selected_targets, organism, endonuclease):
        org_files = self.get_organism_to_files()
        if organism in org_files and endonuclease in org_files[organism]:
            cspr_file = org_files[organism][endonuclease][0]
            cspr_path = os.path.join(self.global_settings.get_db_path(), cspr_file)
            self.cspr_parser = CSPRparser(cspr_path, self.global_settings.get_casper_info_path())
            
            self.targets = []
            self.available_genes = []
            for target in selected_targets:
                self.chromosome = target['chromosome']  # Store the chromosome
                chrom = target['chromosome']
                start, end = map(int, target['location'].split('-'))
                pos_tuple = (chrom, start, end)
                cspr_targets = self.cspr_parser.read_targets(target['feature_name'], pos_tuple, endonuclease)
                self.targets.extend(cspr_targets)
                if target['feature_name'] not in self.available_genes:
                    self.available_genes.append(target['feature_name'])
            
            if self.available_genes:
                self._load_gene_data(self.available_genes[0])

    def _load_gene_data(self, gene_name):
        annotation_files = self.get_annotation_files()
        if annotation_files:
            self.annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_files[0])
            self.annotation_parser.set_annotation_file(self.annotation_path)
            gene_data = self.annotation_parser.get_gene_data(gene_name)
            print("gene_data", gene_data)
            if gene_data:
                self.gene_sequence = gene_data.get('sequence', '')
                self.gene_info = gene_data.get('info', {})
                
                location = self.gene_info.get('feature_location', '')
                if location:
                    start, end, strand = self._parse_location(location)
                    print("strand", strand)
                    self.current_gene_start = start + 1
                    self.current_gene_end = end
                    
                    # Fetch the extended sequence
                    extended_start = max(0, start - 30)
                    extended_end = end + 30
                    extended_seq = self.get_sequence_from_annotation(self.chromosome, extended_start, extended_end)
                    
                    if strand == '-':
                        # Reverse complement the entire extended sequence
                        self.gene_sequence = str(Seq(self.gene_sequence).reverse_complement())
                        print("self.gene_sequence", self.gene_sequence)
                    
                    # Create the extended sequence with lowercase placeholders
                    left_placeholder = extended_seq[:30].lower()
                    right_placeholder = extended_seq[-30:].lower()
                    self.gene_sequence = extended_seq[30:-30]
                    self.extended_sequence = f"{left_placeholder}{self.gene_sequence}{right_placeholder}"
                
                self.global_settings.logger.debug(f"Loaded gene sequence for {gene_name}: {self.extended_sequence[:50]}...")
            else:
                self.global_settings.logger.warning(f"Gene {gene_name} not found in annotation file")

    def _parse_location(self, location):
        parts = location.split(':')
        start = parts[0]
        end_strand = parts[1]
        end = end_strand[:-3]  # Remove the last 3 characters (strand info)
        strand = end_strand[-2:-1]  # Get the strand info (last 2 characters, excluding the closing parenthesis)
        return int(start), int(end), strand

    def get_targets(self):
        return self.targets

    def perform_off_target_analysis(self, selected_targets):
        # Implement off-target analysis here
        pass

    def perform_cotargeting(self, selected_targets, endonucleases):
        # Implement cotargeting analysis here
        pass

    def load_gene_viewer_data(self):
        annotation_files = self.get_annotation_files()
        if annotation_files:
            self.annotation_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_files[0])
            self.annotation_parser.set_annotation_file(self.annotation_path)
            self.available_genes = [target['feature_name'] for target in self.targets]
        else:
            self.global_settings.logger.warning("No annotation files found.")
            self.available_genes = []

    def get_gene_data(self, gene_name):
        self._load_gene_data(gene_name)
        return {
            'sequence': self.extended_sequence,
            'info': self.gene_info,
            'start': self.current_gene_start,
            'end': self.current_gene_end
        }

    def update_gene_viewer_indices(self, start, end):
        if start < end and start >= self.current_gene_start and end <= self.current_gene_end:
            relative_start = start - self.current_gene_start + 30  # Add 30 to account for the left placeholder
            relative_end = end - self.current_gene_start + 30
            self.extended_sequence = self.extended_sequence[:30] + self.extended_sequence[relative_start:relative_end] + self.extended_sequence[-30:]
            return True
        return False

    def reset_gene_viewer_indices(self):
        if self.available_genes:
            self._load_gene_data(self.available_genes[0])

    def get_available_genes(self):
        return self.available_genes

    def get_gene_sequence(self):
        return self.gene_sequence

    def get_highlighted_gene_sequence(self):
        return self.highlighted_sequence

    def highlight_targets_in_gene_viewer(self, selected_targets):
        highlighted_sequence = self.extended_sequence
        for target in selected_targets:
            sequence = target['sequence']
            strand = target['strand']
            
            if strand == '+':
                index = highlighted_sequence.upper().find(sequence.upper())
                if index != -1:
                    highlighted_sequence = (
                        highlighted_sequence[:index] +
                        f"<span style='background-color: green;'>{highlighted_sequence[index:index+len(sequence)]}</span>" +
                        highlighted_sequence[index+len(sequence):]
                    )
            else:  # strand == '-'
                rev_comp_sequence = str(Seq(sequence).reverse_complement())
                index = highlighted_sequence.upper().find(rev_comp_sequence.upper())
                if index != -1:
                    highlighted_sequence = (
                        highlighted_sequence[:index] +
                        f"<span style='background-color: red;'>{highlighted_sequence[index:index+len(rev_comp_sequence)]}</span>" +
                        highlighted_sequence[index+len(rev_comp_sequence):]
                    )
        
        self.highlighted_sequence = highlighted_sequence
        return self.highlighted_sequence

    def get_filter_options(self):
        return self.filter_options

    def set_filter_options(self, options):
        self.filter_options = options

    def get_scoring_options(self):
        return self.scoring_options

    def set_scoring_options(self, options):
        self.scoring_options = options

    def get_filtered_targets(self):
        # Implement filtering logic based on self.filter_options
        return self.targets

    def export_targets(self, selected_targets, file_path):
        # Implement export logic here
        pass

    def get_sequence_from_annotation(self, chrom, start, end):
        for record in SeqIO.parse(self.annotation_path, 'genbank'):
            if record.id == chrom:
                return str(record.seq[start:end])
        return ""
