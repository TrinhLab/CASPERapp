from PyQt6.QtWidgets import QMessageBox
from Bio import SeqIO
import os
import traceback

class AnnotationParser:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.annotation_file_name = ""
        self.available_genes = []

    def set_annotation_file(self, file_path):
        self.annotation_file_name = file_path
        self.logger.debug(f"Set annotation file to: {file_path}")
        self._parse_available_genes()

    def get_max_chrom(self):
        try:
            parser = SeqIO.parse(self.annotation_file_name, 'genbank')
            max_chrom = sum(1 for _ in parser)
            return max_chrom
        except Exception as e:
            self.logger.error(f"Error in get_max_chrom: {str(e)}")
            self._show_error("Error in get_max_chrom", str(e))
            return 0

    def get_sequence_info(self, query):
        # Implement this method if needed
        pass

    def genbank_search(self, queries):
        try:
            if not self.annotation_file_name:
                raise ValueError("Annotation file not set")
            
            self.logger.debug(f"Searching in annotation file: {self.annotation_file_name}")
            results_list = []
            
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                for feature in record.features:
                    if feature.type in ['CDS', 'gene']:
                        feature_info = self._get_feature_info(feature)
                        for query in queries:
                            if query.lower() in feature_info['feature_name'].lower() or \
                               query.lower() in feature_info['feature_id'].lower() or \
                               query.lower() in feature_info['feature_description'].lower():
                                results_list.append((record.id, feature))
            
            self.logger.debug(f"Found {len(results_list)} results")
            return results_list
        except Exception as e:
            self.logger.error(f"Error in genbank_search: {str(e)}")
            raise

    def find_which_file_version(self):
        try:
            if not self.annotation_file_name or os.path.basename(self.annotation_file_name) == "None":
                return -1
            if self.annotation_file_name.endswith(('.gbff', '.gbk')):
                return "gbff"
            else:
                return -1
        except Exception as e:
            self.logger.error(f"Error in find_which_file_version: {str(e)}")
            self._show_error("Error in find_which_file_version", str(e))
            return -1

    def _show_error(self, title, message):
        QMessageBox.critical(None, title, f"{message}\n\nFor more information, check the log file.")

    @staticmethod
    def flatten_list(t):
        return [item.lower() for sublist in t for item in sublist]

    def _get_feature_info(self, feature):
        return {
            'feature_id': self._get_feature_id(feature),
            'feature_name': self._get_feature_name(feature),
            'feature_location': self._get_feature_location(feature),
            'feature_description': self._get_feature_description(feature)
        }

    def _get_feature_id(self, feature):
        for key in ['locus_tag']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def _get_feature_name(self, feature):
        for key in ['gene']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"
    
    def _get_feature_location(self, feature):
        if feature.location:
            start = feature.location.start
            end = feature.location.end
            strand = '+' if feature.location.strand == 1 else '-'
            return f"{start}:{end}({strand})"
        return "N/A"

    def _get_feature_description(self, feature):
        for key in ['product', 'note']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def get_available_genes(self):
        return self.available_genes

    def get_gene_data(self, gene_name):
        try:
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                for feature in record.features:
                    if feature.type == 'gene' and self._get_feature_name(feature) == gene_name:
                        self.logger.debug(f"Found gene {gene_name} in annotation file")
                        self.logger.debug(f"Its sequence is {feature.extract(record.seq)}")
                        self.logger.debug(f"Its feature is {feature}")
                        self.logger.debug(f"Its feature info is {self._get_feature_info(feature)}")
                        return {
                            'sequence': str(feature.extract(record.seq)),
                            'info': self._get_feature_info(feature)
                        }
            self.logger.warning(f"Gene {gene_name} not found in annotation file")
            return None
        except Exception as e:
            self.logger.error(f"Error in get_gene_data: {str(e)}")
            return None

    def _parse_available_genes(self):
        self.available_genes = []
        try:
            for record in SeqIO.parse(self.annotation_file_name, "genbank"):
                for feature in record.features:
                    if feature.type == 'gene':
                        self.available_genes.append(self._get_feature_name(feature))
        except Exception as e:
            self.logger.error(f"Error parsing available genes: {str(e)}")

    def get_full_gene_sequence(self):
        # Implement this method if needed
        pass
