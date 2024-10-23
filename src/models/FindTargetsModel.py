from models.HomeWindowModel import HomeWindowModel
from models.CSPRparser import CSPRparser
from models.AnnotationParser import AnnotationParser
import os

class FindTargetsModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.results = {}
        self.annotation_parser = AnnotationParser(global_settings)

    def find_targets(self, input_data):
        self.global_settings.logger.debug(f"Received input data: {input_data}")
        print(f"Received input data: {input_data}")
        
        organism = input_data['organism']
        endo = input_data['endonuclease']
        
        org_files = self.get_organism_to_files()

        if organism not in org_files:
            error_msg = f"Organism '{organism}' not found in the database. Available organisms: {list(org_files.keys())}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)
        
        if endo not in org_files[organism]:
            error_msg = f"Endonuclease '{endo}' not found for organism '{organism}'. Available endonucleases: {list(org_files[organism].keys())}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)
        
        file_path = os.path.join(self.global_settings.get_db_path(), org_files[organism][endo][0])
        
        parser = CSPRparser(file_path, self.global_settings.get_casper_info_path())
        
        if input_data['search_type'] == 'feature':
            self.results = self.find_targets_by_feature(parser, input_data)
        elif input_data['search_type'] == 'position':
            self.results = self.find_targets_by_position(parser, input_data)
        elif input_data['search_type'] == 'sequence':
            self.results = self.find_targets_by_sequence(parser, input_data)
        else:
            error_msg = f"Invalid search type: {input_data['search_type']}"
            self.global_settings.logger.error(error_msg)
            raise ValueError(error_msg)
        
        return self.results

    def find_targets_by_feature(self, parser, input_data):
        annotation_file = input_data['annotation_file']
        search_query = input_data['search_query'].strip().lower()
        
        annotation_file_path = os.path.join(self.global_settings.get_db_path(), 'GBFF', annotation_file)
        self.annotation_parser.set_annotation_file(annotation_file_path)
        
        try:
            results_list = self.annotation_parser.genbank_search([search_query])
            self.global_settings.logger.debug(f"Genbank search results: {results_list}")
        except Exception as e:
            self.global_settings.logger.error(f"Error in genbank_search: {str(e)}")
            return []
        
        formatted_results = []
        for chrom, feature in results_list:
            if feature.type in ['CDS']:
                feature_info = self._get_feature_info(feature)
                if (search_query == feature_info['feature_name'].lower() or 
                    search_query == feature_info['feature_id'].lower() or
                    search_query in feature_info['feature_name'].lower() or
                    search_query in feature_info['feature_id'].lower() or
                    search_query in feature_info['feature_description'].lower()):
                    formatted_results.append({
                        'feature_type': feature.type,
                        'chromosome': chrom,
                        'feature_id': feature_info['feature_id'],
                        'feature_name': feature_info['feature_name'],
                        'feature_description': feature_info['feature_description'],
                        'location': f"{feature.location.start}-{feature.location.end}",
                        'strand': '+' if feature.strand == 1 else '-'
                    })
        
        self.global_settings.logger.debug(f"Total features found: {len(formatted_results)}")
        return formatted_results

    def _get_feature_info(self, feature):
        return {
            'feature_id': self._get_feature_id(feature),
            'feature_name': self._get_feature_name(feature),
            'feature_description': self._get_feature_description(feature)
        }

    def _find_associated_cds(self, gene_feature):
        for feature in gene_feature.parent.features:
            if feature.type == 'CDS' and feature.location == gene_feature.location:
                return feature
        return None

    def _get_feature_id(self, feature):
        for key in ['locus_tag', 'protein_id', 'id']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def _get_feature_name(self, feature):
        for key in ['gene', 'product']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def _get_feature_description(self, feature):
        for key in ['product', 'note']:
            if key in feature.qualifiers:
                return feature.qualifiers[key][0]
        return "N/A"

    def find_targets_by_position(self, parser, input_data):
        search_query = input_data['search_query']
        chrom, start, end = map(int, search_query.split(','))
        pos_tuple = (chrom, start - 1, end)  # Adjust for 0-based indexing
        
        targets = parser.read_targets(f"position_{chrom}_{start}_{end}", pos_tuple, input_data['endonuclease'])
        return self._format_results(targets)

    def find_targets_by_sequence(self, parser, input_data):
        search_query = input_data['search_query'].upper()
        annotation_file = input_data['annotation_file']
        
        self.annotation_parser.annotationFileName = os.path.join(self.global_settings.get_db_path(), annotation_file)
        sequence_info = self.annotation_parser.get_sequence_info(search_query)
        
        if sequence_info:
            chrom, start, end = sequence_info
            pos_tuple = (chrom, start - 1, end)
            targets = parser.read_targets(f"sequence_{start}_{end}", pos_tuple, input_data['endonuclease'])
            return self._format_results(targets)
        else:
            return []

    def _format_results(self, targets):
        formatted_results = []
        for target in targets:
            formatted_results.append({
                'location': target[0],
                'sequence': target[1],
                'pam': target[2],
                'score': target[3],
                'strand': target[4],
                'endonuclease': target[5]
            })
        return formatted_results
