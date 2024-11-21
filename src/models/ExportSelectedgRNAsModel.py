import os
from typing import Dict, List
from utils.sequence_utils import get_table_headers

class ExportSelectedgRNAsModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.data = {
            'selected_items': [],
            'window_type': "", 
            'has_locus_tag': False,
            'has_gene_name': False
        }

    def set_export_data(self, items: List, window_type: str) -> None:
        """Set the selected items and window type"""
        try:
            self.data['selected_items'] = items
            self.data['window_type'] = window_type
            self.logger.debug(f"Set export data: {len(items)} items, window type: {window_type}")
        except Exception as e:
            self.logger.error(f"Error setting export data: {str(e)}")
            raise

    def get_headers(self) -> List[str]:
        """Get appropriate headers based on window type"""
        try:
            if self.data['window_type'] == "Multitargeting":
                return self._get_multitargeting_headers()
            elif self.data['window_type'] == "Population Analysis":
                return self._get_population_analysis_headers()
            else:
                return self._get_view_targets_headers()
        except Exception as e:
            self.logger.error(f"Error getting headers: {str(e)}")
            raise

    def _get_multitargeting_headers(self) -> List[str]:
        """Get headers for multitargeting window"""
        try:
            headers = [
                "Seed", "Total Repeats", "Avg. Repeats/Scaffold",
                "Consensus Sequence", "% Consensus", "Score", "PAM", "Strand"
            ]
            
            insertion_index = headers.index("Consensus Sequence")
            headers.insert(insertion_index + 1, "Full Sequence")
            
            self.logger.debug(f"Multitargeting headers: {headers}")
            return headers
            
        except Exception as e:
            self.logger.error(f"Error getting multitargeting headers: {str(e)}")
            raise

    def _get_population_analysis_headers(self) -> List[str]:
        """Get headers for population analysis window"""
        try:
            headers = [
                "Seed", "% Coverage", "Total Repeats", "Avg. Repeats/Scaffold",
                "Consensus Sequence", "% Consensus", "Score", "PAM", "Strand"
            ]
            
            insertion_index = headers.index("Consensus Sequence")
            headers.insert(insertion_index + 1, "Full Sequence")
            
            self.logger.debug(f"Population Analysis headers: {headers}")
            return headers
            
        except Exception as e:
            self.logger.error(f"Error getting population analysis headers: {str(e)}")
            raise

    def _get_view_targets_headers(self) -> List[str]:
        """Get headers for view targets window"""
        try:
            headers = [
                "Location", "Endonuclease", "Sequence", "Strand", 
                "PAM", "Score", "Off-Target" 
            ]
            
            sequence_index = headers.index("Sequence")
            headers.insert(sequence_index + 1, "Full Sequence")
            
            if self.data['selected_items']:
                first_target = self.data['selected_items'][0]
                
                if 'locus_tag' in first_target or 'feature_id' in first_target:
                    headers.append("Locus_Tag")
                    self.data['has_locus_tag'] = True
                    
                if 'gene_name' in first_target or 'feature_name' in first_target:
                    headers.append("Gene_Name")
                    self.data['has_gene_name'] = True
                    
            return headers
            
        except Exception as e:
            self.logger.error(f"Error getting view targets headers: {str(e)}")
            raise

    def get_file_extension(self, delimiter: str) -> str:
        """Get appropriate file extension based on delimiter"""
        if delimiter == ",":
            return ".csv"
        elif delimiter == r"\t":
            return ".tsv"
        return ".txt"

    def get_full_path(self, directory: str, filename: str, delimiter: str) -> str:
        """Construct full file path"""
        if '.' in filename:
            return os.path.join(directory, filename)
        extension = self.get_file_extension(delimiter)
        return os.path.join(directory, filename + extension)
