import os
import glob
from typing import Dict, List
from utils.ui import show_error
from models.DatabaseManager import FileChangeType

class HomeWindowModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.data = {
            'organism_to_files': {},
            'organism_to_endonuclease': {},
            'annotation_files': set()  # Using set for efficient updates
        }
        self.load_data()

    def load_data(self) -> None:
        """Load all required data"""
        try:
            self.load_organisms_and_endonuclease()
            self.load_annotation_files()
            self.logger.debug("Successfully loaded all home window data")
        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")
            show_error(self.global_settings, "Error Loading Data", str(e))

    def update_for_file_changes(self, changes: Dict[FileChangeType, List[str]]) -> None:
        """Update model data based on file changes"""
        try:
            needs_organism_reload = False
            needs_annotation_reload = False

            # Check what types of files changed
            if (FileChangeType.CSPR_ADDED in changes or 
                FileChangeType.CSPR_REMOVED in changes):
                needs_organism_reload = True
                
            if (FileChangeType.GBFF_ADDED in changes or 
                FileChangeType.GBFF_REMOVED in changes):
                needs_annotation_reload = True

            # Reload only what's necessary
            if needs_organism_reload:
                self.load_organisms_and_endonuclease()
                
            if needs_annotation_reload:
                self.load_annotation_files()
                
            self.logger.debug(f"Updated model data for changes: {changes}")
            
        except Exception as e:
            self.logger.error(f"Error updating for file changes: {str(e)}")

    def load_organisms_and_endonuclease(self) -> None:
        """Load organism and endonuclease data from CSPR files"""
        try:
            # Clear existing data
            self.data["organism_to_files"] = {}
            self.data["organism_to_endonuclease"] = {}
            
            cspr_files = glob.glob(os.path.join(self.global_settings.get_db_path(), "*.cspr"))
            
            for file in cspr_files:
                file_name = os.path.basename(file)
                file_name_no_ext = file_name[:-5]
                endonuclease = file_name_no_ext[file_name_no_ext.rfind("_")+1:]

                with open(file, 'r') as f:
                    organism = f.readline().strip().replace("GENOME: ", '')
                
                # Update organism to files mapping
                if organism not in self.data["organism_to_files"]:
                    self.data["organism_to_files"][organism] = {}
                self.data["organism_to_files"][organism][endonuclease] = [
                    file_name, 
                    file_name.replace(".cspr", "_repeats.db")
                ]
                
                # Update organism to endonuclease mapping
                if organism not in self.data["organism_to_endonuclease"]:
                    self.data["organism_to_endonuclease"][organism] = []
                if endonuclease not in self.data["organism_to_endonuclease"][organism]:
                    self.data["organism_to_endonuclease"][organism].append(endonuclease)
            
            self.logger.debug(f"Loaded data for {len(self.data['organism_to_files'])} organisms")
            
        except Exception as e:
            self.logger.error(f"Error loading organisms and endonucleases: {str(e)}")
            raise

    def load_annotation_files(self) -> None:
        """Load annotation files from the database directory"""
        try:
            # Get all .gb* files recursively
            annotation_files = glob.glob(
                os.path.join(self.global_settings.get_db_path(), "**", "*.gb*"), 
                recursive=True
            )
            
            # Process files
            self.data["annotation_files"] = {
                os.path.basename(file) for file in annotation_files 
                if not file.endswith('.index')  # Exclude index files
            }
            
            self.logger.debug(f"Loaded {len(self.data['annotation_files'])} annotation files")
            
        except Exception as e:
            self.logger.error(f"Error loading annotation files: {str(e)}")
            raise

    def get_organism_to_files(self) -> Dict[str, Dict[str, List[str]]]:
        """Get mapping of organisms to their files"""
        return self.data['organism_to_files']

    def get_organism_to_endonuclease(self) -> Dict[str, List[str]]:
        """Get mapping of organisms to their endonucleases"""
        return self.data.get("organism_to_endonuclease", {})

    def get_annotation_files(self) -> List[str]:
        """Get list of annotation files"""
        return sorted(self.data.get("annotation_files", set()), key=str.lower)
    
    def find_targets(self, input_data: dict) -> None:
        pass

    # Add other methods that handle data processing, validation, and storage
