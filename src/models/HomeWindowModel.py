import os
import glob
from typing import Dict, List
from utils.ui import show_error

class HomeWindowModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.data = {}
        self.settings = {}
        self.dbpath = ""
        self.inputstring = ""
        self.anno_name = ""
        self.endo_name = ""
        self.org = ""
        self.genlib_list = []
        self.results_list = []

    def load_data(self) -> None:
        try:
            self.load_organisms_and_endonuclease()
            self.load_annotation_files()
        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")
            show_error(self.global_settings, "Error Loading Data", str(e))

    def load_organisms_and_endonuclease(self) -> None:
        try:
            self.data["organism_to_files"] = {}
            self.data["organism_to_endonuclease"] = {}
            cspr_files = glob.glob(os.path.join(self.global_settings.get_db_path(), "*.cspr"))

            for file in cspr_files:
                file_name = os.path.basename(file)
                file_name_no_ext = file_name[:-5]
                endonuclease = file_name_no_ext[file_name_no_ext.rfind("_")+1:]

                with open(file, 'r') as f:
                    organism = f.readline().strip().replace("GENOME: ", '')
                
                if organism in self.data["organism_to_files"]:
                    self.data["organism_to_files"][organism][endonuclease] = [file_name, file_name.replace(".cspr", "_repeats.db")]
                else:
                    self.data["organism_to_files"][organism] = {endonuclease: [file_name, file_name.replace(".cspr", "_repeats.db")]}
                
                if organism in self.data["organism_to_endonuclease"]:
                    self.data["organism_to_endonuclease"][organism].append(endonuclease)
                else:
                    self.data["organism_to_endonuclease"][organism] = [endonuclease]
            
            self.logger.debug("Successfully loaded genome and endonuclease data.")
        except Exception as e:
            self.logger.error(f"Error in load_organisms_and_endonuclease: {str(e)}")

    def load_annotation_files(self) -> None:
        try:
            self.data["annotation_files"] = glob.glob(os.path.join(self.global_settings.CSPR_DB, "**", "*.gb*"), recursive=True)
            self.data["annotation_files"] = [os.path.basename(file) for file in self.data["annotation_files"]]
            self.data["annotation_files"].sort(key=str.lower)
            self.data["annotation_files"].append("None")
            self.logger.debug("Successfully loaded annotation files.")
        except Exception as e:
            self.logger.error(f"Error in load_annotation_files: {str(e)}")

    def get_organism_to_files(self) -> Dict[str, Dict[str, List[str]]]:
        return self.data.get("organism_to_files", {})

    def get_organism_to_endonuclease(self) -> Dict[str, List[str]]:
        return self.data.get("organism_to_endonuclease", {})

    def get_annotation_files(self) -> List[str]:
        return self.data.get("annotation_files", [])

    # Add other methods that handle data processing, validation, and storage
