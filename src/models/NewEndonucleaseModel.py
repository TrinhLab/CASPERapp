import os
from utils.ui import show_error
from PyQt6.QtCore import QObject, pyqtSignal

class NewEndonucleaseModel(QObject):
    endonuclease_updated = pyqtSignal()

    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.app_dir_path = global_settings.get_app_dir_path()
        self.casper_info_path = global_settings.get_casper_info_path()

    def get_on_off_data(self):
        try:
            retList_on = []
            retList_off = []
            with open(self.casper_info_path, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    line = str(line)
                    if "ON-TARGET DATA" in line:
                        index = i
                        while index < len(lines) and "-----" not in lines[index]:
                            if "DATA:" in lines[index]:
                                retList_on.append(lines[index].split("DATA:")[-1].strip())
                            index += 1
                    elif "OFF-TARGET MATRICES" in line:
                        index = i
                        while index < len(lines) and "-----" not in lines[index]:
                            if "MATRIX:" in lines[index]:
                                retList_off.append(lines[index].split("MATRIX:")[-1].strip())
                            index += 1
            return retList_on, retList_off
        except Exception as e:
            show_error(self.global_settings, "Error in get_on_off_data() in New Endonuclease Model", str(e))
            return [], []  # Return empty lists if there's an error

    def create_new_endonuclease(self, new_endonuclease_str):
        """Add a new endonuclease to CASPERinfo file"""
        try:
            found_section = False
            new_lines = []
            
            with open(self.casper_info_path, 'r') as f:
                for line in f:
                    new_lines.append(line)
                    if line.strip() == 'ENDONUCLEASES':
                        found_section = True
                        new_lines.append(new_endonuclease_str + '\n')
            
            if not found_section:
                new_lines.append('ENDONUCLEASES\n')
                new_lines.append(new_endonuclease_str + '\n')
                
            with open(self.casper_info_path, 'w') as f:
                f.writelines(new_lines)
                
            self.global_settings.config_manager.load_endonucleases_data()
            self.endonuclease_updated.emit()
            
        except Exception as e:
            show_error(self.global_settings, "Error in create_new_endonuclease()", str(e))

    def is_duplicate_abbreviation(self, abbr):
        try:
            # Get the organism_to_endonuclease data from GlobalSettings
            organism_to_endonuclease = self.global_settings.get_endonucleases()
            
            # Check if the abbreviation exists in any of the endonucleases
            for endonucleases in organism_to_endonuclease.values():
                if abbr == endonucleases.get('endonuclease_abbreviation'):
                    return True
            return False
        except Exception as e:
            show_error(self.global_settings, "Error checking for duplicate abbreviation", str(e))
            return False

    def create_endonuclease_string(self, form_data):
        """Create a properly formatted endonuclease string for CASPERinfo"""
        # Format: abbr;PAM;5_prime_length;seed_length;3_prime_length;direction;organism;CRISPR_type;on_target;off_target
        argument_list = [
            form_data['endonuclease_abbreviation'],
            form_data['endonuclease_pam_sequence'],
            form_data['endonuclease_five_prime_length'],
            form_data['endonuclease_seed_length'],
            form_data['endonuclease_three_prime_length'],
            '3' if form_data['endonuclease_direction'] == '3' else '5',
            form_data['endonuclease_organism'],
            form_data['endonuclease_CRISPR_type'],
            form_data['endonuclease_on_target_scoring'],
            form_data['endonuclease_off_target_scoring']
        ]
        return ";".join(str(arg) for arg in argument_list)

    def update_endonuclease(self, selected, form_data):
        """Update an existing endonuclease in CASPERinfo file"""
        try:
            new_endonuclease_str = self.create_endonuclease_string(form_data)
            new_lines = []
            updated = False
            selected_abbr = selected.split(' - ')[0]  # Get abbreviation part
            
            with open(self.casper_info_path, 'r') as f:
                in_endo_section = False
                for line in f:
                    if line.strip() == 'ENDONUCLEASES':
                        in_endo_section = True
                        new_lines.append(line)
                        continue
                        
                    if in_endo_section and line.strip():
                        fields = line.strip().split(';')
                        if fields[0] == selected_abbr:
                            new_lines.append(new_endonuclease_str + '\n')
                            updated = True
                        else:
                            new_lines.append(line)
                    else:
                        new_lines.append(line)
            
            if not updated:
                raise ValueError(f"Endonuclease '{selected}' not found")
                
            with open(self.casper_info_path, 'w') as f:
                f.writelines(new_lines)
                
            self.global_settings.config_manager.load_endonucleases_data()
            self.endonuclease_updated.emit()
            
        except Exception as e:
            show_error(self.global_settings, "Error updating endonuclease", str(e))

    def delete_endonuclease(self, selected):
        """Delete an endonuclease from CASPERinfo file"""
        try:
            new_lines = []
            deleted = False
            selected_abbr = selected.split(' - ')[0]
            
            with open(self.casper_info_path, 'r') as f:
                in_endo_section = False
                for line in f:
                    if line.strip() == 'ENDONUCLEASES':
                        in_endo_section = True
                        new_lines.append(line)
                        continue
                        
                    if in_endo_section and line.strip():
                        fields = line.strip().split(';')
                        if fields[0] == selected_abbr:
                            deleted = True
                            continue
                    new_lines.append(line)
            
            if not deleted:
                raise ValueError(f"Endonuclease '{selected}' not found")
                
            with open(self.casper_info_path, 'w') as f:
                f.writelines(new_lines)
                
            self.global_settings.config_manager.load_endonucleases_data()
            self.endonuclease_updated.emit()
            
        except Exception as e:
            show_error(self.global_settings, "Error deleting endonuclease", str(e))
