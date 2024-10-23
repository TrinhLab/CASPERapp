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
        try:
            new_file_path = os.path.join(self.app_dir_path, "new_file")
            with open(self.casper_info_path, 'r') as f, open(new_file_path, 'w+') as f1:
                for line in f:
                    f1.write(line)
                    if 'ENDONUCLEASES' in line:
                        f1.write(new_endonuclease_str + '\n')
            os.remove(self.casper_info_path)
            os.rename(new_file_path, self.casper_info_path)
            self.global_settings.config_manager.load_endonucleases_data()
            self.endonuclease_updated.emit()
        except Exception as e:
            show_error(self.global_settings, "Error in create_new_endonuclease() in New Endonuclease Model.", str(e))

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
        pam = form_data['endonuclease_pam_sequence']
        if ',' in pam:
            pam = ','.join([x.strip() for x in pam.split(',')])

        argument_list = [
            form_data['endonuclease_organism'],
            form_data['endonuclease_abbreviation'],
            form_data['endonuclease_CRISPR_type'],
            pam,
            form_data['endonuclease_five_prime_length'],
            form_data['endonuclease_seed_length'],
            form_data['endonuclease_three_prime_length'],
            form_data['endonuclease_direction'],
            form_data['endonuclease_on_target_scoring'],
            form_data['endonuclease_off_target_scoring']
        ]

        return ";".join(str(arg) for arg in argument_list)

    def update_endonuclease(self, selected, form_data):
        try:
            new_endonuclease_str = self.create_endonuclease_string(form_data)
            with open(self.casper_info_path, 'r') as f:
                lines = f.readlines()
            
            updated = False
            selected_abbr = selected.split(' - ')[0]
            for i, line in enumerate(lines):
                if line.startswith('ENDONUCLEASES'):
                    continue
                fields = line.strip().split(';')
                if len(fields) >= 2 and fields[1] == selected_abbr:
                    lines[i] = new_endonuclease_str + '\n'
                    updated = True
                    break
            
            if updated:
                with open(self.casper_info_path, 'w') as f:
                    f.writelines(lines)
                self.global_settings.config_manager.load_endonucleases_data()
                self.endonuclease_updated.emit()
            else:
                raise ValueError(f"Endonuclease '{selected}' not found in CASPERinfo file")
        except Exception as e:
            show_error(self.global_settings, "Error updating endonuclease", str(e))

    def delete_endonuclease(self, selected):
        try:
            with open(self.casper_info_path, 'r') as f:
                lines = f.readlines()
            
            deleted = False
            new_lines = []
            selected_abbr = selected.split(' - ')[0]
            for line in lines:
                if line.startswith('ENDONUCLEASES'):
                    new_lines.append(line)
                    continue
                fields = line.strip().split(';')
                if len(fields) >= 2 and fields[1] == selected_abbr:
                    deleted = True
                    continue  # Skip this line to delete it
                new_lines.append(line)  # Keep all other lines
            
            if deleted:
                with open(self.casper_info_path, 'w') as f:
                    f.writelines(new_lines)
                self.global_settings.config_manager.load_endonucleases_data()
                self.endonuclease_updated.emit()
            else:
                raise ValueError(f"Endonuclease '{selected}' not found in CASPERinfo file")
        except Exception as e:
            show_error(self.global_settings, "Error deleting endonuclease", str(e))
