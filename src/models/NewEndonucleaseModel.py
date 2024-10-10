import os
from utils.ui import show_error

class NewEndonucleaseModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.app_dir = global_settings.get_app_dir()

    def get_on_off_data(self):
        try:
            filename = self.global_settings.get_casper_info_path()
            retList_on = []
            retList_off = []
            with open(filename, 'r') as f:
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

    def write_new_endonuclease(self, new_endonuclease_str):
        try:
            casper_info_path = os.path.join(self.app_dir, 'CASPERinfo')
            new_file_path = os.path.join(self.app_dir, "new_file")
            with open(casper_info_path, 'r') as f, open(new_file_path, 'w+') as f1:
                for line in f:
                    f1.write(line)
                    if 'ENDONUCLEASES' in line:
                        f1.write(new_endonuclease_str + '\n')
            os.remove(casper_info_path)
            os.rename(new_file_path, casper_info_path)
        except Exception as e:
            show_error("Error in write_new_endonuclease() in New Endonuclease Model.", str(e))

    def is_duplicate_abbreviation(self, abbr):
        for key in self.global_settings.main_window.organisms_to_endos:
            endo = self.global_settings.main_window.organisms_to_endos[key]
            if abbr in endo:
                return True
        return False

    def create_endonuclease_string(self, form_data):
        pam = form_data['pam']
        if len(pam.split(',')) > 0:
            pam = [x.strip() for x in pam.split(',')]
            pam = ",".join(pam)

        argument_list = [
            form_data['abbr'],
            pam,
            form_data['five_len'],
            form_data['seed_len'],
            form_data['three_len'],
            form_data['pam_dir'],
            form_data['name'],
            form_data['crisprtype'],
            form_data['on_scoring'],
            form_data['off_scoring']
        ]

        return ";".join(str(arg) for arg in argument_list)
