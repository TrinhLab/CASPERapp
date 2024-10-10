import os
import platform

class NCBIRenameWindowModel:
    def __init__(self, settings):
        self.settings = settings
        self.logger = settings.get_logger()
        self.files = []

    def set_files(self, files):
        self.files = files

    def get_files(self):
        return self.files

    def rename_file(self, old_name, new_name):
        try:
            file_type = 'GBFF' if old_name.endswith('.gbff') else 'FNA'
            old_path = os.path.join(self.settings.CSPR_DB, file_type, os.path.basename(old_name))
            new_path = os.path.join(self.settings.CSPR_DB, file_type, new_name)

            if platform.system() == "Windows":
                new_path = new_path.replace("/", "\\")
            else:
                new_path = new_path.replace(" ", "")  # Remove spaces for Unix-like systems

            if not os.path.isfile(new_path):
                os.rename(old_path, new_path)
                return True, ""
            else:
                return False, f"The filename: {new_name} already exists. Please use a different name."
        except Exception as e:
            self.logger.error(f"Error renaming file {old_name} to {new_name}: {str(e)}")
            return False, f"An error occurred while renaming the file: {str(e)}"

    def validate_new_name(self, old_name, new_name):
        if not new_name:
            return False, "New filename cannot be empty."
        
        if old_name.endswith('.gbff'):
            if not new_name.endswith('.gbff'):
                new_name += '.gbff'
        elif old_name.endswith('.fna'):
            if not new_name.endswith('.fna'):
                new_name += '.fna'
        
        # Check for invalid characters in the filename
        invalid_chars = r'<>:"/\|?*'
        if any(char in new_name for char in invalid_chars):
            return False, f"Filename contains invalid characters. Avoid using: {invalid_chars}"
        
        return True, new_name

    def process_rename_batch(self, new_names):
        results = []
        for old_name, new_name in zip(self.files, new_names):
            if new_name:
                is_valid, message_or_name = self.validate_new_name(old_name, new_name)
                if is_valid:
                    success, message = self.rename_file(old_name, message_or_name)
                    results.append((success, old_name, message_or_name, message))
                else:
                    results.append((False, old_name, new_name, message_or_name))
            else:
                results.append((True, old_name, old_name, "No change"))
        
        return results

    def clear_files(self):
        self.files.clear()
