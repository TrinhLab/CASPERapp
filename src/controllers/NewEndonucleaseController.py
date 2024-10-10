import os
from PyQt6 import QtWidgets
from utils.ui import show_message, show_error
from views.NewEndonucleaseView import NewEndonucleaseView
from models.NewEndonucleaseModel import NewEndonucleaseModel

class NewEndonucleaseController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        try:
            self.model = NewEndonucleaseModel(global_settings)
            self.view = NewEndonucleaseView(global_settings)
        
            self._connect_ui_signals()
            self._init_ui()
        except Exception as e:
            show_error(self.global_settings, "Error initializing NewEndonucleaseController", str(e))
            raise

    def _connect_ui_signals(self):
        try:
            # boxlayhbotButtons
            self.view.push_button_submit.clicked.connect(self._handle_submit)
            self.view.push_button_reset_form.clicked.connect(self._handle_reset_form)
            print(f"Push Button reset form found", self.view.push_button_reset_form)
            self.view.show()
        except Exception as e:
            show_error(self.global_settings, "Error setting up connections in NewEndonucleaseController", str(e))

    
    def _init_ui(self):
        try:
            self.load_combo_box_data()
        except Exception as e:
            show_error(self.global_settings, "Error initializing UI in NewEndonucleaseController", str(e))

    def load_combo_box_data(self):
        on_list, off_list = self.model.get_on_off_data()
        self.view.update_combo_box_on_target_matrix(on_list)
        self.view.update_combo_box_off_target_matrix(off_list)

    def _handle_submit(self):
        try:
            form_data = self.view.get_form_data()
            
            if self.validate_form_data(form_data):
                new_endonuclease_str = self.model.create_endonuclease_string(form_data)
                self.model.write_new_endonuclease(new_endonuclease_str)
                
                # Refresh endonuclease dropdown in New Genome
                self.global_settings.get_new_genome().fillEndo()
                
                self.view.clear_form()
            
        except Exception as e:
            show_error(self.global_settings, "Error in _handle_submit() in New Endonuclease Controller", str(e))

    def _handle_reset_form(self):
        try:
            print("Button pressed")
            self.view.line_edit_organism.clear()
            self.view.line_edit_abbreviation.clear()
            self.view.line_edit_CRISPR_type.clear()
            self.view.line_edit_seed_length.clear()
            self.view.line_edit_five_length.clear()
            self.view.line_edit_three_length.clear()
            self.view.line_edit_pam_sequence.clear()
        except Exception as e:
            show_error(self.global_settings, "Error in _handle_reset_form() in New Endonuclease Controller", str(e))

    def validate_form_data(self, form_data):
        # Check for empty fields
        for key, value in form_data.items():
            if value == "":
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Empty Field",
                    message="Please fill in all fields."
                )
                return False

        # Check for semicolons
        if any(';' in str(value) for value in form_data.values()):
            show_message(
                fontSize=12,
                icon=QtWidgets.QMessageBox.Icon.Critical,
                title="Invalid Semicolon",
                message="Invalid character used: ';'."
            )
            return False

        # Check PAM sequence
        valid_pam_chars = set('ACDEFGHIKLMNPQRSTVWY')
        if not set(form_data['pam']).issubset(valid_pam_chars):
            show_message(
                fontSize=12,
                icon=QtWidgets.QMessageBox.Icon.Critical,
                title="Invalid PAM",
                message="Invalid characters in PAM Sequence."
            )
            return False

        # Check for duplicate endo abbreviations
        if self.model.is_duplicate_abbreviation(form_data['abbr']):
            show_message(
                fontSize=12,
                icon=QtWidgets.QMessageBox.Icon.Critical,
                title="Duplicate endo name.",
                message="The given abbreviation already exists. Please choose a unique identifier."
            )
            return False

        return True
