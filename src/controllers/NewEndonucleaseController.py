import os
from PyQt6 import QtWidgets
from utils.ui import show_message, show_error
from views.NewEndonucleaseView import NewEndonucleaseView
from models.NewEndonucleaseModel import NewEndonucleaseModel

class NewEndonucleaseController:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        try:
            self.model = NewEndonucleaseModel(self.settings)
            self.view = NewEndonucleaseView(self.settings)
            self.model.endonuclease_updated.connect(self._on_endonuclease_updated)
        
            self._setup_connections()
            self._init_ui()
        except Exception as e:
            show_error(self.settings, "Error initializing NewEndonucleaseController", str(e))
            raise

    def _setup_connections(self):
        try:
            self.view.push_button_submit.clicked.connect(self._handle_submit)
            self.view.push_button_reset_form.clicked.connect(self._handle_reset_form)
            self.view.push_button_delete.clicked.connect(self._handle_delete)
            self.view.combo_box_select_endonuclease.currentIndexChanged.connect(self._handle_endonuclease_selection)
            self.view.show()
        except Exception as e:
            show_error(self.settings, "Error setting up connections in NewEndonucleaseController", str(e))

    def _init_ui(self):
        try:
            self.load_combo_box_data()
            self._populate_endonuclease_combo_box()
        except Exception as e:
            show_error(self.settings, "Error initializing UI in NewEndonucleaseController", str(e))

    def _populate_endonuclease_combo_box(self):
        endonucleases = list(self.settings.get_endonucleases().keys())
        self.view.update_combo_box_select_endonuclease(endonucleases)

    def _handle_endonuclease_selection(self):
        selected = self.view.get_selected_endonuclease()
        if selected == "":
            self.view.disable_form_elements()
            self.view.clear_form()
            self.view.disable_delete_button()
        elif selected == "Define New Endonuclease":
            self.view.enable_form_elements()
            self.view.clear_form()
            self.view.disable_delete_button()
        elif selected:
            self.view.enable_form_elements()
            endonuclease_data = self.settings.get_endonucleases()[selected]
            self.view.populate_form(endonuclease_data)
            self.view.enable_delete_button()

    def load_combo_box_data(self):
        on_list, off_list = self.model.get_on_off_data()
        self.view.update_combo_box_on_target_matrix(on_list)
        self.view.update_combo_box_off_target_matrix(off_list)

    def _handle_submit(self):
        try:
            form_data = self.view.get_form_data()
            
            selected = self.view.get_selected_endonuclease()
            is_new_endonuclease = selected == "Define New Endonuclease"

            if not self.validate_form_data(form_data, is_new_endonuclease):
                return

            if not is_new_endonuclease:
                # Update existing endonuclease
                self.model.update_endonuclease(selected, form_data)
                show_message("Success", "Endonuclease updated successfully.")
            else:
                # Add new endonuclease
                new_endonuclease_str = self.model.create_endonuclease_string(form_data)
                self.model.create_new_endonuclease(new_endonuclease_str)
                show_message("Success", "New endonuclease added successfully.")
            
            self._populate_endonuclease_combo_box()
            self.view.clear_form()
        
        except Exception as e:
            show_error(self.settings, "Error in _handle_submit() in New Endonuclease Controller", str(e))

    def _handle_reset_form(self):
        try:
            print("Button pressed")
            self.view.line_edit_organism.clear()
            self.view.line_edit_abbreviation.clear()
            self.view.line_edit_CRISPR_type.clear()
            self.view.line_edit_seed_length.clear()
            self.view.line_edit_PAM_sequence.clear()
        except Exception as e:
            show_error(self.settings, "Error in _handle_reset_form() in New Endonuclease Controller", str(e))

    def _handle_delete(self):
        try:
            selected = self.view.get_selected_endonuclease()
            if selected and selected != "Define New Endonuclease":
                self.model.delete_endonuclease(selected)
                show_message("Success", f"Endonuclease '{selected}' deleted successfully.")
                self._populate_endonuclease_combo_box()
                self.view.clear_form()
            else:
                show_message("Error", "Please select an existing endonuclease to delete.")
        except Exception as e:
            show_error(self.global_settings, "Error in _handle_delete() in New Endonuclease Controller", str(e))

    def validate_form_data(self, form_data, is_new_endonuclease):
        # Check for empty fields
        if any(value.strip() == "" for value in form_data.values()):
            show_message("Empty Field", "Please fill in all fields.")
            return False

        # Check for semicolons
        if any(';' in str(value) for value in form_data.values()):
            show_message("Invalid Character", "Invalid character used: ';'.")
            return False

        # Check PAM sequence
        valid_pam_chars = set('ACDEFGHIKLMNPQRSTVWY')
        if not set(form_data['endonuclease_pam_sequence']).issubset(valid_pam_chars):
            show_message("Invalid PAM", "Invalid characters in PAM Sequence.")
            return False

        # Check for duplicate endo abbreviations only for new endonucleases
        if is_new_endonuclease and self.model.is_duplicate_abbreviation(form_data['endonuclease_abbreviation']):
            show_message("Duplicate Abbreviation", "The given abbreviation already exists. Please choose a unique identifier.")
            return False

        return True

    def _on_endonuclease_updated(self):
        print("on_endonuclease_updated")
        self._populate_endonuclease_combo_box()
        self.settings.endonuclease_updated.emit()
