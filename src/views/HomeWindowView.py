from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QComboBox, QPlainTextEdit, QProgressBar, QRadioButton
from PyQt6 import uic, QtWidgets
from utils.ui import show_error
from typing import Optional
import os
class HomeWindowView(QWidget):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = self.global_settings.logger
        self._init_ui()

    def _init_ui(self) -> None:
        try:
            uic.loadUi(os.path.join(self.global_settings.get_ui_dir_path(), "home_window.ui"), self)
            self._init_ui_elements()
        except Exception as e:
            self._handle_init_error(e)

    def _init_ui_elements(self) -> None:
        # Create a main layout to hold everything
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Create a widget to hold the original content
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        for child in self.children():
            if isinstance(child, QWidget):
                content_layout.addWidget(child)

        main_layout.addWidget(content_widget)

        self._init_grpNavigationMenu()
        self._init_grpStep1()
        self._init_grpStep2()
        self._init_grpStep3()

        # Connect to database manager signals
        self.global_settings.db_manager.db_files_changed.connect(self._handle_db_files_changed)
        self.global_settings.db_manager.db_state_changed.connect(self._handle_db_state_changed)

    def _init_grpNavigationMenu(self) -> None:
        self.push_button_new_genome = self._find_widget("pbtnNewGenome", QPushButton)
        self.push_button_new_endonuclease = self._find_widget("pbtnNewEndonuclease", QPushButton)
        self.push_button_multitargeting_analysis = self._find_widget("pbtnMultitargetingAnalysis", QPushButton)
        self.push_button_population_analysis = self._find_widget("pbtnPopulationAnalysis", QPushButton)
        self.push_button_combine_files = self._find_widget("pbtnCombineFiles", QPushButton)

    def _init_grpStep1(self) -> None:
        self.combo_box_organism = self._find_widget("cmbOrganism", QComboBox)
        self.combo_box_endonuclease = self._find_widget("cmbEndonuclease", QComboBox)

    def _init_grpStep2(self) -> None:
        self.push_button_ncbi_file_search = self._find_widget("pbtnNCBIFileSearch", QPushButton)
        self.combo_box_local_annotation_files = self._find_widget("cmbLocalAnnotationFiles", QComboBox)

    def _init_grpStep3(self) -> None:
        self.radio_button_feature = self._find_widget("rbtnFeature", QRadioButton)
        self.radio_button_position = self._find_widget("rbtnPosition", QRadioButton)
        self.radio_button_sequence = self._find_widget("rbtnSequence", QRadioButton)
        self.text_edit_gene_entry = self._find_widget("txtedGeneEntry", QPlainTextEdit)
        self.push_button_find_view_targets = self._find_widget("pbtnFindViewTargets", QPushButton)

        placeholder_text = ("Example Inputs: \n\n"
                            "Option 1: Feature (ID, Locus Tag, or Name)\n"
                            "Example: 854068/YOL086C/ADH1 for S. cerevisiae alcohol dehydrogenase 1\n\n"
                            "Option 2: Position (chromosome,start,stop)\n"
                            "Example: 1,1,1000 for targeting chromosome 1, base pairs 1 to 1000\n\n"
                            "Option 3: Sequence (must be within the selected organism)\n"
                            "Example: Any nucleotide sequence between 100 and 10,000 base pairs.\n\n"
                            "*Note: to multiplex, separate multiple queries by new lines*\n"
                            "Example:\n"
                            "1,1,1000\n"
                            "5,1,500\n"
                            "etc.")
        self.text_edit_gene_entry.setPlaceholderText(placeholder_text)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.global_settings.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def _handle_init_error(self, e: Exception) -> None:
        error_msg = f"Error initializing HomeWindowView: {str(e)}"
        self.global_settings.logger.error(error_msg, exc_info=True)
        show_error(self.global_settings, "Initialization Error", error_msg)
        raise

    def update_combo_box_endonuclease(self, endonuclease: list) -> None:
        self.combo_box_endonuclease.clear()
        self.combo_box_endonuclease.addItems(endonuclease)

    def update_combo_box_organism(self, organisms: list) -> None:
        self.combo_box_organism.clear()
        self.combo_box_organism.addItems(organisms)

    # def update_combo_box_annotation_files(self, annotation_files: list) -> None:
    #     self.combo_box_local_annotation_files.clear()
    #     self.combo_box_local_annotation_files.addItems(annotation_files)

    def get_find_targets_input(self) -> dict:
        return {
            "organism": self.combo_box_organism.currentText(),
            "endonuclease": self.combo_box_endonuclease.currentText(),
            "annotation_file": self.combo_box_local_annotation_files.currentText(),
            "search_type": self.get_search_type(),
            "search_query": self.text_edit_gene_entry.toPlainText()
        }

    def get_search_type(self) -> str:
        if self.radio_button_feature.isChecked():
            return "feature"
        elif self.radio_button_position.isChecked():
            return "position"
        elif self.radio_button_sequence.isChecked():
            return "sequence"
        else:
            return "feature"  # Default to feature if somehow none are selected
        
    def get_annotation_file(self) -> str:
        return self.combo_box_local_annotation_files.currentText()

    def update_combo_box_annotation_files(self, files):
        """Update local annotation files combo box, excluding .index files"""
        try:
            # Clear existing items
            self.combo_box_local_annotation_files.clear()
            
            # Filter out .index files
            filtered_files = [f for f in files if not f.endswith('.index')]
            
            # Add filtered files to combo box
            if filtered_files:
                self.combo_box_local_annotation_files.addItems(filtered_files)
                self.combo_box_local_annotation_files.setCurrentIndex(0)
                self.logger.debug(f"Added {len(filtered_files)} local annotation files to combo box")
            else:
                self.logger.debug("No local annotation files found")
                
        except Exception as e:
            self.logger.error(f"Error updating local annotation files: {str(e)}")

    def show_warning(self, title: str, message: str) -> None:
        """Show a warning message dialog"""
        QtWidgets.QMessageBox.warning(self, title, message)

    def _update_cspr_related_ui(self) -> None:
        """Update UI elements that depend on CSPR files"""
        try:
            # Store current selections
            current_organism = self.combo_box_organism.currentText()
            current_endo = self.combo_box_endonuclease.currentText()
            
            # Get fresh data from controller
            controller = self.global_settings.main_window.controller
            organism_to_endonuclease = controller.get_organism_to_endonuclease()
            
            # Update organism combo box
            self.combo_box_organism.clear()
            self.combo_box_organism.addItems(sorted(organism_to_endonuclease.keys()))
            
            # Restore organism selection if still valid
            if current_organism in organism_to_endonuclease:
                self.combo_box_organism.setCurrentText(current_organism)
                # Restore endonuclease selection if still valid for this organism
                if current_endo in organism_to_endonuclease[current_organism]:
                    self.combo_box_endonuclease.setCurrentText(current_endo)
                
        except Exception as e:
            self.logger.error(f"Error updating CSPR-related UI: {str(e)}")

    def _update_gbff_related_ui(self) -> None:
        """Update UI elements that depend on GBFF files"""
        try:
            # Store current selection
            current_file = self.combo_box_local_annotation_files.currentText()
            
            # Update annotation files
            annotation_files = self.global_settings.main_window.controller.get_annotation_files()
            self.update_combo_box_annotation_files(annotation_files)
            
            # Restore selection if still valid
            if current_file in annotation_files:
                self.combo_box_local_annotation_files.setCurrentText(current_file)
                
        except Exception as e:
            self.logger.error(f"Error updating GBFF-related UI: {str(e)}")

    def _handle_db_files_changed(self, changes):
        """Handle database file changes"""
        try:
            if (FileChangeType.CSPR_ADDED in changes or 
                FileChangeType.CSPR_REMOVED in changes):
                self._update_cspr_related_ui()
                
            if (FileChangeType.GBFF_ADDED in changes or 
                FileChangeType.GBFF_REMOVED in changes):
                self._update_gbff_related_ui()
                
        except Exception as e:
            self.logger.error(f"Error handling database file changes: {str(e)}")

    def _handle_db_state_changed(self, is_valid, message, changes):
        """Handle database state changes"""
        try:
            if not is_valid:
                self.show_warning("Database Warning", message)
                return
            
            if changes:  # If there are any changes
                if any(change in changes for change in 
                      [FileChangeType.CSPR_ADDED, FileChangeType.CSPR_REMOVED]):
                    self._update_cspr_related_ui()
                    
                if any(change in changes for change in 
                      [FileChangeType.GBFF_ADDED, FileChangeType.GBFF_REMOVED]):
                    self._update_gbff_related_ui()
                    
        except Exception as e:
            self.logger.error(f"Error handling database state change: {str(e)}")