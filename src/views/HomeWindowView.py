from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QComboBox, QPlainTextEdit, QProgressBar, QRadioButton
from PyQt6 import uic, QtWidgets
from utils.ui import show_error
from typing import Optional

class HomeWindowView(QWidget):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._init_ui()

    def _init_ui(self) -> None:
        try:
            uic.loadUi(self.global_settings.get_ui_dir() + "/home_window_copy.ui", self)
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
        self.push_button_find_targets = self._find_widget("pbtnFindTargets", QPushButton)
        self.progress_bar_find_targets = self._find_widget("progBarFindTargets", QProgressBar)
        self.push_button_view_targets = self._find_widget("pbtnViewTargets", QPushButton)
        self.push_button_generate_library = self._find_widget("pbtnGenerateLibrary", QPushButton)

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

    def update_combo_box_annotation_files(self, annotation_files: list) -> None:
        self.combo_box_local_annotation_files.clear()
        self.combo_box_local_annotation_files.addItems(annotation_files)

    def set_progress_bar(self, value: int) -> None:
        self.progress_bar_find_targets.setValue(value)

    def reset_progress_bar(self) -> None:
        self.set_progress_bar(0)
