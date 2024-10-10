from PyQt6 import QtWidgets, QtGui, QtCore, uic
import os

class MainWindowUI(QtWidgets.QMainWindow):
    def __init__(self, settings):
        super(MainWindowUI, self).__init__()
        self.settings = settings
        self.setup_ui()

    def setup_ui(self):
        # Load the UI file
        uic.loadUi(os.path.join(self.settings.get_ui_dir(), 'main_window.ui'), self)

        # Set window properties
        self.setWindowTitle("CASPER")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.settings.get_assets_dir(), "cas9image.ico")))

        # Initialize UI components
        self.init_ui_components()

        # Set up styles
        self.set_styles()

        # Initialize progress bar to 0
        self.progressBar.setValue(0)

        # Add the theme toggle button
        self.theme_toggle_button = QtWidgets.QPushButton(self)
        self.theme_toggle_button.setFixedSize(32, 32)
        self.theme_toggle_button.setStyleSheet("border: none;")
        self.update_theme_icon()

        # Position the button in the top right corner
        self.theme_toggle_button.setGeometry(self.width() - 40, 10, 32, 32)

        # Connect the button to a slot (to be implemented in the controller)
        self.theme_toggle_button.clicked.connect(self.on_theme_toggle)

    def init_ui_components(self):
        # Initialize and find all the UI components
        self.org_choice = self.findChild(QtWidgets.QComboBox, 'orgChoice')
        self.endo_choice = self.findChild(QtWidgets.QComboBox, 'endoChoice')
        self.annotation_files = self.findChild(QtWidgets.QComboBox, 'annotation_files')
        self.gene_entry_field = self.findChild(QtWidgets.QPlainTextEdit, 'gene_entry_field')

        self.push_button_find_targets = self.findChild(QtWidgets.QPushButton, 'pushButton_FindTargets')
        self.push_button_view_targets = self.findChild(QtWidgets.QPushButton, 'pushButton_ViewTargets')
        self.generate_library = self.findChild(QtWidgets.QPushButton, 'GenerateLibrary')
        
        self.radio_button_gene = self.findChild(QtWidgets.QRadioButton, 'radioButton_Gene')
        self.radio_button_position = self.findChild(QtWidgets.QRadioButton, 'radioButton_Position')
        self.radio_button_sequence = self.findChild(QtWidgets.QRadioButton, 'radioButton_Sequence')
        
        self.new_genome_button = self.findChild(QtWidgets.QPushButton, 'newGenome_button')
        self.new_endo_button = self.findChild(QtWidgets.QPushButton, 'newEndo_button')
        self.multitargeting_button = self.findChild(QtWidgets.QPushButton, 'multitargeting_button')
        self.population_analysis_button = self.findChild(QtWidgets.QPushButton, 'populationAnalysis_button')
        self.combine_files_button = self.findChild(QtWidgets.QPushButton, 'combineFiles_button')
        
        self.progress_bar = self.findChild(QtWidgets.QProgressBar, 'progressBar')
        
        self.step1 = self.findChild(QtWidgets.QGroupBox, 'Step1')
        self.step2 = self.findChild(QtWidgets.QGroupBox, 'Step2')
        self.step3 = self.findChild(QtWidgets.QGroupBox, 'Step3')
        self.casper_navigation = self.findChild(QtWidgets.QGroupBox, 'CASPER_Navigation')
        
        self.ncbi_button = self.findChild(QtWidgets.QPushButton, 'ncbi_button')

        # Connect the actionChange_Directory to a slot
        self.actionChange_Directory.triggered.connect(self.on_change_directory)

    def set_styles(self):
        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        margin-top: 10px;
                        font: bold 14pt 'Arial';}
        """
        self.step1.setStyleSheet(groupbox_style)
        self.step2.setStyleSheet(groupbox_style.replace("Step1", "Step2"))
        self.step3.setStyleSheet(groupbox_style.replace("Step1", "Step3"))
        self.casper_navigation.setStyleSheet(groupbox_style.replace("Step1", "CASPER_Navigation")
                                             .replace("solid","dashed")
                                             .replace("rgb(111,181,110)","rgb(88,89,91)"))

    def set_gene_entry_placeholder(self):
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
        self.gene_entry_field.setPlaceholderText(placeholder_text)

    def enable_view_targets(self, enable):
        self.push_button_view_targets.setEnabled(enable)

    def enable_generate_library(self, enable):
        self.generate_library.setEnabled(enable)

    def set_progress(self, value):
        if self.progressBar:
            self.progressBar.setValue(value)

    def reset_progress(self):
        if self.progressBar:
            self.progressBar.setValue(0)

    def toggle_annotation(self, gene_checked):
        self.step2.setEnabled(True)

    def update_endo_choice(self, endos):
        self.endo_choice.clear()
        self.endo_choice.addItems(endos)

    def bring_to_front(self):
        self.show()
        self.setWindowState(self.windowState() & ~QtCore.Qt.WindowState.WindowMinimized | QtCore.Qt.WindowState.WindowActive)
        self.raise_()
        self.activateWindow()
        QtWidgets.QApplication.setActiveWindow(self)

    def on_change_directory(self):
        # This method will be connected to the controller
        pass

    def update_theme_icon(self):
        # Swap the icons: use dark_mode.png for light mode, and light_mode.png for dark mode
        icon = QtGui.QIcon(os.path.join(self.settings.get_assets_dir(), "dark_mode.png" if self.settings.is_dark_mode() else "light_mode.png"))
        self.theme_toggle_button.setIcon(icon)
        self.theme_toggle_button.setIconSize(QtCore.QSize(24, 24))
        
        # Update the entire application's theme
        self.settings.apply_theme()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Reposition the theme toggle button when the window is resized
        self.theme_toggle_button.setGeometry(self.width() - 40, 10, 32, 32)

    def on_theme_toggle(self):
        # This method will be connected to the controller
        pass
