from PyQt6.QtWidgets import QMainWindow, QPushButton, QRadioButton, QComboBox, QPlainTextEdit, QProgressBar, QMenuBar, QMenu, QStackedWidget, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QFrame
from PyQt6.QtGui import QIcon, QAction, QFont
from PyQt6.QtCore import Qt, QPoint
from PyQt6 import uic, QtWidgets, QtCore
from utils.ui import scale_ui, show_error
import os
from typing import Optional

class MainWindowView(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._init_ui()

    def _init_ui(self) -> None:
        try:
            # self._load_ui_file()
            self._init_window_properties()
            self._init_custom_title_bar()
            self._init_ui_elements()
            self._scale_ui()
        except Exception as e:
            self._handle_init_error(e)

    def _load_ui_file(self) -> None:
        ui_file = os.path.join(self.global_settings.get_ui_dir(), "home_window.ui")
        uic.loadUi(ui_file, self)

    def _init_window_properties(self) -> None:
        self.setWindowFlags(Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)

    def _init_ui_elements(self) -> None:
        # Create a main widget to hold everything
        main_widget = QWidget()
        main_layout = QVBoxLayout(main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Add the custom title bar
        main_layout.addWidget(self.title_bar)

        # Create and add the divider
        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)
        divider.setStyleSheet("background-color: #c0c0c0;")  # Light gray color
        divider.setFixedHeight(1)  # 1 pixel height
        main_layout.addWidget(divider)

        # Create a widget to hold the original content
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        # Move the existing widgets to the content layout
        for child in self.children():
            if isinstance(child, (QMenuBar, QWidget)) and child != self.title_bar:
                content_layout.addWidget(child)

        # Create the stacked widget
        self.stacked_widget = QStackedWidget()

        # Add the content widget and stacked widget to the main layout
        main_layout.addWidget(content_widget)
        main_layout.addWidget(self.stacked_widget)

        # Set the main widget as the central widget
        self.setCentralWidget(main_widget)

        # Initialize other UI elements
        self._init_menuBar()
        self._init_grpNavigationMenu()
        self._init_grpStep1()
        self._init_grpStep2()
        self._init_grpStep3()

    def _init_custom_title_bar(self) -> None:
        self.title_bar = QWidget(self)
        self.title_bar.setObjectName("custom_title_bar")
        self.title_bar.setFixedHeight(32)  # Reduced height

        # Create the main horizontal layout for the title bar
        layout = QHBoxLayout(self.title_bar)
        layout.setContentsMargins(10, 0, 10, 0)  # Equal margins on left and right
        layout.setSpacing(5)  # Reduced spacing between items

        # ----- Window Control Buttons -----
        button_font = QFont("Arial", 8)  # Smaller font size for button text
        
        self.minimize_button = QPushButton("-", self.title_bar)
        self.minimize_button.setObjectName("minimize_button")
        self.minimize_button.setFixedSize(20, 20)  # Reduced size from 24x24 to 20x20
        self.minimize_button.setFont(button_font)

        self.maximize_button = QPushButton("⛶", self.title_bar)
        self.maximize_button.setObjectName("maximize_button")
        self.maximize_button.setFixedSize(20, 20)  # Reduced size from 24x24 to 20x20
        self.maximize_button.setFont(button_font)

        self.close_button = QPushButton("✕", self.title_bar)
        self.close_button.setObjectName("close_button")
        self.close_button.setFixedSize(20, 20)  # Reduced size from 24x24 to 20x20
        self.close_button.setFont(button_font)

        # Apply a style to center the text vertically and horizontally
        button_style = """
        QPushButton {
            padding: 0px;
            margin: 0px;
            line-height: 20px;
            text-align: center;
        }
        """
        self.minimize_button.setStyleSheet(button_style)
        self.maximize_button.setStyleSheet(button_style)
        self.close_button.setStyleSheet(button_style)

        # ----- Left Widget (Minimize, Maximize, Close) -----
        left_widget = QWidget()
        left_layout = QHBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(5)
        left_layout.addWidget(self.close_button)
        left_layout.addWidget(self.minimize_button)
        left_layout.addWidget(self.maximize_button)

        # ----- Theme Toggle Button -----
        self.theme_toggle_button = QPushButton(self.title_bar)
        self.theme_toggle_button.setObjectName("theme_toggle_button")
        self.theme_toggle_button.setFixedSize(20, 20)  # Reduced size from 24x24 to 20x20
        self.theme_toggle_button.setStyleSheet("border: none;")
        self.update_theme_icon()

        # ----- Right Widget (Theme Toggle + Stretch) -----
        right_widget = QWidget()
        right_layout = QHBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(5)

        # Add a stretch to push the toggle button to the far right within right_widget
        right_layout.addStretch()
        right_layout.addWidget(self.theme_toggle_button)

        # ----- Synchronize Widths of Left and Right Widgets -----
        # Adjust left_widget to calculate its required width
        left_widget.adjustSize()
        left_width = left_widget.sizeHint().width()

        # Set right_widget's fixed width to match left_widget's width
        right_widget.setFixedWidth(left_width)

        # ----- Title Label -----
        self.title_label = QLabel("CASPER", self.title_bar)
        self.title_label.setObjectName("title_label")
        self.title_label.setFont(QFont("Arial", 10, QFont.Weight.Bold))  # Reduced font size
        self.title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center the text in the label
        self.title_label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Preferred)

        # ----- Add Widgets to the Main Title Bar Layout -----
        layout.addWidget(left_widget)       # Left side buttons
        layout.addStretch(1)                 # Stretchable space
        layout.addWidget(self.title_label)  # Centered title
        layout.addStretch(1)                 # Stretchable space
        layout.addWidget(right_widget)      # Right side buttons (theme toggle + stretch)

        # Optional: Ensure the title_label is truly centered
        self.title_label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Preferred)
    def _init_menuBar(self) -> None:
        self.action_change_directory = self._find_widget("actChangeDirectory", QAction)
        self.action_open_genome_browser = self._find_widget("actOpenGenomeBrowser", QAction)
        self.action_open_repository = self._find_widget("actOpenRepository", QAction)
        self.action_open_NCBI_BLAST = self._find_widget("actOpenNCBIBLAST", QAction)
        self.action_open_NCBI = self._find_widget("actOpenNCBI", QAction)

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

    def _scale_ui(self) -> None:
        scale_ui(self, custom_scale_width=1000, custom_scale_height=350)

    def _handle_init_error(self, e: Exception) -> None:
        error_msg = f"Error initializing MainWindowView: {str(e)}"
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

    def update_theme_icon(self) -> None:
        icon_name = "dark_mode.png" if self.global_settings.is_dark_mode() else "light_mode.png"
        icon_path = os.path.join(self.global_settings.get_assets_dir(), icon_name)
        icon = QIcon(icon_path)
        self.theme_toggle_button.setIcon(icon)
        self.theme_toggle_button.setIconSize(QtCore.QSize(16, 16))  # Reduced icon size from 18x18 to 16x16

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton and self.title_bar.underMouse():
            self.drag_position = event.globalPosition().toPoint() - self.frameGeometry().topLeft()
            event.accept()

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.MouseButton.LeftButton and self.drag_position:
            self.move(event.globalPosition().toPoint() - self.drag_position)
            event.accept()

    def mouseReleaseEvent(self, event):
        self.drag_position = None