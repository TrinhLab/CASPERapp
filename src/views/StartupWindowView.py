from PyQt6 import QtWidgets, QtGui, QtCore, uic
from typing import Optional
from utils.ui import show_error
import os

class StartupWindowView(QtWidgets.QMainWindow):
    db_path_text_changed = QtCore.pyqtSignal(str)
    open_new_genome_requested = QtCore.pyqtSignal()

    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()

        self._init_ui()

    def _init_ui(self):
        try:
            uic.loadUi(os.path.join(self.settings.get_ui_dir_path(), 'startup_window.ui'), self)
            self._init_ui_components()
        except Exception as e:
            show_error(self.settings, "Error initializing StartupWindowView", str(e))
            raise

    def _init_ui_components(self):
        self._init_boxlayvTop()
        self._init_boxlayvMiddle()
        self._init_boxlayvBottom()

    def _init_boxlayvTop(self):
        self.label_logo = self._find_widget('lblLogo', QtWidgets.QLabel)
        self._set_startup_image()

    def _init_boxlayvMiddle(self):
        self.line_edit_database_directory = self._find_widget('ledDatabaseDirectory', QtWidgets.QLineEdit)
        self.push_button_change_directory = self._find_widget('pbtnChangeDirectory', QtWidgets.QPushButton)
        self.label_db_status = QtWidgets.QLabel(self)
        self.label_db_status.setWordWrap(True) 
        self.boxlayvMiddle.addWidget(self.label_db_status)

        self.line_edit_database_directory.textChanged.connect(self._on_db_path_text_changed)

    def _init_boxlayvBottom(self):
        self.push_button_go_to_home_or_new_genome = self._find_widget('pbtnGoToHomeOrNewGenome', QtWidgets.QPushButton)
        self.push_button_go_to_home_or_new_genome.clicked.connect(self._on_go_to_home_or_new_genome_clicked)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def _set_startup_image(self):
        image_path = os.path.join(self.settings.get_assets_dir_path(), "startup_image.jpg")
        if os.path.exists(image_path):
            pixmap = QtGui.QPixmap(image_path)
            self.label_logo.setPixmap(pixmap)
            self.label_logo.setScaledContents(True)
        else:
            self.logger.error(f"Startup image not found at {image_path}")
            self.label_logo.setText("Image not found")
    
    def set_db_path(self, directory_path):
        self.line_edit_database_directory.setText(directory_path)

    def get_db_path(self):
        return self.line_edit_database_directory.text()

    def _on_db_path_text_changed(self, new_path):
        self.db_path_text_changed.emit(new_path)

    def set_db_status(self, is_valid, message):
        if not hasattr(self, 'label_db_status'):
            return  # Exit if the label doesn't exist anymore
        
        if is_valid:
            self.label_db_status.hide()
            self.push_button_go_to_home_or_new_genome.setText("Go to Home")
        else:
            self.label_db_status.setText(message)
            self.label_db_status.show()
            self.label_db_status.setStyleSheet("color: red;")
            self.push_button_go_to_home_or_new_genome.setText("Analyze a New Genome")
        self.push_button_go_to_home_or_new_genome.setEnabled(True)

    def _on_go_to_home_or_new_genome_clicked(self):
        if self.push_button_go_to_home_or_new_genome.text() == "Analyze a New Genome":
            self.open_new_genome_requested.emit()
