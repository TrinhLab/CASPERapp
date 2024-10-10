from PyQt6 import QtWidgets, QtGui, QtCore, uic
from typing import Optional
import os

class StartupWindowView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()

        self._init_ui()

    def _init_ui(self):
        uic.loadUi(os.path.join(self.global_settings.get_ui_dir(), 'startup_window.ui'), self)

        self._init_ui_components()
        self.setMinimumSize(self.sizeHint())

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

    def _init_boxlayvBottom(self):
        self.push_button_go_to_home_or_new_genome = self._find_widget('pbtnGoToHomeOrNewGenome', QtWidgets.QPushButton)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def _set_startup_image(self):
        image_path = os.path.join(self.global_settings.get_assets_dir(), "startup_image.jpg")
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

    def set_valid_db_state(self):
        self.logger.debug("Setting valid DB state")
        self.label_db_status.setText("")
        self.label_db_status.hide()
        self.push_button_go_to_home_or_new_genome.setText("Go to Home")
        self.push_button_go_to_home_or_new_genome.setEnabled(True)
        self.adjustSize()

    def set_invalid_db_state(self):
        self.logger.debug("Setting invalid DB state")
        self.label_db_status.setText("You must select a directory with CSPR files. Define a new genome or change directory.")
        self.label_db_status.show()
        self.push_button_go_to_home_or_new_genome.setText("Define a New Genome")
        self.push_button_go_to_home_or_new_genome.setEnabled(True)
        self.adjustSize()

    def adjustSize(self):
        super().adjustSize()
        self.setFixedSize(self.sizeHint())
