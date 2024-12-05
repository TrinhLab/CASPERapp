from typing import Optional
from PyQt6.QtWidgets import QMainWindow 
from PyQt6 import uic, QtWidgets

class ExportSelectedgRNAsView(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self._init_ui()

    def _init_ui(self) -> None:
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/export_selected_gRNAs.ui', self)
            self.setWindowTitle("Export Selected gRNAs")
            
            # Set fixed size for the window
            self.setFixedSize(500, 300)  # Width: 500px, Height: 300px
            
            self._init_ui_components()
        except Exception as e:
            self.logger.error(f"Error initializing ExportSelectedgRNAsView: {str(e)}", exc_info=True)
            raise

    def _init_ui_components(self) -> None:
        self._init_grpExportSettings()
        self._init_grpGuideOptions()

        self.push_button_cancel = self._find_widget('pbtnCancel', QtWidgets.QPushButton)
        self.push_button_export = self._find_widget('pbtnExport', QtWidgets.QPushButton)

    def _init_grpExportSettings(self) -> None:
        self.line_edit_file_path = self._find_widget('ledFilePath', QtWidgets.QLineEdit)
        self.push_button_browse = self._find_widget('pbtnBrowse', QtWidgets.QPushButton)
        self.line_edit_file_name = self._find_widget('ledFileName', QtWidgets.QLineEdit)
        self.combo_box_delimiter = self._find_widget('cmbDelimiter', QtWidgets.QComboBox)

    def _init_grpGuideOptions(self) -> None:
        self.line_edit_leading_sequence = self._find_widget('ledLeadingsequence', QtWidgets.QLineEdit)
        self.line_edit_trailing_sequence = self._find_widget('ledTrailingSequence', QtWidgets.QLineEdit)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget 

    def show_dialog(self) -> None:
        # Center the window on screen
        screen = self.screen().availableGeometry()
        self.move(
            screen.center().x() - self.width() // 2,
            screen.center().y() - self.height() // 2
        )
        
        self.show()
        self.activateWindow()

    def get_export_settings(self) -> dict:
        return {
            'file_path': self.line_edit_file_path.text(),
            'file_name': self.line_edit_file_name.text(),
            'delimiter': self.combo_box_delimiter.currentText(),
            'leading_sequence': self.line_edit_leading_sequence.text().strip(),
            'trailing_sequence': self.line_edit_trailing_sequence.text().strip()
        }

    def set_file_path(self, path: str) -> None:
        self.line_edit_file_path.setText(path)
