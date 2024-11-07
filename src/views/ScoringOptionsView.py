from PyQt6 import QtWidgets, uic
from PyQt6.QtCore import pyqtSignal
import traceback

class ScoringOptionsView(QtWidgets.QMainWindow):
    # Define signals
    fasta_selected = pyqtSignal(str)  # Signal when FASTA file is selected
    submit_clicked = pyqtSignal()  # Signal when submit button is clicked
    
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self._init_ui()

    def _init_ui(self):
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/scoring_options.ui', self)
            
            # Get UI elements
            self.push_button_browse = self.findChild(QtWidgets.QPushButton, 'pbtnBrowse')
            self.push_button_submit = self.findChild(QtWidgets.QPushButton, 'pbtnSubmit')
            self.line_edit_fasta = self.findChild(QtWidgets.QLineEdit, 'ledInputFASTA')
            self.radio_button_azimuth = self.findChild(QtWidgets.QRadioButton, 'rbtnAzimuth')
            
            # Connect signals
            self.push_button_browse.clicked.connect(self._browse_fasta)
            self.push_button_submit.clicked.connect(self.submit_clicked.emit)
            
            # Set window title
            self.setWindowTitle("Select Scoring Algorithm")
            
            # Apply theme
            self.apply_theme()
            
        except Exception as e:
            self.logger.error(f"Error initializing ScoringOptionsView: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            raise

    def apply_theme(self):
        """Apply the current theme to the window"""
        try:
            current_theme = self.settings.get_theme()
            themes = {
                "dark": {
                    "bg_color": "#2b2b2b",
                    "fg_color": "#ffffff",
                    "button_bg_color": "#3a3a3a",
                    "button_border_color": "#5a5a5a",
                    "button_hover_bg_color": "#4a4a4a",
                    "input_bg_color": "#3a3a3a",
                    "input_border_color": "#5a5a5a",
                    "progress_bar_bg": "#3a3a3a",
                    "progress_bar_chunk": "#51b85e"
                },
                "light": {
                    "bg_color": "#f0f0f0",
                    "fg_color": "#000000",
                    "button_bg_color": "#e0e0e0",
                    "button_border_color": "#c0c0c0",
                    "button_hover_bg_color": "#d0d0d0",
                    "input_bg_color": "#ffffff",
                    "input_border_color": "#c0c0c0",
                    "progress_bar_bg": "#e0e0e0",
                    "progress_bar_chunk": "#51b85e"
                }
            }
            
            theme = themes["dark"] if current_theme == "dark" else themes["light"]
            
            # Set the stylesheet
            self.setStyleSheet(f"""
                QMainWindow, QWidget {{ 
                    background-color: {theme['bg_color']}; 
                    color: {theme['fg_color']}; 
                }}
                QPushButton {{ 
                    background-color: {theme['button_bg_color']}; 
                    border: 1px solid {theme['button_border_color']}; 
                    padding: 5px;
                    min-width: 80px;
                }}
                QPushButton:hover {{ 
                    background-color: {theme['button_hover_bg_color']}; 
                }}
                QLineEdit {{ 
                    background-color: {theme['input_bg_color']}; 
                    border: 1px solid {theme['input_border_color']}; 
                    padding: 5px;
                }}
                QRadioButton {{ 
                    color: {theme['fg_color']}; 
                }}
                QProgressBar {{
                    border: 1px solid {theme['button_border_color']};
                    background-color: {theme['progress_bar_bg']};
                    text-align: center;
                }}
                QProgressBar::chunk {{
                    background-color: {theme['progress_bar_chunk']};
                }}
                QGroupBox {{ 
                    border: 1px solid {theme['button_border_color']};
                    margin-top: 0.5em;
                    padding-top: 0.5em;
                }}
                QGroupBox::title {{
                    color: {theme['fg_color']};
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 3px 0 3px;
                }}
            """)
            
        except Exception as e:
            self.logger.error(f"Error applying theme: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def _browse_fasta(self):
        try:
            # Get database directory path
            db_path = self.settings.get_db_path()
            
            file_dialog = QtWidgets.QFileDialog()
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(
                file_dialog,
                "Choose FASTA File",
                db_path,  # Set initial directory to database path
                "FASTA Files (*.fa *.fasta *.fna)"
            )
            
            if file_path:
                self.line_edit_fasta.setText(file_path)
                self.fasta_selected.emit(file_path)
                self.logger.debug(f"Selected FASTA file: {file_path}")
                
        except Exception as e:
            self.logger.error(f"Error browsing FASTA file: {str(e)}")
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                f"Error selecting FASTA file: {str(e)}"
            )

    def get_selected_algorithm(self):
        """Get the currently selected scoring algorithm"""
        if self.radio_button_azimuth.isChecked():
            return "Azimuth 2.0"
        return None

    def get_fasta_path(self):
        """Get the selected FASTA file path"""
        return self.line_edit_fasta.text()

    def show_error(self, title, message):
        """Show error message box"""
        QtWidgets.QMessageBox.critical(self, title, message)

    def show_info(self, title, message):
        """Show info message box"""
        QtWidgets.QMessageBox.information(self, title, message)
