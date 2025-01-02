from PyQt6 import QtWidgets, uic
from PyQt6.QtWidgets import QMainWindow, QFileDialog, QWidget
from PyQt6.QtCore import pyqtSignal, Qt
import os
import platform

class GenerateLibraryView(QMainWindow):
    # Define signals
    submit_clicked = pyqtSignal(dict)  # Signal to emit settings dict when submit is clicked
    progress_updated = pyqtSignal(int)
    
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.logger
        self._init_ui()
        self._connect_signals()
        
        # Set window properties
        self.setWindowModality(Qt.WindowModality.ApplicationModal)  # Make window modal
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)  # Clean up on close
        
    def _init_ui(self):
        try:
            # Load UI file
            ui_file = os.path.join(self.global_settings.get_ui_dir_path(), 'generate_library.ui')
            uic.loadUi(ui_file, self)
            
            # Set window properties
            self.setWindowTitle("Generate Library")
            self.setMinimumSize(800, 600)  # Set minimum size
            
            # Initialize comboboxes
            self._init_guides_per_gene_combo()
            self._init_min_score_combo()
            
            # Set default values
            self.ledTargetRangeStart.setText('0')
            self.ledTargetRangeEnd.setText('100')
            self.ledSpaceBetweenGuides.setText('15')
            
            # Set default file path
            default_path = self.global_settings.get_db_path()
            if platform.system() == "Windows":
                self.ledFilePath.setText(default_path + "\\")
            else:
                self.ledFilePath.setText(default_path + "/")
                
            # Apply styles
            self._set_styles()
                
            # Center the window
            self._center_window()
            
        except Exception as e:
            self.logger.error(f"Error initializing GenerateLibraryView UI: {str(e)}")
            raise
            
    def _center_window(self):
        """Center the window on the screen"""
        try:
            screen = QtWidgets.QApplication.primaryScreen().geometry()
            size = self.geometry()
            x = (screen.width() - size.width()) // 2
            y = (screen.height() - size.height()) // 2
            self.move(x, y)
        except Exception as e:
            self.logger.error(f"Error centering window: {str(e)}")
            
    def _init_guides_per_gene_combo(self):
        """Initialize guides per gene combobox"""
        try:
            self.cmbGuidesPerGene.clear()
            for i in range(1, 11):
                self.cmbGuidesPerGene.addItem(str(i))
        except Exception as e:
            self.logger.error(f"Error initializing guides per gene combo: {str(e)}")
            
    def _init_min_score_combo(self):
        """Initialize minimum on-target score combobox"""
        try:
            self.cmbMinimumOnTargetScore.clear()
            for i in range(20, 71):
                self.cmbMinimumOnTargetScore.addItem(str(i))
        except Exception as e:
            self.logger.error(f"Error initializing min score combo: {str(e)}")
            
    def _connect_signals(self):
        """Connect button signals"""
        try:
            self.pbtnBrowse.clicked.connect(self._browse_output_path)
            self.pbtnCancel.clicked.connect(self.close)
            self.pbtnSubmit.clicked.connect(self._on_submit)
        except Exception as e:
            self.logger.error(f"Error connecting signals: {str(e)}")
            
    def showEvent(self, event):
        """Override showEvent to ensure proper window display"""
        super().showEvent(event)
        self.raise_()  # Bring window to front
        self.activateWindow()  # Activate the window
        
    def closeEvent(self, event):
        """Override closeEvent to ensure proper cleanup"""
        try:
            self.logger.debug("Closing GenerateLibraryView")
            super().closeEvent(event)
        except Exception as e:
            self.logger.error(f"Error in closeEvent: {str(e)}")
        
    def _browse_output_path(self):
        """Handle browse button click"""
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            self.global_settings.get_db_path(),
            QFileDialog.Option.ShowDirsOnly
        )
        
        if folder:
            if platform.system() == "Windows":
                self.ledFilePath.setText(folder + "\\")
            else:
                self.ledFilePath.setText(folder + "/")
                
    def get_library_settings(self):
        """Get all settings from the UI"""
        try:
            settings = {
                'guides_per_gene': int(self.cmbGuidesPerGene.currentText()),
                'target_range_start': float(self.ledTargetRangeStart.text()),
                'target_range_end': float(self.ledTargetRangeEnd.text()),
                'space_between_guides': int(self.ledSpaceBetweenGuides.text() or '15'),
                'min_score': int(self.cmbMinimumOnTargetScore.currentText()),
                'find_off_targets': self.chkFindOffTargets.isChecked(),
                'modify_params': self.chkModifyParameters.isChecked(),
                'five_prime_seq': self.led5PrimeSpecificity.text(),
                'output_file': os.path.join(
                    self.ledFilePath.text(),
                    self.ledFileName.text()
                )
            }
            
            if settings['find_off_targets']:
                max_off_target_score = self.cmbMaximumOffTargetScore.text().strip()
                if not max_off_target_score:
                    raise ValueError("Please enter a maximum off-target score when Find Off Targets is enabled")
                try:
                    score = float(max_off_target_score)
                    if not 0 <= score <= 0.5:
                        raise ValueError("Maximum off-target score must be between 0 and 0.5 (inclusive)")
                    settings['max_off_target_score'] = score
                except ValueError as e:
                    if str(e).startswith("Maximum"):
                        raise
                    raise ValueError("Invalid maximum off-target score - please enter a valid number")
                    
            return settings
            
        except ValueError as e:
            raise ValueError(f"Invalid input: {str(e)}")
        except Exception as e:
            raise ValueError(f"Error getting settings: {str(e)}")
        
    def update_progress(self, value):
        """Update progress bar"""
        self.progBar.setValue(value)
        
    def _on_submit(self):
        """Validate and emit submit signal"""
        try:
            settings = self.get_library_settings()
            # Emit the settings through the signal
            self.submit_clicked.emit(settings)
        except ValueError as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Invalid Input",
                str(e)
            )
            
    def show_error(self, title, message):
        """Show error message"""
        QtWidgets.QMessageBox.critical(self, title, message)
        
    def show_success(self, message):
        """Show success message"""
        QtWidgets.QMessageBox.information(
            self,
            "Success",
            message
        )

    def _set_styles(self):
        """Apply the global groupbox style"""
        try:
            style = self.global_settings.get_groupbox_style()
            for groupbox in self.findChildren(QtWidgets.QGroupBox):
                groupbox.setStyleSheet(style)
        except Exception as e:
            self.logger.error(f"Error setting styles: {str(e)}")
