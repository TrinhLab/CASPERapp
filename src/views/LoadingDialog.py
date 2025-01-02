from PyQt6.QtWidgets import QDialog, QProgressBar, QVBoxLayout, QLabel
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

class LoadingDialog(QDialog):
    def __init__(self, parent=None, message="Loading..."):
        super().__init__(parent)
        self.setWindowTitle("Please Wait")
        self.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.setFixedSize(300, 100)
        
        # Remove window decorations and set dialog flags
        self.setWindowFlags(Qt.WindowType.Dialog | Qt.WindowType.FramelessWindowHint | Qt.WindowType.WindowStaysOnTopHint)
        
        # Create layout
        layout = QVBoxLayout()
        
        # Add message label
        self.label = QLabel(message)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)  # Set range for percentage
        layout.addWidget(self.progress_bar)
        
        self.setLayout(layout)

    def showEvent(self, event):
        """Override show event to ensure dialog is centered when shown"""
        super().showEvent(event)
        self.center_on_parent()

    def center_on_parent(self):
        """Center the dialog on the parent window"""
        if self.parent():
            parent_geometry = self.parent().geometry()
            x = parent_geometry.x() + (parent_geometry.width() - self.width()) // 2
            y = parent_geometry.y() + (parent_geometry.height() - self.height()) // 2
            self.move(x, y)

    def set_message(self, message, progress=None):
        """Update the loading message and optionally the progress"""
        if progress is not None:
            self.progress_bar.setValue(progress)
        self.label.setText(message)
        self.center_on_parent()

    def set_progress(self, value):
        """Set progress value (0-100)"""
        self.progress_bar.setValue(value)
        self.label.setText("Loading...")

    def set_indeterminate(self):
        """Set indeterminate progress"""
        self.progress_bar.setRange(0, 0)
        self.label.setText("Loading...")
        