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
        
        # Center on main window
        self.center_on_parent()

    def center_on_parent(self):
        """Center the dialog on the main window or parent"""
        parent = self.parent()
        if parent:
            # Get the main window from global settings if available
            main_window = None
            if hasattr(parent, 'global_settings'):
                main_window = parent.global_settings.main_window
            elif hasattr(parent, 'settings'):
                main_window = parent.settings.main_window

            # Get geometry of the window to center on
            if main_window and main_window.view:
                geometry = main_window.view.geometry()
            else:
                geometry = parent.geometry()

            # Calculate center position
            x = geometry.x() + (geometry.width() - self.width()) // 2
            y = geometry.y() + (geometry.height() - self.height()) // 2
            
            # Ensure dialog stays within screen bounds
            screen = QApplication.primaryScreen().geometry()
            x = max(screen.left(), min(x, screen.right() - self.width()))
            y = max(screen.top(), min(y, screen.bottom() - self.height()))
            
            self.move(x, y)

    def set_message(self, message, progress=None):
        """Update the loading message and optionally the progress"""
        if progress is not None:
            self.progress_bar.setValue(progress)
        self.label.setText(message)
        
        # Recenter after updating message
        self.center_on_parent()

    def set_progress(self, value):
        """Set progress value (0-100)"""
        self.progress_bar.setValue(value)
        self.label.setText("Loading...")

    def set_indeterminate(self):
        """Set indeterminate progress"""
        self.progress_bar.setRange(0, 0)
        self.label.setText("Loading...")

    def showEvent(self, event):
        """Override show event to ensure dialog is centered when shown"""
        super().showEvent(event)
        self.center_on_parent()
        