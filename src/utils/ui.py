from PyQt6 import QtWidgets 
import traceback
import sys
from PyQt6.QtWidgets import QMessageBox

def show_message(title, message, fontSize=12, icon=QtWidgets.QMessageBox.Icon.Information, button=QtWidgets.QMessageBox.StandardButton.Close):
    try:
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet(f"font: {fontSize}pt 'Arial'")
        msgBox.setIcon(icon)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.addButton(button)
        msgBox.exec()
    except Exception as e:
        print(f"Error showing message: {e}") 

def show_error(global_settings, message, exception=None):
    """Show error dialog and log the error"""
    logger = global_settings.get_logger() if global_settings else None
    
    if logger:
        logger.critical(message)
        if exception:
            if isinstance(exception, str):
                logger.critical(exception)
            else:
                logger.critical(str(exception))
                logger.critical(''.join(traceback.format_tb(exception.__traceback__)))

    error_box = QMessageBox()
    error_box.setIcon(QMessageBox.Icon.Critical)
    error_box.setText(message)
    if exception:
        if isinstance(exception, str):
            error_box.setDetailedText(exception)
        else:
            error_box.setDetailedText(f"{str(exception)}\n\n{''.join(traceback.format_tb(exception.__traceback__))}")
    error_box.exec()
    
    sys.exit(1) 
