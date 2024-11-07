from PyQt6 import QtWidgets, QtGui, QtCore
import traceback
import models.GlobalSettings as GlobalSettings


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

def show_error(settings, message, e):
    try:
        logger = settings.get_logger()
        logger.critical(message)
        logger.critical(e)
        logger.critical(traceback.format_exc())

        show_message(
            fontSize=12,
            icon=QtWidgets.QMessageBox.Icon.Critical,
            title="Fatal Error",
            message=f"Fatal Error:\n{str(e)}\n\nFor more information on this error, look at CASPER.log in the application folder."
        )
    except Exception as e:
        print(f"Error showing error message: {e}") 

    exit(-1)
