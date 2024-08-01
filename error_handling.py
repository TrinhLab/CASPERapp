from PyQt5 import QtWidgets
import traceback
import GlobalSettings
from common_utils import show_message

#global logger
logger = GlobalSettings.logger

def show_error(message, e):
    logger.critical(message)
    logger.critical(e)
    logger.critical(traceback.format_exc())

    show_message(
        fontSize=12,
        icon=QtWidgets.QMessageBox.Icon.Critical,
        title="Fatal Error",
        message=f"Fatal Error:\n{str(e)}\n\nFor more information on this error, look at CASPER.log in the application folder."
    )

    exit(-1)