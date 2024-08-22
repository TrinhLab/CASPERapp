import sys
import os
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import models.GlobalSettings as GlobalSettings
import controllers.multitargeting as multitargeting
import controllers.populationAnalysis as populationAnalysis
import platform
import logging
from utils.ui import show_error, center_ui
from views.annotation_functions import *
from views.StartupWindow import StartupWindow
from views.CMainWindow import CMainWindow

logger = GlobalSettings.logger

fontSize = 12

def setup_logger():
    logger.info(f"System OS: {platform.system()}")

    if hasattr(sys, 'frozen'):
        logger.info("Running a packaged version of CASPER.")
        GlobalSettings.appdir = sys.executable
        if platform.system() == 'Windows':
            GlobalSettings.appdir = sys._MEIPASS + "\\"
        else:
            GlobalSettings.appdir = GlobalSettings.appdir[:GlobalSettings.appdir.rfind("Contents/") + 9] + "Resources/"

    else:
        logger.info("Running a non-packaged version of CASPER.")
        GlobalSettings.appdir = os.path.dirname(os.path.abspath(__file__)) + ('\\' if platform.system() == 'Windows' else '/')

    fh = logging.FileHandler(GlobalSettings.appdir + 'logs/app.log', mode='w')
    fh_formatter = logging.Formatter('%(asctime)s %(levelname)s %(lineno)d:%(filename)s(%(process)d) - %(message)s')
    fh.setFormatter(fh_formatter)
    fh.setLevel(logging.DEBUG)
    GlobalSettings.logger.addHandler(fh)

# def update_directory(current_directory, fontSize, update_ui_callback):
#     try:
#         filed = QtWidgets.QFileDialog()
#         new_directory = QtWidgets.QFileDialog.getExistingDirectory(
#             filed, "Open a folder...", current_directory, QtWidgets.QFileDialog.ShowDirsOnly)

#         if not os.path.isdir(new_directory):
#             show_message("Not a directory", "The directory you selected does not exist.", fontSize)
#             return

#         if not any(file.endswith(".cspr") for file in os.listdir(new_directory)):
#             show_message("Directory is invalid!", "You must select a directory with CSPR Files!", fontSize)
#             return

#         if platform.system() == "Windows":
#             new_directory = new_directory.replace("/", "\\")

#         os.chdir(new_directory)
#         GlobalSettings.CSPR_DB = new_directory
#         update_ui_callback(new_directory)

#     except Exception as e:
#         show_critical_error("Error in updating directory", e, fontSize)

def main():
    setup_logger()
    app = QtWidgets.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")

    #load startup window
    try:
        # Initialize windows
        startup = StartupWindow()
        logger.debug("Successfully initialized Startup Window.")

        GlobalSettings.mainWindow = CMainWindow(os.getcwd())
        logger.debug("Successfully initialized Main Window.")

        GlobalSettings.MTWin = multitargeting.Multitargeting()
        logger.debug("Successfully initialized Multi-targeting Window.")

        GlobalSettings.pop_Analysis = populationAnalysis.Pop_Analysis()
        logger.debug("Successfully initialized Population Analysis Window.")

        center_ui(startup)
        startup.show()
        sys.exit(app.exec_())
    except Exception as e:
        show_error("An error occurred during application initialization", e)

if __name__ == '__main__':
    main()
