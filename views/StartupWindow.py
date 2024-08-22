import os
import platform
import traceback
import logging
from PyQt5 import QtGui, QtWidgets, QtCore, uic, Qt
import models.GlobalSettings as GlobalSettings
from utils.ui import show_message, show_error, scale_ui

logger = GlobalSettings.logger

class StartupWindow(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(StartupWindow, self).__init__()
            try:
                uic.loadUi(GlobalSettings.appdir + 'ui/startupCASPER.ui', self)
                self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
            except Exception as e:
                show_error("Unable to load UX files for Startup Window.", e)

            #set "Main" button to be the default highlighted button on startup
            self.goToMain.setDefault(True)

            #get current directory, and update based on current operating system
            self.currentDirectory = os.getcwd()
            self.databaseDirectory = self.loadDatabaseDirectory()
            GlobalSettings.CSPR_DB = self.databaseDirectory
            if platform.system() == "Windows":
                GlobalSettings.CSPR_DB = GlobalSettings.CSPR_DB.replace("/","\\")
            else:
                GlobalSettings.CSPR_DB = GlobalSettings.CSPR_DB.replace("\\","/")

            #setup event handlers for startup buttons
            self.currentDirText.setText(self.databaseDirectory)
            self.changeDir.clicked.connect(self.change_directory)
            self.goToMain.clicked.connect(self.launchMainWindow)
            self.goToNewGenome.clicked.connect(self.launchNewGenome)

            self.setWindowTitle("CASPER")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))

            scale_ui(self, custom_scale_width=1150, custom_scale_height=650)

        except Exception as e:
            show_error("Error initializing StartupWindow class.", e)

    #event handler for user clicking the "Change..." button - used for changing CASPER database directory
    def change_directory(self):
        try:
            # Launch OS file browser
            newDirectory = QtWidgets.QFileDialog.getExistingDirectory(
                self, "Open a folder...", self.databaseDirectory, QtWidgets.QFileDialog.ShowDirsOnly)

            # Check if selected path is a directory in the system
            if not os.path.isdir(newDirectory):
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Not a directory",
                    message="The directory you selected does not exist.",
                )
                return

            # Ensure directory contains correct filepath format based on OS
            newDirectory = newDirectory.replace("/", "\\") if platform.system() == "Windows" else newDirectory.replace("\\", "/")

            # Update text edit showing the current selected database directory
            self.currentDirText.setText(newDirectory)

            # Update CASPER database directories
            self.databaseDirectory = newDirectory
            GlobalSettings.CSPR_DB = newDirectory

        except Exception as e:
            show_error("change_directory() in startup window", e)

    #function for loading the default database directory specified in CASPERinfo
    #returns: default database parsed from CASPERinfo
    def loadDatabaseDirectory(self):
        casperInfoPath = os.path.join(GlobalSettings.appdir, "CASPERinfo")
        defaultDirectory = "Where would you like to store CASPER database files?"  # Default message if directory not found

        try:
            with open(casperInfoPath, 'r') as file:
                for line in file:
                    if 'DIRECTORY:' in line:
                        defaultDirectory = line.strip().replace("DIRECTORY:", "").strip()
                        break

            # Ensure the directory path is formatted correctly based on the operating system
            if platform.system() == "Windows":
                defaultDirectory = defaultDirectory.replace("/", "\\")
            else:
                defaultDirectory = defaultDirectory.replace("\\", "/")

            logger.debug("Successfully parsed CASPERinfo for default database directory.")
            return defaultDirectory

        except Exception as e:
            show_error(f"Error reading {casperInfoPath}: {e}", e)

        return defaultDirectory

    #function for saving the currently selected database directory to CASPERinfo to be the new default value on startup
    def saveDatabaseDirectory(self):
        try:
            #variable to hold the CASPERinfo data with new default directory change
            CASPERInfoNewData = ""

            #new default directory string for CASPERinfo
            newDefaultDirectory = "DIRECTORY:" + str(self.databaseDirectory)

            #open CASPERinfo file to read in the files data and add in new change
            try:
                CASPERInfo = open(GlobalSettings.appdir + "CASPERinfo", 'r+')
                CASPERinfoData = CASPERInfo.read()
                CASPERinfoData = CASPERinfoData.split('\n')
                for line in CASPERinfoData:
                    #if directory line found, use new default directory string instead
                    if 'DIRECTORY:' in line:
                        CASPERInfoNewData = CASPERInfoNewData + "\n" + newDefaultDirectory
                    else:
                        CASPERInfoNewData = CASPERInfoNewData + "\n" + line
                CASPERInfoNewData = CASPERInfoNewData[1:]

                #close CASPERinfo
                CASPERInfo.close()

                #re-open the file and re-write it with current changes
                CASPERInfo = open(GlobalSettings.appdir + "CASPERinfo", 'w+')
                CASPERInfo.write(CASPERInfoNewData)
                CASPERInfo.close()
                logger.debug("Successfully updated CASPERinfo with new default database directory.")
            except Exception as e:
                show_error("Unable to write to CASPERinfo file to update database directory.", e)
        except Exception as e:
            show_error("Error in saveDatabaseDirectory() in startup window.", e)

    # Event handler for user clicking the "New Genome" button - used for launching New Genome
    def launchNewGenome(self):
        try:
            # Make sure database directory variable is up-to-date based on what the user has in the text edit
            self.databaseDirectory = str(self.currentDirText.text())

            if not os.path.isdir(self.databaseDirectory):
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Not a directory",
                    message="The directory you selected does not exist.",
                )
                return

            # Change directories to the specified database directory provided
            os.chdir(self.databaseDirectory)

            # Write out the database directory to CASPERinfo to be the new default loaded value
            self.saveDatabaseDirectory()            

            # Update global database variable
            GlobalSettings.CSPR_DB = self.databaseDirectory

            # Create app directories
            initialize_app_directories()

            # Launch New Genome window
            self.launch_new_genome()

            self.close()
        except Exception as e:
            show_error("launchNewGenome() in startup window", e)

    def launch_new_genome(self):
        try:
            GlobalSettings.mainWindow.launch_newGenome()
            logger.debug("Successfully initialized New Genome in startup window.")
        except Exception as e:
            show_error("launch_new_genome() in startup window", e)

    # Event handler for user clicking "Main Program" button - used to launch Main Window
    def launchMainWindow(self):
        try:
            # Make sure database directory variable is up-to-date based on what the user has in the text edit
            self.databaseDirectory = str(self.currentDirText.text())

            # Make sure the path is a valid path before launching New Genome
            if not os.path.isdir(self.databaseDirectory):
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Not a directory",
                    message="The directory you selected does not exist.",
                )
                return

            # Check if database directory has CSPR files in it
            if not any(file.endswith(".cspr") for file in os.listdir(self.databaseDirectory)):
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="Directory is invalid!",
                    message="You must select a directory with CSPR Files!",
                )
                return

            # Change directory to database directory
            os.chdir(self.databaseDirectory)

            # Update database directory global variable
            GlobalSettings.CSPR_DB = self.databaseDirectory

            # Save database directory to CASPERinfo
            self.saveDatabaseDirectory()

            initialize_app_directories()

            # Fill in organism/endo/annotation dropdown information for main, mulit-targeting, and populatin analysis
            self.load_dropdown_data()

            # Show main window
            if GlobalSettings.mainWindow.first_show:
                GlobalSettings.mainWindow.first_show = False
            GlobalSettings.mainWindow.show()
            self.close()
            
        except Exception as e:
            show_error("launchMainWindow() in startup window", e)

    def load_dropdown_data(self):
        try:
            GlobalSettings.mainWindow.getData()
            GlobalSettings.mainWindow.fill_annotation_dropdown()
            logger.debug("Successfully loaded organism/endo/annotation drop down information in Main.")
        except Exception as e:
            show_error("load_dropdown_data() in Main", e)

        try:
            GlobalSettings.MTWin.launch()
            logger.debug("Successfully loaded organism/endo drop down information in Multi-targeting.")
        except Exception as e:
            show_error("load_dropdown_data() in Multi-targeting", e)

        try:
            GlobalSettings.pop_Analysis.launch()
            logger.debug("Successfully loaded organism/endo drop down information in Population Analysis.")
        except Exception as e:
            show_error("load_dropdown_data() in Population Analysis", e)

def initialize_app_directories():
    required_dirs = ["FNA", "GBFF"]
    for directory in required_dirs:
        path = os.path.join(GlobalSettings.CSPR_DB, directory)
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logging.info(f"Directory created: {path}")