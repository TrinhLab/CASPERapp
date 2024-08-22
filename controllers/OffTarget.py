import os, platform
from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
from functools import partial
import models.GlobalSettings as GlobalSettings
from utils.ui import show_message, show_error, scale_ui, center_ui

logger = GlobalSettings.logger

class OffTarget(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(OffTarget, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'ui/off_target.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Off-Target Analysis")
            self.progressBar.setMinimum(0)
            self.progressBar.setMaximum(100)
            self.progressBar.setValue(0)
            self.run_clicked = False
            self.Run.clicked.connect(self.run_analysis)
            # self.tolerancehorizontalSlider.valueChanged.connect(self.tol_change)
            # self.tolerancehorizontalSlider.setMaximum(100)
            # self.tolerancehorizontalSlider.setMinimum(0)
            self.tolerance = 0.0

            self.cancelButton.clicked.connect(self.exit)
            self.fill_data_dropdown()
            self.perc = False
            self.bool_temp = False
            self.running = False
            self.process = QtCore.QProcess()
            self.output_path = ''

            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 5px;}
            QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            font: bold 14pt 'Arial';
                            margin-top: 10px;}"""

            self.Step1.setStyleSheet(groupbox_style)
            self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2"))
            self.Step3.setStyleSheet(groupbox_style.replace("Step1", "Step3"))

            scale_ui(self, custom_scale_width=400, custom_scale_height=450)

        except Exception as e:
            show_error("Error initializing OffTarget class.", e)

    #copied from MT to fill in the chromo and endo dropdowns based on CSPR files user provided at the startup
    def fill_data_dropdown(self):
        try:
            try:
                self.EndocomboBox.diconnect()
            except:
                pass
            try:
                self.OrgcomboBox.diconnect()
            except:
                pass

            self.OrgcomboBox.clear()
            self.EndocomboBox.clear()
            self.mismatchcomboBox.clear()

            self.organisms_to_files = {}
            self.organisms_to_endos = {}

            #fill in chromosome and endo dropdowns
            onlyfiles = [f for f in os.listdir(GlobalSettings.CSPR_DB) if os.path.isfile(os.path.join(GlobalSettings.CSPR_DB , f))]
            self.orgsandendos = {}
            self.shortName = {}
            for file in onlyfiles:
                if file.find('.cspr') != -1:
                    newname = file[0:-4]
                    endo = newname[newname.rfind("_") + 1:-1]
                    hold = open(file, 'r')
                    buf = (hold.readline())
                    hold.close()
                    buf = str(buf)
                    buf = buf.strip()
                    species = buf.replace("GENOME: ", "")

                    if species in self.organisms_to_files:
                        self.organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]
                    else:
                        self.organisms_to_files[species] = {}
                        self.organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]

                    if species in self.organisms_to_endos:
                        self.organisms_to_endos[species].append(endo)
                    else:
                        self.organisms_to_endos[species] = [endo]
                        if self.OrgcomboBox.findText(species) == -1:
                            self.OrgcomboBox.addItem(species)

            # fill in endos dropdown based on current organism
            endos = self.organisms_to_endos[str(self.OrgcomboBox.currentText())]
            self.EndocomboBox.addItems(endos)
            self.OrgcomboBox.currentIndexChanged.connect(self.update_endos)
            self.EndocomboBox.currentIndexChanged.connect(self.change_endos)

            # update file names for current org/endo combo
            self.cspr_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][endos[0]][0]
            self.db_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][endos[0]][1]

            #fill in Max Mismatch dropdown
            mismatch_list = ['1','2','3','4','5','6','7','8','9','10']
            self.mismatchcomboBox.addItems(mismatch_list)
            self.mismatchcomboBox.setCurrentIndex(3) ### Max number of mismatches is 4 by default
        except Exception as e:
            show_error("Error in fill_data_dropdown() in OffTarget.", e)

    def change_endos(self):
        try:
            #update file names based on current org/endo combo
            self.cspr_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][str(self.EndocomboBox.currentText())][0]
            self.db_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][str(self.EndocomboBox.currentText())][1]
        except Exception as e:
            show_error("Error in change_endos() in OffTarget.", e)
            
    def update_endos(self):
        try:
            #try to disconnect index changed signal on endo dropdown if there is one
            try:
                self.EndocomboBox.currentIndexChanged.disconnect()
            except:
                pass

            #clear endo dropdown and fill in with endos relative to the current organism
            self.EndocomboBox.clear()
            endos = self.organisms_to_endos[str(self.OrgcomboBox.currentText())]
            self.EndocomboBox.addItems(endos)
            self.cspr_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][endos[0]][0]
            self.db_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][endos[0]][1]

            #reconnect index changed signal on endo dropdown
            self.EndocomboBox.currentIndexChanged.connect(self.change_endos)
        except Exception as e:
            show_error("Error in update_endos() in OffTarget.", e)
            
    #tolerance slider / entry box. Allows for slider to update, or the user to input in text box
    # def tol_change(self):
    #     try:
    #         if(self.tolerance == float(self.tolerancelineEdit.text())):
    #             self.tolerance = self.tolerancehorizontalSlider.value() / 100 * 0.5
    #             self.tolerance = round(self.tolerance, 3)
    #             self.tolerancelineEdit.setText(str(self.tolerance))
    #         else:
    #             self.tolerance = float(self.tolerancelineEdit.text())
    #             self.tolerance = round(self.tolerance, 3)
    #             self.tolerancehorizontalSlider.setValue(round(self.tolerance/0.5 * 100))
    #     except Exception as e:
    #         logger.critical("Error in tol_change() in OffTarget.")
    #         logger.critical(e)
    #         logger.critical(traceback.format_exc())
    #         exit(-1)

    #run button linked to run_analysis, which is linked to the run button
    def run_command(self):
        try:
            output_file_name = self.outputFileName.text().strip()
            save_internally = self.outputCheckBox.isChecked()
            if not output_file_name and not save_internally:
                show_message(
                    fontSize=self.fontSize,
                    icon=QtWidgets.QMessageBox.Icon.Warning,
                    title="Input Error",
                    message="Please enter a valid output file name."
                )
                return

            if save_internally:
                full_output_path = os.path.join(GlobalSettings.appdir, 'local/local_output.txt')  
            else:
                full_output_path = os.path.join(GlobalSettings.CSPR_DB, output_file_name)
                if os.path.isfile(full_output_path):
                    show_message(
                        fontSize=self.fontSize,
                        icon=QtWidgets.QMessageBox.Icon.Critical,
                        title="Error",
                        message="Output file already exists. Please choose a new output file name."
                    )
                    return

            self.tolerance = self.toleranceSpinBox.value()
            self.perc = False
            self.bool_temp = False
            self.running = False

            app_path = GlobalSettings.appdir.replace('\\','/')
            exe_path = os.path.join(app_path, 'OffTargetFolder', 'OT_Win.exe' if platform.system() == 'Windows' else 'OT_Lin' if platform.system() == 'Linux' else 'OT_Mac')
            data_path = os.path.join(app_path, 'OffTargetFolder', 'temp.txt')
            endo = self.EndocomboBox.currentText()
            cspr_path = os.path.join(GlobalSettings.CSPR_DB, self.cspr_file)
            db_path = os.path.join(GlobalSettings.CSPR_DB, self.db_file)
            CASPER_info_path = os.path.join(app_path, 'CASPERinfo')
            num_of_mismatches = int(self.mismatchcomboBox.currentText())
            hsu = GlobalSettings.mainWindow.Results.endo_data[self.EndocomboBox.currentText()][2]

            self.output_path = ' "' + full_output_path + '"'
            
            cmd = f'"{exe_path}" "{data_path}" "{endo}" "{cspr_path}" "{db_path}" "{full_output_path}" "{CASPER_info_path}" {num_of_mismatches} {self.tolerance} {"TRUE" if self.AVG.isChecked() else "FALSE"} {"FALSE" if self.AVG.isChecked() else "TRUE"} "{hsu}"'   
            cmd = cmd.replace('/', '\\') if platform.system() == 'Windows' else cmd

            def finished():
                self.running = False
                self.run_clicked = True
                self.progressBar.setValue(100)

            #used to know when data is ready to read from stdout
            def dataReady():
                #filter the data from stdout, bools used to know when the .exe starts outputting the progress
                #percentages to be able to type cast them as floats and update the progress bar. Also, must
                #split the input read based on '\n\ characters since the stdout read can read multiple lines at
                #once and is all read in as raw bytes
                raw_data = self.process.readAllStandardOutput().data().decode()
                lines = raw_data.replace('\r', '').split('\n')
                for line in lines:
                    if "Parsing Input Arguments" in line:
                        self.progressBar.setValue(10)
                    elif "Loading data for algorithm" in line:
                        self.progressBar.setValue(25)
                    elif "Running OffTarget Analysis" in line:
                        self.progressBar.setValue(50)
                
            #connect QProcess to the dataReady func, and finished func, reset progressBar only if the outputfile name
            #given does not already exist
            self.process.readyReadStandardOutput.connect(partial(dataReady))
            self.process.readyReadStandardError.connect(partial(dataReady))
            self.progressBar.setValue(1)
            QtCore.QTimer.singleShot(100, partial(self.process.start, cmd))
            self.process.finished.connect(finished)
        except Exception as e:
            show_error("Error in run_command() in OffTarget.", e)

    def run_analysis(self):
        try:
            #make sure an analysis isn't already running before starting
            if(self.running == False):
                self.running = True
                self.run_command()
        except Exception as e:
            show_error("Error in run_analysis() in OffTarget.", e)
            
    #exit linked to user clicking cancel, resets bools, and kills process if one was running
    def exit(self):
        try:
            self.perc = False
            self.bool_temp = False
            self.running = False
            self.process.kill()
            self.hide()

            local_output_path = os.path.join(GlobalSettings.appdir, 'local/local_output.txt')
        
            if os.path.exists(local_output_path):
                os.remove(local_output_path)
                print("Local output file deleted successfully.")
            else:
                print("Local output file does not exist.")
        except Exception as e:
            show_error("Error in exit() in OffTarget.", e)
            
    #closeEvent linked to user pressing the x in the top right of windows, resets bools, and
    #kills process if there was one running
    def closeEvent(self, event):
        try:
            self.process.kill()
            self.perc = False
            self.bool_temp = False
            self.running = False
            event.accept()
        except Exception as e:
            show_error("Error in closeEvent() in OffTarget.", e)
