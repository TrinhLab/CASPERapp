import os, platform
from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
from functools import partial
import GlobalSettings
import gzip
import traceback
import math

#global logger
logger = GlobalSettings.logger

class OffTarget(QtWidgets.QMainWindow):

    def __init__(self):
        try:
            super(OffTarget, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'off_target.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Off-Target Analysis")
            self.progressBar.setMinimum(0)
            self.progressBar.setMaximum(100)
            self.progressBar.setValue(0)
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

            # make sure to intialize the class variable in init. That way elsewhere and other classes can access it
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

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing OffTarget class.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #scale UI based on current screen
    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.fontSize = fontSize
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            # CASPER header scaling
            fontSize = 20
            self.title.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # window scaling
            # 1920x1080 => 850x750
            scaledWidth = int((width * 400) / 1920)
            scaledHeight = int((height * 450) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(scaledWidth / 2))
            y = y - (math.ceil(scaledHeight / 2))
            self.setGeometry(x, y, scaledWidth, scaledHeight)

            self.repaint()
            QtWidgets.QApplication.processEvents()

        except Exception as e:
            logger.critical("Error in scaleUI() in Off-Target.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #center UI on current screen
    def centerUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            # center window on current screen
            width = self.width()
            height = self.height()
            screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
            centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
            x = centerPoint.x()
            y = centerPoint.y()
            x = x - (math.ceil(width / 2))
            y = y - (math.ceil(height / 2))
            self.setGeometry(x, y, width, height)

            self.repaint()
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            logger.critical("Error in centerUI() in Off-Target.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

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
                    hold = gzip.open(file, 'r')
                    buf = (hold.readline())
                    hold.close()
                    buf = str(buf)
                    buf = buf.strip("'b")
                    buf = buf[:len(buf) - 2]
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
        except Exception as e:
            logger.critical("Error in fill_data_dropdown() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    def change_endos(self):
        try:
            #update file names based on current org/endo combo
            self.cspr_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][str(self.EndocomboBox.currentText())][0]
            self.db_file = self.organisms_to_files[str(self.OrgcomboBox.currentText())][str(self.EndocomboBox.currentText())][1]
        except Exception as e:
            logger.critical("Error in change_endos() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

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
            logger.critical("Error in update_endos() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

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
            #get tolerance value
            self.tolerance = self.toleranceSpinBox.value()

            #reset bools for new command to run
            self.perc = False
            self.bool_temp = False
            self.running = False

            if (self.AVG.isChecked()):
                avg_output = r'TRUE'
                detailed_output = r' FALSE '
            else:
                avg_output = r'FALSE'
                detailed_output = r' TRUE '

            #setup arguments for C++ .exe
            app_path = GlobalSettings.appdir.replace('\\','/')
            if platform.system() == 'Windows':
                exe_path = app_path + r'OffTargetFolder/OT_Win.exe'
            elif platform.system() == 'Linux':
                exe_path = app_path + r'OffTargetFolder/OT_Lin'
            else:
                exe_path = app_path + r'OffTargetFolder/OT_Mac'
            exe_path = '"' +  exe_path + '"'
            data_path = ' "' + app_path + 'OffTargetFolder/temp.txt' + '"' ##
            cspr_path = ' "' + GlobalSettings.CSPR_DB + '/' + self.cspr_file + '"'
            db_path = ' "' + GlobalSettings.CSPR_DB + '/' + self.db_file + '"'
            self.output_path = ' "' + GlobalSettings.CSPR_DB + '/' + self.FileName.text() + '"'
            filename = self.output_path
            filename = filename[:len(filename) - 1]
            filename = filename[1:]
            filename = filename.replace('"', '')
            exists = os.path.isfile(filename)
            CASPER_info_path = r' "' + app_path + 'CASPERinfo' + '" '
            num_of_mismathes = int(self.mismatchcomboBox.currentText())
            tolerance = self.tolerance
            endo = ' "' + self.EndocomboBox.currentText() + '"'
            hsu = ' "' + GlobalSettings.mainWindow.Results.endo_data[self.EndocomboBox.currentText()][2] + '"'

            #create command string
            cmd = exe_path + data_path + endo + cspr_path + db_path + self.output_path + CASPER_info_path + str(num_of_mismathes) + ' ' + str(tolerance) + detailed_output + avg_output + hsu
            if platform.system() == 'Windows':
                cmd = cmd.replace('/', '\\')

            #used to know when the process is done
            def finished():
                self.running = False
                self.progressBar.setValue(100)

            #used to know when data is ready to read from stdout
            def dataReady():
                #filter the data from stdout, bools used to know when the .exe starts outputting the progress
                #percentages to be able to type cast them as floats and update the progress bar. Also, must
                #split the input read based on '\n\ characters since the stdout read can read multiple lines at
                #once and is all read in as raw bytes
                line = str(self.process.readAllStandardOutput())

                line = line[2:]
                line = line[:len(line)-1]
                if platform.system() == 'Windows':
                    for lines in filter(None, line.split(r'\r\n')):
                        if(lines.find("Running Off Target Algorithm for") != -1 and self.perc == False):
                            self.perc = True
                        if(self.perc == True and self.bool_temp == False and lines.find("Running Off Target Algorithm for") == -1):
                            lines = lines[32:]
                            lines = lines.replace("%","")
                            if(float(lines) <= 99.5):
                                num = float(lines)
                                self.progressBar.setValue(num)
                            else:
                                self.bool_temp = True
                else:
                    for lines in filter(None, line.split(r'\n')):
                        if(lines.find("Running Off Target Algorithm for") != -1 and self.perc == False):
                            self.perc = True
                        if(self.perc == True and self.bool_temp == False and lines.find("Running Off Target Algorithm for") == -1):
                            lines = lines[32:]
                            lines = lines.replace("%","")
                            if(float(lines) <= 99.5):
                                num = float(lines)
                                self.progressBar.setValue(num)
                            else:
                                self.bool_temp = True


            #connect QProcess to the dataReady func, and finished func, reset progressBar only if the outputfile name
            #given does not already exist
            if(exists == False):
                self.process.readyReadStandardOutput.connect(partial(dataReady))
                self.process.readyReadStandardError.connect(partial(dataReady))
                self.progressBar.setValue(0)
                QtCore.QTimer.singleShot(100, partial(self.process.start, cmd))
                self.process.finished.connect(finished)

            else: #error message about file already being created
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Error")
                msgBox.setText("Output file already exists. Please choose a new output file name.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

        except Exception as e:
            logger.critical("Error in run_command() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #linked to run button
    def run_analysis(self):
        try:
            #make sure an analysis isn't already running before starting
            if(self.running == False):
                self.running = True
                self.run_command()
        except Exception as e:
            logger.critical("Error in run_analysis() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

    #exit linked to user clicking cancel, resets bools, and kills process if one was running
    def exit(self):
        try:
            self.perc = False
            self.bool_temp = False
            self.running = False
            self.process.kill()
            self.hide()
        except Exception as e:
            logger.critical("Error in exit() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

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
            logger.critical("Error in closeEvent() in OffTarget.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            exit(-1)

