from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import GlobalSettings
import matplotlib
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
from matplotlib.ticker import MaxNLocator
import os
import sqlite3
import gzip
from collections import Counter
import statistics
import traceback
import math
import numpy as np
from matplotlib.widgets import Slider
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

#global logger
logger = GlobalSettings.logger

class Multitargeting(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            self.count = 0
            super(Multitargeting, self).__init__()
            uic.loadUi(GlobalSettings.appdir + 'mt.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.multitargeting_statistics = Multitargeting_Statistics()

            self.sq = SeqTranslate()  # SeqTranslate object used in class
            self.line_bool = False # Used to check if VBoxLayout already has canvas in it
            self.bar_bool = False # Used to check if VBoxLayout already has canvas in it
            self.seed_bar_bool = False # Used to check if VBoxLayout already has canvas in it

            # GroupBox Styling
            groupbox_style = """
            QGroupBox:title{subcontrol-origin: margin;
                            left: 10px;
                            padding: 0 5px 0 15px;}
            QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                            border-radius: 9px;
                            font: bold;
                            margin-top: 10px;}"""
            self.groupBox.setStyleSheet(groupbox_style)
            self.groupBox_2.setStyleSheet(groupbox_style.replace("groupBox","groupBox_2"))
            self.groupBox_3.setStyleSheet(groupbox_style.replace("groupBox","groupBox_3"))

            #layout for table
            self.table.setColumnCount(8)
            self.table.setShowGrid(False)
            # for keeping track of where we are in the sorting clicking for each column
            self.switcher_table = [1, 1, 1, 1, 1, 1, 1, 1]
            self.table.setHorizontalHeaderLabels(
                ["Seed", "Total Repeats", "Avg. Repeats/Scaffold", "Consensus Sequence", "% Consensus",
                 "Score", "PAM", "Strand"])
            self.table.horizontalHeader().setSectionsClickable(True)
            self.table.horizontalHeader().sectionClicked.connect(self.table_sorting)
            self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.table.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
            self.table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.table.resizeColumnsToContents()

            # Initializes layouts for the graphs
            self.global_line = QtWidgets.QVBoxLayout()
            self.global_bar = QtWidgets.QVBoxLayout()
            self.seed_bar = QtWidgets.QVBoxLayout()
            self.global_line.setContentsMargins(0,0,0,0)
            self.global_bar.setContentsMargins(0,0,0,0)
            self.seed_bar.setContentsMargins(0,0,0,0)
            self.chrom_viewer_layout = QtWidgets.QVBoxLayout()
            self.chrom_viewer_layout.setContentsMargins(0,0,0,0)

            self.data = ""
            self.shortHand = ""
            self.chromo_length = list()

            # Listeners for changing the seed sequence or the .cspr file
            self.Analyze_Button.clicked.connect(self.make_graphs)
            self.statistics_overview.clicked.connect(self.show_statistics)
            self.export_button.clicked.connect(self.export_tool)
            self.selectAll.stateChanged.connect(self.select_all)
            self.selectAll.setEnabled(False)

            # go back to main button
            self.back_button.clicked.connect(self.go_back)

            # Statistics storage variables
            self.max_repeats = 1
            self.average = 0
            self.median = 0
            self.mode = 0
            self.average_unique = 0
            self.average_rep = 0
            self.bar_coords = []
            self.seed_id_seq_pair = {}

            # parser object
            self.parser = CSPRparser("")

            self.ready_chromo_min_max = True
            self.ready_chromo_make_graph = True
            self.info_path = os.getcwd()

            ##################################
            self.scene = QtWidgets.QGraphicsScene()
            self.scene2 = QtWidgets.QGraphicsScene()
            self.graphicsView_2.setScene(self.scene2)
            self.scrollArea.viewport().installEventFilter(self)
            self.graphicsView_2.viewport().installEventFilter(self)

            self.loading_window = loading_window()

            #sql query settings
            self.row_limit = 1000
            self.sql_query_settings.clicked.connect(self.update_sql_query_settings)
            self.sql_settings = sql_query_settings()
            self.sql_settings.row_count.textChanged.connect(self.sql_row_count_value_changed)

            #scale UI
            self.first_show = True
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing Multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()
            exit(-1)

    def select_all(self):
        try:
            if self.selectAll.isChecked():
                self.table.selectAll()
            else:
                self.table.clearSelection()
        except Exception as e:
            logger.critical("Error in selectAll() in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()


            exit(-1)


    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            self.dpi = dpi
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            fontSize = 12
            self.fontSize = fontSize
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )
            self.menuBar().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';" )

            #CASPER header scaling
            fontSize = 30
            self.title.setStyleSheet("font: bold " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            #make sure chromosome viewer doesnt get too small
            self.groupBox_2.setMinimumHeight(int(0.3 * height))

            # window scaling
            scaledWidth = int((width * 1400) / 1920)
            scaledHeight = int((height * 900) / 1080)

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            #set min width of table
            self.table.setMinimumWidth(int(0.5 * scaledWidth))

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
            logger.critical("Error in scaleUI() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def centerUI(self):
        self.repaint()
        QtWidgets.QApplication.processEvents()

        width = self.width()
        height = self.height()
        # scale/center window
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        x = centerPoint.x()
        y = centerPoint.y()
        x = x - (math.ceil(width / 2))
        y = y - (math.ceil(height / 2))
        self.setGeometry(x, y, width, height)

        self.Analyze_Button.resize(200, 200)

        self.repaint()
        QtWidgets.QApplication.processEvents()

    def export_tool(self):
        try:
            select_items = self.table.selectedItems()
            if len(select_items) <= 0:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Nothing Selected")
                msgBox.setText("No targets were highlighted. Please highlight the targets you want to be exported to a CSV File!")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return
            GlobalSettings.mainWindow.export_tool_window.launch(select_items,"mt")
        except Exception as e:
            logger.critical("Error in export_tool() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def show_statistics(self):
        try:
            if (self.line_bool and self.bar_bool):
                self.multitargeting_statistics.centerUI()
                self.multitargeting_statistics.show()
                self.multitargeting_statistics.activateWindow()
            else:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("No analysis run.")
                msgBox.setText('Multitargeting Analysis must be performed before viewing statistics.\n\nSelect an organism and endonuclease and click "Analyze" then try again.')
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return True
        except Exception as e:
            logger.critical("Error in show_statistics() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #event handler to show details of targets in chromosome viewer while hovering over canvases
    def chromosome_event_handler(self, event):
        try:
            #get current moust location
            x = event.xdata
            y = event.y

            #get event data relative to the canvas (chromosome) the viewer is hovering at
            curr_chromosome = self.canvas_chromosome_map[event.canvas]
            chromosome_seed_data = self.event_data[curr_chromosome]

            #get targets within small range of the mouse location
            local_targets = []
            for entry in chromosome_seed_data:
                try:
                    if x >= entry[0] - 0.001 and x <= entry[0] + 0.001:
                        local_targets.append(entry)
                except:
                    pass

            #make sure targets are found before overwriting the viewers details
            if local_targets != []:
                #prep the viewer to show the target details
                self.scene2 = QtWidgets.QGraphicsScene()
                self.graphicsView_2.setScene(self.scene2)
                output = str()
                for target in local_targets:
                    output += f"Location: {target[1]} | Seq: {target[2]} | PAM: {target[3]} | SCR: {target[4]} | DIRA: {target[5]}\n"

                text = self.scene2.addText(output)
                font = QtGui.QFont()
                font.setPointSize(self.fontSize-2)
                text.setFont(font)

        except Exception as e:
            logger.critical("Error in event_data() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def launch(self):
        try:
            self.get_data()
        except Exception as e:
            logger.critical("Error in launch() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #button trigger for sql settings
    def update_sql_query_settings(self):
        try:
            self.sql_settings.centerUI()
            self.sql_settings.show()
            self.sql_settings.activateWindow()
        except Exception as e:
            logger.critical("Error in update_sql_query_settings() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #trigger for if sql line edit value has changed
    def sql_row_count_value_changed(self):
        try:
            self.row_limit = int(self.sql_settings.row_count.text())
        except Exception as e:
            logger.critical("Error in sql_row_count_value_changed() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            pass

    def get_data(self):
        try:
            #disconnect index changed signal on endo dropdown if there is one
            try:
                self.organism_drop.currentIndexChanged.disconnect()
            except:
                pass
            try:
                self.endo_drop.currentIndexChanged.disconnect()
            except:
                pass


            self.organism_drop.clear()
            self.endo_drop.clear()

            onlyfiles = [f for f in os.listdir(GlobalSettings.CSPR_DB) if os.path.isfile(os.path.join(GlobalSettings.CSPR_DB, f))]
            self.organisms_to_files = {}
            self.organisms_to_endos = {}
            # shortName = {}
            self.endo_drop.clear()
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
                        if self.organism_drop.findText(species) == -1:
                            self.organism_drop.addItem(species)

            #fill in endos dropdown based on current organism
            endos = self.organisms_to_endos[str(self.organism_drop.currentText())]
            self.endo_drop.addItems(endos)
            self.organism_drop.currentIndexChanged.connect(self.update_endos)

            #update file names for current org/endo combo
            self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][endos[0]][0]
            self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][endos[0]][1]
        except Exception as e:
            logger.critical("Error in get_data() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def update_endos(self):
        try:
            #try to disconnect index changed signal on endo dropdown if there is one
            try:
                self.endo_drop.currentIndexChanged.disconnect()
            except:
                pass

            #clear endo dropdown and fill in with endos relative to the current organism
            self.endo_drop.clear()
            endos = self.organisms_to_endos[str(self.organism_drop.currentText())]
            self.endo_drop.addItems(endos)
        except Exception as e:
            logger.critical("Error in update_endos() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def make_graphs(self):
        try:
            self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][str(self.endo_drop.currentText())][0]
            self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][str(self.endo_drop.currentText())][1]

            self.loading_window.loading_bar.setValue(0)
            self.loading_window.centerUI()
            self.loading_window.show()
            QtCore.QCoreApplication.processEvents()
            self.chromo_length.clear()
            self.loading_window.loading_bar.setValue(10)
            self.plot_repeats_vs_seeds()
            self.loading_window.loading_bar.setValue(30)
            self.bar_seeds_vs_repeats()
            self.loading_window.loading_bar.setValue(50)
            self.fill_table()
            self.table.selectRow(0) # Set the index to first row by default so
            # all graphs are generated
            self.selectAll.setEnabled(True) # Enable select all checkbox only
            # after analysis has been performed
            self.loading_window.loading_bar.setValue(100)
            self.multitargeting_statistics.avg_rep.setText(str(round(float(self.average),1)))
            self.multitargeting_statistics.med_rep.setText(str(round(float(self.median),1)))
            self.multitargeting_statistics.mode_rep.setText(str(round(float(self.mode),1)))
            self.multitargeting_statistics.nbr_seq.setText(str(round(float(self.repeat_count),1)))
            self.loading_window.hide()
            self.repaint()
            QtWidgets.QApplication.processEvents()

        except Exception as e:
            logger.critical("Error in make_graphs() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #function to fill table in UI
    def fill_table(self):
        try:
            #disable row triggers
            try:
                self.table.itemSelectionChanged.disconnect()
            except:
                pass

            #empty table
            self.table.setRowCount(0)

            #query db file for data
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            row_cnt = 0

            if self.row_limit == -1:
                sql_query = "SELECT * FROM repeats ORDER BY count DESC;"
            else:
                sql_query = "SELECT * FROM repeats ORDER BY count DESC LIMIT 0, " + str(self.row_limit) + ";"

            for repeat in c.execute(sql_query):
                #expand table by 1 row
                self.table.setRowCount(row_cnt + 1)

                #extract repeat info
                seed = repeat[0]
                chroms = repeat[1].split(",")
                locs = repeat[2].split(",")
                threes = repeat[3].split(",")
                fives = repeat[4].split(",")
                pams = repeat[5].split(",")
                scores = repeat[6].split(",")
                count = repeat[7]

                #push seed
                table_seed = QtWidgets.QTableWidgetItem()
                table_seed.setData(QtCore.Qt.EditRole, seed)
                table_seed.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 0, table_seed)

                if len(threes) < len(fives):
                    for i in range(len(fives) - len(threes)):
                        threes.append('')

                elif len(fives) < len(threes):
                    for i in range(len(threes) - len(fives)):
                        fives.append('')

                majority_index = 0
                three_prime, five_prime, both_prime = False, False, False
                if threes[0] == '':
                    majority = max(set(fives), key=fives.count)
                    majority_index = fives.index(majority)
                    five_prime = True
                elif fives[0] == '':
                    majority = max(set(threes), key=threes.count)
                    majority_index = threes.index(majority)
                    three_prime = True
                else:
                    # account for both 3 and 5 present
                    threes_and_fives = []
                    for i in range(len(threes)):
                        threes_and_fives.append(threes[i] + fives[i])
                    majority = max(set(threes_and_fives), key=threes_and_fives.count)
                    majority_index = threes_and_fives.index(majority)
                    both_prime = True

                # push total count
                table_count = QtWidgets.QTableWidgetItem()
                table_count.setData(QtCore.Qt.EditRole, count)
                table_count.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 1, table_count)

                # push avg repeat
                location_repeat_counts = Counter(chroms)
                avg_rep_per_scaff = sum(location_repeat_counts.values()) / len(location_repeat_counts.values())
                avg_rep = QtWidgets.QTableWidgetItem()
                avg_rep_per_scaff = float("%.2f" % avg_rep_per_scaff)
                avg_rep.setData(QtCore.Qt.EditRole, avg_rep_per_scaff)
                avg_rep.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 2, avg_rep)

                # push seq
                seq = QtWidgets.QTableWidgetItem()
                seq.setData(QtCore.Qt.EditRole, fives[majority_index] + seed + threes[majority_index])
                seq.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 3, seq)

                # push percent consensus
                perc_cons = QtWidgets.QTableWidgetItem()
                percent_consensus = 0
                if five_prime == True:
                    percent_consensus = (fives.count(fives[majority_index]) / len(fives)) * 100
                elif three_prime == True:
                    percent_consensus = (threes.count(threes[majority_index]) / len(threes)) * 100
                elif both_prime:
                    percent_consensus = (threes_and_fives.count(threes_and_fives[majority_index]) / len(threes_and_fives)) * 100

                percent_consensus = float("%.2f" % percent_consensus)
                perc_cons.setData(QtCore.Qt.EditRole, percent_consensus)
                perc_cons.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 4, perc_cons)

                # push score
                score = QtWidgets.QTableWidgetItem()
                score.setData(QtCore.Qt.EditRole, scores[majority_index])
                score.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 5, score)

                # push PAM
                pam = QtWidgets.QTableWidgetItem()
                pam.setData(QtCore.Qt.EditRole, pams[majority_index])
                pam.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 6, pam)

                # push strand
                strand_val = ""
                if int(locs[majority_index]) < 0:
                    strand_val = "-"
                else:
                    strand_val = "+"

                strand = QtWidgets.QTableWidgetItem()
                strand.setData(QtCore.Qt.EditRole, strand_val)
                strand.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table.setItem(row_cnt, 7, strand)

                #increment row count
                row_cnt += 1

            #close db connection
            c.close()
            conn.close()

            #resize columns in table to match contents
            self.table.resizeColumnsToContents()

            #reconnect row trigger
            self.table.itemSelectionChanged.connect(self.row_selection_trigger)
        except Exception as e:
            logger.critical("Error in fill_table() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #function for triggering graph updates when user selects row in table
    def row_selection_trigger(self):
        try:
            item = self.table.currentItem()
            if item.isSelected() == True:
                row_num = item.row()
                seed = self.table.item(row_num,0).text()
                self.fill_Chromo_Text(seed)
                self.chro_bar_create(seed)
        except Exception as e:
            logger.critical("Error in row_selection_trigger() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # sorting to table
    def table_sorting(self, logicalIndex):
        try:
            self.switcher_table[logicalIndex] *= -1
            if self.switcher_table[logicalIndex] == -1:
                self.table.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
            else:
                self.table.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)
        except Exception as e:
            logger.critical("Error in table_sorting() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #fill in chromo bar visualization
    def fill_Chromo_Text(self, seed):
        try:
            #global dictionary to map canvases to chromosomes
            self.canvas_chromosome_map = {}

            # get kstats
            kstats = []
            with open(self.cspr_file, "r") as f:
                for line in f:
                    buf = str(line)
                    if buf.find("KARYSTATS") != -1:
                        buf = buf.replace("KARYSTATS: ", "")
                        kstats = buf.split(',')
                        kstats = kstats[:-1]
                        break

            #get chromosomes/locations of repeats for current seed
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            data = c.execute("SELECT chromosome, location, pam, score, five, three FROM repeats WHERE seed = ? ", (seed,)).fetchone()
            c.close()

            #make sure there is data on current seed
            if data != None:

                # build out dictionary (seed_data) mapping chromosome numbers to list of repeat locations
                # normalize the location values based on kstat value of chromosome
                seed_data = {}
                self.event_data = {}
                data = list(data)
                chromo = data[0].split(',')
                pos = data[1].split(',')
                pam = data[2].split(',')
                score = data[3].split(',')
                five = data[4].split(',')
                three = data[5].split(',')
                #get x-cordinates for graphs and event data
                for i in range(len(chromo)):
                    curr_chromo = int(chromo[i])
                    if int(pos[i]) >= 0:
                        dir = "+"
                    else:
                        dir = "-"
                    normalized_location = abs(float(pos[i]) / float(kstats[curr_chromo - 1]))
                    if int(chromo[i]) in seed_data.keys():
                        seed_data[int(chromo[i])].append(normalized_location)
                        self.event_data[int(chromo[i])].append([normalized_location, pos[i], five[i] + seed + three[i], pam[i], score[i], dir])
                    else:
                        seed_data[int(chromo[i])] = [normalized_location]
                        self.event_data[int(chromo[i])] = [[normalized_location, pos[i], five[i] + seed + three[i] ,pam[i], score[i], dir]]

                # graph the locations for each chromosome
                # Clear out old widgets in layout
                for i in reversed(range(self.chrom_viewer_layout.count())):
                    self.chrom_viewer_layout.itemAt(i).widget().setParent(None)

                top_widget = QtWidgets.QWidget()
                top_layout = QtWidgets.QVBoxLayout()

                chromo_keys = sorted(list(seed_data.keys()))

                screen = self.screen()
                height = screen.geometry().height()
                groupbox_height = int((height * 100) / 1080)

                for i in range(len(chromo_keys)):
                    curr_chromo = chromo_keys[i]
                    group_box = QtWidgets.QGroupBox()

                    group_box.setTitle(f"Chromosome {curr_chromo}")
                    group_box.setMinimumHeight(groupbox_height)
                    group_box.setMaximumHeight(groupbox_height)
                    layout = QtWidgets.QVBoxLayout(group_box)

                    canvas = MplCanvas()
                    canvas.axes.eventplot(seed_data[curr_chromo])
                    canvas.mpl_connect("motion_notify_event", self.chromosome_event_handler)

                    # surrounding border
                    canvas.axes.hlines(1.5, -0.01, 1.01, colors="Black", linewidth=1.5)
                    canvas.axes.hlines(0.5, -0.01, 1.01, colors="Black", linewidth=1.5)
                    canvas.axes.vlines(-0.01, 0.5, 1.5, colors="Black", linewidth=1.5)
                    canvas.axes.vlines(1.01, 0.5, 1.5, colors="Black", linewidth=1.5)
                    canvas.axes.set_ylim(0.45, 1.55)
                    canvas.axes.set_xlim(-0.05, 1.05)

                    canvas.axes.axis('off')
                    canvas.draw()

                    self.canvas_chromosome_map[canvas] = curr_chromo

                    layout.addWidget(canvas)
                    top_layout.addWidget(group_box)

                top_widget.setLayout(top_layout)
                self.scrollArea.setWidget(top_widget)

                return True
            else:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
                msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
                msgBox.setWindowTitle("Seed Error")
                msgBox.setText("No such seed exists in the repeats section of this organism.")
                msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
                msgBox.exec()

                return False
        except Exception as e:
            logger.critical("Error in fill_Chromo_text() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # creates bar graph num of repeats vs. chromosome
    def chro_bar_create(self, seed):
        try:
            ###Clear out old widgets in layout
            for i in reversed(range(self.seed_bar.count())):
                self.seed_bar.itemAt(i).widget().setParent(None)
            self.seed_canvas = MplCanvas(self, width=5, height=3, dpi=self.dpi) ###Initialize new Canvas
            self.seed_bar.addWidget(self.seed_canvas) ### Add canvas to global line layout
            self.repeats_vs_chromo.setLayout(self.seed_bar) ### Add global line layout to repeats vs. seeds line plot widget
            y = []
            x_labels = []
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            data = c.execute("SELECT chromosome FROM repeats WHERE seed = ? ", (seed,)).fetchone()
            c.close()
            data = data[0].split(',')
            for i in range(len(data)):
                data[i] = int(data[i])
            bar_data = Counter(data)
            for chromo in sorted(bar_data):
                x_labels.append(chromo)
                y.append(bar_data[chromo])
            x = list(range(0, len(x_labels)))

            #the following statements are plottings / formatting for the graph
            self.seed_canvas.axes.bar(x, y, align='center')
            self.seed_canvas.axes.yaxis.set_major_locator(MaxNLocator(integer=True))
            self.seed_canvas.axes.set_ylim(0, max(y) + 1)
            self.seed_canvas.axes.set_xticks(x)
            self.seed_canvas.axes.set_xticklabels(x_labels)
            if len(x_labels) > 10:
                tick_spacing = round(len(x_labels)/10)
                for i, t in enumerate(self.seed_canvas.axes.get_xticklabels()):
                    if (i % tick_spacing) != 0:
                        t.set_visible(False)
            self.seed_canvas.axes.set_xlabel('Chromosome', fontsize = 10)
            self.seed_canvas.axes.set_ylabel('Number of Repeats', fontsize=10)
            self.line_canvas.axes.set_title('Repeats per ID Number',fontsize=10)
            self.line_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
            self.line_canvas.draw()
        except Exception as e:
            logger.critical("Error in chro_bar_create() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def bar_seeds_vs_repeats(self):
        try:
            ###Clear out old widgets in layout
            for i in reversed(range(self.global_bar.count())):
                self.global_bar.itemAt(i).widget().setParent(None)
            self.bar_canvas = MplCanvas(self, width=5, height=3, dpi=self.dpi) ###Initialize new Canvas
            self.global_bar.addWidget(self.bar_canvas) ### Add canvas to global line layout
            self.seeds_vs_repeats_bar.setLayout(self.global_bar) ### Add global line layout to repeats vs. seeds line plot widget

            """ Get the data """
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            x_labels = []
            y = []
            for obj in c.execute("select count, COUNT(count) as cnt from repeats group by count order by cnt DESC;"):
                x_labels.append(obj[0])
                y.append(obj[1])

            x_labels = sorted(x_labels)

            # the following are plotting / formatting for the graph
            self.bar_canvas.axes.scatter(x_labels, y, s=10)
            self.bar_canvas.axes.set_xlim(x_labels[0]-0.5, x_labels[-1] + 0.5)
            self.bar_canvas.axes.set_yscale('log')
            self.bar_canvas.axes.set_xlabel('Number of Repeats', fontsize=10)
            self.bar_canvas.axes.set_ylabel('Number of Sequences', fontsize=10)
            self.bar_canvas.axes.set_title('Number of Sequences per Number of Repeats',fontsize=10)
            self.bar_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
            self.bar_canvas.draw()

            self.bar_bool = True
        except Exception as e:
            logger.critical("Error in bar_seeds_vs_repeats() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # # plots the repeats per ID number graph as line graph
    def plot_repeats_vs_seeds(self):
        try:
            ###Clear out old widgets in layout
            for i in reversed(range(self.global_line.count())):
                self.global_line.itemAt(i).widget().setParent(None)
            self.line_canvas = MplCanvas(self, width=5, height=3, dpi=self.dpi) ###Initialize new Canvas
            self.global_line.addWidget(self.line_canvas) ### Add canvas to global line layout
            self.repeats_vs_seeds_line.setLayout(self.global_line) ### Add global line layout to repeats vs. seeds line plot widget

            """ Fetch all the data """
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()

            y1 = []
            for obj in c.execute("SELECT count from repeats;"):
                y1.append(obj[0])

            c.close()

            #get stats
            self.average = statistics.mean(y1)
            self.mode = statistics.mode(y1)
            self.median = statistics.median(y1)
            self.repeat_count = len(y1)

            # clear axes
            self.line_canvas.axes.clear()
            #Plotting / formatting
            self.line_canvas.axes.plot(y1)
            self.line_canvas.axes.set_xlabel('Seed ID Number',fontsize=10)
            self.line_canvas.axes.set_ylabel('Number of Repeats',fontsize=10)
            self.line_canvas.axes.set_title('Number of Repeats per Seed ID Number',fontsize=10)
            self.line_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
            # always redraw at the end
            self.line_bool = True
            self.line_canvas.draw()
        except Exception as e:
            logger.critical("Error in plot_repeats_vs_seeds() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def seed_chromo_changed(self):
        try:
            self.loading_window.loading_bar.setValue(5)
            self.loading_window.show()
            QtCore.QCoreApplication.processEvents()
            self.seed.setText('')
            self.loading_window.loading_bar.setValue(50)
            self.loading_window.loading_bar.setValue(100)
            self.loading_window.hide()
        except Exception as e:
            logger.critical("Error in seed_chromo_changed() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    #connects to go back button in bottom left to switch back to the main CASPER window
    def go_back(self):
        try:
            self.sql_settings.hide()
            self.multitargeting_statistics.hide()
            GlobalSettings.mainWindow.show()
            self.hide()
        except Exception as e:
            logger.critical("Error in go_back() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        try:
            GlobalSettings.mainWindow.closeFunction()
            self.sql_settings.hide()
            self.multitargeting_statistics.hide()
            event.accept()
        except Exception as e:
            logger.critical("Error in closeEvent() in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)


class loading_window(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(loading_window, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "loading_data_form.ui", self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.loading_bar.setValue(0)
            self.setWindowTitle("Loading Data")

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing loading_window class in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    def scaleUI(self):
        try:
            self.repaint()
            QtWidgets.QApplication.processEvents()

            screen = self.screen()
            dpi = screen.physicalDotsPerInch()
            width = screen.geometry().width()
            height = screen.geometry().height()

            # font scaling
            # 14px is used for 92 dpi
            fontSize = 12
            self.centralWidget().setStyleSheet("font: " + str(fontSize) + "pt 'Arial';")

            self.adjustSize()

            currentWidth = self.size().width()
            currentHeight = self.size().height()

            # scale/center window
            scaledWidth = int((width * 450) / 1920)
            scaledHeight = int((height * 125) / 1080)

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
            logger.critical("Error in scaleUI() in loading_window() class in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

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
            logger.critical("Error in centerUI() in loading window in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        try:
            # Create Figure instance
            #fig = Figure(figsize=(width, height), dpi=dpi)

            fig = Figure(dpi=dpi, tight_layout=True)
            self.axes = fig.add_subplot(111)
            self.axes.clear()
            super(MplCanvas, self).__init__(fig)

            # Create axes
            "this is the one need to be changed and be dynamic for any number of subplots"
            #self.axes = fig.add_subplot(111)

            # Display the figure
            # FigureCanvasQTAgg.__init__(self, fig)
            # self.setParent(parent)
            # FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
            #                                 QtWidgets.QSizePolicy.Expanding)
            # FigureCanvasQTAgg.updateGeometry(self)
        except Exception as e:
            logger.critical("Error initializing MplCanvas class in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)


class Multitargeting_Statistics(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        try:
            super(Multitargeting_Statistics, self).__init__(parent)
            uic.loadUi(GlobalSettings.appdir + 'multitargeting_stats.ui', self)
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Statistics")

            #scale UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing Multitargeting_Statistics class in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

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
            # 16px is used for 92 dpi / 1920x1080
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
            scaledWidth = int((width * 275) / 1920)
            scaledHeight = int((height * 185) / 1080)

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
            logger.critical("Error in scaleUI() in multitargeting statistics in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)

    # center UI on current screen
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
            logger.critical("Error in centerUI() in multitargeting statistics in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)


class sql_query_settings(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(sql_query_settings, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "multitargeting_sql_settings.ui", self)
            self.setWindowTitle("SQL Settings")
            self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.row_count.setValidator(QtGui.QIntValidator())

            #scale the UI
            self.scaleUI()

        except Exception as e:
            logger.critical("Error initializing sql_query_settings class in multi-targeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

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
            # 16px is used for 92 dpi / 1920x1080
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
            scaledWidth = int((width * 375) / 1920)
            scaledHeight = int((height * 140) / 1080)

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
            logger.critical("Error in scaleUI() in sql settings in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

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
            logger.critical("Error in centerUI() in sql settings in multitargeting.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            msgBox = QtWidgets.QMessageBox()
            msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            msgBox.setWindowTitle("Fatal Error")
            msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
            msgBox.exec()

            exit(-1)
