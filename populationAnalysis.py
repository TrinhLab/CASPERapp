from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import GlobalSettings
import os
from CSPRparser import CSPRparser
import Algorithms
from functools import partial
import numpy as np
from PyQt5.QtWidgets import *
import gzip
import sqlite3
from collections import Counter
import statistics
import itertools
from matplotlib import cm
import matplotlib
import mplcursors
import time
import copy


class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'populationanalysis.ui', self)
        self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
        self.goBackButton.clicked.connect(self.go_back)
        self.analyze_button.clicked.connect(self.pre_analyze)
        self.clear_Button.clicked.connect(self.clear)
        self.export_button.clicked.connect(self.export_to_csv)
        self.sq = Algorithms.SeqTranslate()
        self.ref_para_list = list()
        self.mode = 0
        self.find_locs_button.clicked.connect(self.find_locations)
        self.clear_loc_button.clicked.connect(self.clear_loc_table)
        self.directory = ""
        self.names = []
        self.names_venn = []

        """ Colormap Graph initialization """
        self.colormap_layout = QtWidgets.QVBoxLayout()
        self.colormap_layout.setContentsMargins(0,0,0,0)
        self.colormap_canvas = MplCanvas(self, width=5, height=4, dpi=100) ###Initialize new Canvas


        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 15px;}
        QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        font: 15pt "Arial";
                        font: bold;
                        margin-top: 10px;}"""
        self.groupBox.setStyleSheet(groupbox_style)
        self.groupBox_2.setStyleSheet(groupbox_style.replace("groupBox","groupBox_2"))

        #organism table
        self.org_Table.setColumnCount(1)
        self.org_Table.setShowGrid(False)
        self.org_Table.setHorizontalHeaderLabels(["Organism"])
        self.org_Table.horizontalHeader().setSectionsClickable(True)
        self.org_Table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.org_Table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.org_Table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.org_Table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.org_Table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)

        #top-right table
        self.table2.setColumnCount(9)
        self.table2.setShowGrid(False)
        self.table2.setHorizontalHeaderLabels(["Seed","% Coverage","Total Repeats","Avg. Repeats/Scaffold", "Consensus Sequence", "% Consensus", "Score","PAM", "Strand"])
        self.table2.horizontalHeader().setSectionsClickable(True)
        self.table2.horizontalHeader().sectionClicked.connect(self.table2_sorting)
        self.table2.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table2.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.table2.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.table2.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.table2.resizeColumnsToContents()

        # Finder table
        self.loc_finder_table.setColumnCount(5)
        self.loc_finder_table.setShowGrid(False)
        self.loc_finder_table.setHorizontalHeaderLabels(["Seed ID", "Sequence", "Organism", "Scaffold", "Location"])
        self.loc_finder_table.horizontalHeader().setSectionsClickable(True)
        self.loc_finder_table.horizontalHeader().sectionClicked.connect(self.loc_table_sorter)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch) #This keeps the organism column from being too small.
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.loc_finder_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.loc_finder_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
#        self.loc_finder_table.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
#        self.loc_finder_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.loc_finder_table.resizeColumnsToContents()

        #custom seed search
        self.query_seed_button.clicked.connect(self.custom_seed_search)

        self.total_org_number = 0

        self.switcher_table2 = [1, 1, 1, 1, 1, 1, 1, 1,
                                1]  # for keeping track of where we are in the sorting clicking for each column
        self.switcher_loc_table = [1, 1, 1, 1, 1]

        #Window centering stuff
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.mwfg.moveCenter(self.cp)  ##Center window
        self.move(self.mwfg.topLeft())  ##Center window
        #screen = QtGui.QGuiApplication.screenAt(QtGui.QCursor().pos())


        #Initialize variables
        self.index_to_cspr = {}
        self.index_to_db = {}
        self.name_to_db = {}
        self.cspr_files = []
        self.db_files = []
        self.Endos = {}
        self.seeds = []

        self.loading_window = loading_window()

    def export_to_csv(self):
        select_items = self.table2.selectedItems()
        if len(select_items) <= 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "No targets were highlighted."
                                           "Please highlight the targets you want to be exported to a CSV File!",
                                           QtWidgets.QMessageBox.Ok)
            return
        GlobalSettings.mainWindow.export_csv_window.launch(select_items,9)

    def prevent_toggle(self):
        self.meta_genomic_cspr_checkbox.setChecked(QtCore.Qt.Checked)


    # this function calls the popParser function and fills all the tables
    def pre_analyze(self):

        #clear saved filenames
        self.cspr_files = []
        self.db_files = []
        self.rows = []

        #get selected indexes
        selected_indexes = self.org_Table.selectionModel().selectedRows()

        # error check
        if len(selected_indexes) == 0:
            QtWidgets.QMessageBox.question(self, "Error",
                                           "Please select CSPR file(s) for analysis.",
                                           QtWidgets.QMessageBox.Ok)
            return

        #get cspr and db filenames
        for index in sorted(selected_indexes):
            self.rows.append(index.row() + 1)
            self.cspr_files.append(self.index_to_cspr[index.row()])
            self.db_files.append(self.index_to_db[index.row()])

        self.get_org_names()
        self.fill_data()


    def launch(self, path):
        os.chdir(path)
        self.directory = path
        self.get_data()


    def get_data(self):
        if self.directory == '':
            return

        self.fillEndo()


    # this function opens CASPERinfo and builds the dropdown menu of selectable endonucleases
    def fillEndo(self):
        try:
            self.endoBox.currentIndexChanged.disconnect()
        except:
            pass

        self.Endos = {}
        self.endoBox.clear()
        f = open(GlobalSettings.appdir + 'CASPERinfo')
        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if (line[0] == "-"):
                        break
                    line_tokened = line.split(";")
                    endo = line_tokened[0]
                    # Checking to see if there is more than one pam sequence in the list
                    if line_tokened[1].find(",") != -1:
                        p_pam = line_tokened[1].split(",")[0]
                    else:
                        p_pam = line_tokened[1]
                    default_five_length = line_tokened[2]
                    default_seed_length = line_tokened[3]
                    default_three_length = line_tokened[4]
                    self.Endos[endo + " PAM: " + p_pam] = (endo, p_pam, default_five_length, default_seed_length, default_three_length)

                break
        f.close()
        self.endoBox.addItems(self.Endos.keys())
        self.endoBox.currentIndexChanged.connect(self.change_endo)
        self.change_endo()

    def change_endo(self):
        onlyfiles = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        index = 0
        for file in onlyfiles:
            if file.find('.cspr') != -1:
                endo = file[file.rfind('_') + 1:file.find('.cspr')]
                if endo == self.Endos[self.endoBox.currentText()][0]:

                    # increase row count
                    self.org_Table.setRowCount(index + 1)

                    # open .cspr file and get genome name
                    f = gzip.open(file, 'r')
                    line = f.readline()
                    f.close()
                    line = str(line)
                    line = line[:len(line) - 3]
                    line = line[2:]
                    colonIndex = line.find(':') + 2
                    orgName = line[colonIndex:]

                    # add genome name to table
                    tabWidget = QtWidgets.QTableWidgetItem(orgName)
                    tabWidget.setTextAlignment(QtCore.Qt.AlignVCenter)
                    self.org_Table.setItem(index, 0, tabWidget)

                    # store cspr and db file information for later
                    self.index_to_cspr[index] = file
                    db_file = file.replace(".cspr", "")
                    db_file += "_repeats.db"
                    self.index_to_db[index] = db_file

                    # incrase row index
                    index += 1

        if index == 0:
            self.org_Table.clearContents()
            self.org_Table.setRowCount(0)

        self.org_Table.resizeColumnsToContents()


    def fill_data(self):
        #update progress bar
        self.loading_window.loading_bar.setValue(5)
        self.loading_window.show()
        QtCore.QCoreApplication.processEvents()

        #prep table
        self.total_org_number = len(self.cspr_files)
        self.table2.setRowCount(0)
        self.loading_window.loading_bar.setValue(10)
        index = 0

        self.seeds = self.get_shared_seeds(self.db_files, True)

        try:
            os.remove(GlobalSettings.appdir + "temp_join.db")
        except:
            pass

        no_seeds = False

        if(len(self.seeds) == 0):
            no_seeds = True

        #QtCore.QCoreApplication.processEvents()

        #retrieve data on shared seeds
        if no_seeds == False:
            increase_val = float(15 / len(self.seeds))
            running_val = self.loading_window.loading_bar.value()
            self.loading_window.info_label.setText("Parsing Seed Data")
            self.counts = []
            for seed in self.seeds:
                # increase row count
                self.table2.setRowCount(index + 1)

                # push seed to table
                table_seed = QtWidgets.QTableWidgetItem()
                table_seed.setData(QtCore.Qt.EditRole, seed)
                table_seed.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 0, table_seed)

                total_count = 0
                org_count = 0
                threes = []
                fives = []
                scores = []
                pams = []
                locs = []

                for db_file in self.db_files:
                    conn = sqlite3.connect(db_file)
                    c = conn.cursor()
                    data = c.execute("SELECT count, three, five, pam, score, location FROM repeats WHERE seed = ? ",(seed,)).fetchone()
                    if data != None:
                        data = list(data)
                        org_count += 1
                        total_count += int(data[0])
                        threes += data[1].split(",")
                        fives += data[2].split(",")
                        pams += data[3].split(",")
                        scores += data[4].split(",")
                        locs += data[5].split(",")

                self.counts.append(total_count)

                if len(threes) < len(fives):
                    for i in range(len(fives) - len(threes)):
                        threes.append('')

                elif len(fives) < len(threes):
                    for i in range(len(threes) - len(fives)):
                        fives.append('')

                majority_index = 0
                if threes[0] == '':
                    majority = max(set(fives), key=fives.count)
                    majority_index = fives.index(majority)
                else:
                    majority = max(set(threes), key=threes.count)
                    majority_index = threes.index(majority)

                # push percent coverage
                perc_cov = QtWidgets.QTableWidgetItem()
                coverage = (org_count / len(self.db_files)) * 100
                coverage = float("%.2f" % coverage)
                perc_cov.setData(QtCore.Qt.EditRole, str(coverage) + '%')
                perc_cov.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 1, perc_cov)

                # push total count
                table_count = QtWidgets.QTableWidgetItem()
                table_count.setData(QtCore.Qt.EditRole, total_count)
                table_count.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 2, table_count)

                # push avg repeat
                avg_rep = QtWidgets.QTableWidgetItem()
                avg_rep_per_scaff = total_count / org_count
                avg_rep_per_scaff = float("%.2f" % avg_rep_per_scaff)
                avg_rep.setData(QtCore.Qt.EditRole, avg_rep_per_scaff)
                avg_rep.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 3, avg_rep)

                # push seq
                seq = QtWidgets.QTableWidgetItem()
                seq.setData(QtCore.Qt.EditRole, fives[majority_index] + seed + threes[majority_index])
                seq.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 4, seq)

                # push percent consensus
                perc_cons = QtWidgets.QTableWidgetItem()
                percent_consensus = (pams.count(pams[majority_index]) / len(pams)) * 100
                percent_consensus = float("%.2f" % percent_consensus)
                perc_cons.setData(QtCore.Qt.EditRole, str(percent_consensus) + "%")
                perc_cons.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 5, perc_cons)

                # push score
                score = QtWidgets.QTableWidgetItem()
                score.setData(QtCore.Qt.EditRole, scores[majority_index])
                score.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 6, score)

                # push PAM
                pam = QtWidgets.QTableWidgetItem()
                pam.setData(QtCore.Qt.EditRole, pams[majority_index])
                pam.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 7, pam)

                # push strand
                strand_val = ""
                if int(locs[majority_index]) < 0:
                    strand_val = "-"
                else:
                    strand_val = "+"

                strand = QtWidgets.QTableWidgetItem()
                strand.setData(QtCore.Qt.EditRole, strand_val)
                strand.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 8, strand)

                index += 1
                running_val += increase_val
                self.loading_window.loading_bar.setValue(running_val)

        self.table2.resizeColumnsToContents()
        self.loading_window.loading_bar.setValue(25)
        if len(self.db_files) > 1:
            self.plot_3D_graph()
        else:
            self.colormap_canvas.figure.set_visible(False)

        self.loading_window.loading_bar.setValue(100)
        self.loading_window.hide()
        QtCore.QCoreApplication.processEvents()


    def custom_seed_search(self):
        seeds = str(self.seed_input.text())
        seeds = seeds.replace(" ","")
        seeds = seeds.upper()
        seeds = seeds.split(",")

        if len(seeds) == 1:
            if seeds[0] == "":
                self.pre_analyze()
                return

        # update progress bar
        self.loading_window.loading_bar.setValue(5)
        self.loading_window.show()
        QtCore.QCoreApplication.processEvents()

        # prep table
        self.total_org_number = len(self.cspr_files)
        self.table2.setRowCount(0)
        self.loading_window.loading_bar.setValue(10)
        index = 0

        if (len(seeds) == 0):
            self.loading_window.hide()
            return

        increase_val = float(15 / len(self.seeds))
        running_val = self.loading_window.loading_bar.value()
        self.loading_window.info_label.setText("Parsing Seed Data")
        # QtCore.QCoreApplication.processEvents()

        # retrieve data on shared seeds
        self.counts = []
        for seed in seeds:
            total_count = 0
            org_count = 0
            threes = []
            fives = []
            scores = []
            pams = []
            locs = []
            none_data = True
            for db_file in self.db_files:
                conn = sqlite3.connect(db_file)
                c = conn.cursor()
                data = c.execute("SELECT count, three, five, pam, score, location FROM repeats WHERE seed = ?",
                                 (seed,)).fetchone()
                if data != None:
                    none_data = False
                    data = list(data)
                    org_count += 1
                    total_count += int(data[0])
                    threes += data[1].split(",")
                    fives += data[2].split(",")
                    pams += data[3].split(",")
                    scores += data[4].split(",")
                    locs += data[5].split(",")

            if none_data == True:
                QtWidgets.QMessageBox.information(self, "Seed Error",
                                                  seed + " : No such seed exists in the repeats section of any organism selected.",
                                                  QtWidgets.QMessageBox.Ok)
                self.loading_window.hide()
            else:

                # increase row count
                self.table2.setRowCount(index + 1)

                # push seed to table
                table_seed = QtWidgets.QTableWidgetItem()
                table_seed.setData(QtCore.Qt.EditRole, seed)
                table_seed.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 0, table_seed)

                self.counts.append(total_count)

                if len(threes) < len(fives):
                    for i in range(len(fives) - len(threes)):
                        threes.append('')

                elif len(fives) < len(threes):
                    for i in range(len(threes) - len(fives)):
                        fives.append('')

                majority_index = 0
                if threes[0] == '':
                    majority = max(set(fives), key=fives.count)
                    majority_index = fives.index(majority)
                else:
                    majority = max(set(threes), key=threes.count)
                    majority_index = threes.index(majority)

                # push percent coverage
                perc_cov = QtWidgets.QTableWidgetItem()
                coverage = (org_count / len(self.db_files)) * 100
                coverage = float("%.2f" % coverage)
                perc_cov.setData(QtCore.Qt.EditRole, str(coverage) + '%')
                perc_cov.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 1, perc_cov)

                # push total count
                table_count = QtWidgets.QTableWidgetItem()
                table_count.setData(QtCore.Qt.EditRole, total_count)
                table_count.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 2, table_count)

                # push avg repeat
                avg_rep = QtWidgets.QTableWidgetItem()
                avg_rep_per_scaff = total_count / org_count
                avg_rep_per_scaff = float("%.2f" % avg_rep_per_scaff)
                avg_rep.setData(QtCore.Qt.EditRole, avg_rep_per_scaff)
                avg_rep.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 3, avg_rep)

                # push seq
                seq = QtWidgets.QTableWidgetItem()
                seq.setData(QtCore.Qt.EditRole, fives[majority_index] + seed + threes[majority_index])
                seq.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 4, seq)

                # push percent consensus
                perc_cons = QtWidgets.QTableWidgetItem()

                percent_consensus = (pams.count(pams[majority_index]) / len(pams)) * 100
                percent_consensus = float("%.2f" % percent_consensus)
                perc_cons.setData(QtCore.Qt.EditRole, str(percent_consensus) + "%")
                perc_cons.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 5, perc_cons)

                # push score
                score = QtWidgets.QTableWidgetItem()
                score.setData(QtCore.Qt.EditRole, scores[majority_index])
                score.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 6, score)

                # push PAM
                pam = QtWidgets.QTableWidgetItem()
                pam.setData(QtCore.Qt.EditRole, pams[majority_index])
                pam.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 7, pam)

                # push strand
                strand_val = ""
                if int(locs[majority_index]) < 0:
                    strand_val = "-"
                else:
                    strand_val = "+"

                strand = QtWidgets.QTableWidgetItem()
                strand.setData(QtCore.Qt.EditRole, strand_val)
                strand.setTextAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
                self.table2.setItem(index, 8, strand)

                index += 1
                running_val += increase_val
                self.loading_window.loading_bar.setValue(running_val)

        self.table2.resizeColumnsToContents()
        self.loading_window.hide()
        QtCore.QCoreApplication.processEvents()


    #db_files is an array of database files for the organisms that will be looked at for shared seeds
    def get_shared_seeds(self, db_files, limit=False):
        #vars
        aliases = []

        #get db attachment aliases
        for i in range(1, len(db_files) + 1):
            aliases.append("main" + str(i))

        #memory connections for inner join on db files to hold what seeds are shared
        new_conn = sqlite3.connect(GlobalSettings.appdir + "temp_join.db")
        new_c = new_conn.cursor()
        new_c.execute("PRAGMA synchronous = OFF;")
        new_c.execute("PRAGMA journal_mode = OFF;")
        new_c.execute("PRAGMA locking_mode = EXCLUSIVE;")
        new_c.execute("DROP TABLE IF EXISTS repeats;")
        new_c.execute("VACUUM;")
        new_c.execute("DROP TABLE IF EXISTS join_results;")
        new_c.execute("CREATE table join_results (seed TEXT PRIMARY KEY);")

        #attach each db file with an alias
        for i in range(len(db_files)):
            new_c.execute("ATTACH DATABASE '" + db_files[i] + "' AS " + aliases[i] + ";")

        # start transaction
        new_c.execute("BEGIN TRANSACTION;")

        sql_inner_join = "INSERT into main.join_results select main1.repeats.seed from main1.repeats "

        for i in range(len(aliases[:-1])):
            sql_inner_join += "inner join " + aliases[i + 1] + ".repeats on "
            sql_inner_join += aliases[i] + ".repeats.seed = " + aliases[i + 1] + ".repeats.seed "

        #execute inner join
        new_c.execute(sql_inner_join)

        #get shared data
        if limit == False:
            shared_seeds = new_c.execute("select count(*) from join_results").fetchall()
            return shared_seeds
        else:
            shared_seeds = new_c.execute("select * from join_results limit 0,100").fetchall()

        #end transaction
        new_c.execute("END TRANSACTION;")

        #close memory db
        new_c.close()
        new_conn.close()

        #parse shared seeds into self.seeds
        seeds = []
        for tup in shared_seeds:
            seeds.append(tup[0])

        return seeds


    def get_org_names(self):
        self.org_names = {}
        for file in self.cspr_files:
            with gzip.open(file, "r") as f:
                line = f.readline()
                buf = str(line)
                buf = buf.strip("'b")
                buf = buf[:len(buf) - 2]
                org_name = buf.replace("GENOME: ", "")

                line = f.readline()
                buf = str(line)
                buf = buf.strip("'b")
                buf = buf[:len(buf) - 2]
                kstats = buf.replace("KARYSTATS: ", "")
                kstats = kstats.split(",")
                self.org_names[org_name] = len(kstats) - 1

    ###This is deprecated
    """
    def plot_repeats_vs_seeds(self):
        self.pop_analysis_repeats_graph.canvas.figure.set_visible(True)
        self.pop_analysis_repeats_graph.canvas.axes.clear()
        x1 = list(range(0, len(self.seeds)))
        y1 = []
        for val in self.counts:
            y1.append(val)
        #y1.sort(reverse=True)

        # get stats
        self.average = statistics.mean(y1)
        self.mode = statistics.mode(y1)
        self.median = statistics.median(y1)
        self.repeat_count = len(self.seeds)

        # clear axes
        self.pop_analysis_repeats_graph.canvas.axes.clear()
        # the following are for plotting / formatting
        self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1)
        self.pop_analysis_repeats_graph.canvas.axes.set_xlabel('Seed ID Number')
        self.pop_analysis_repeats_graph.canvas.axes.set_ylabel('Number of Repeats')
        self.pop_analysis_repeats_graph.canvas.axes.set_title('Number of Repeats per Seed ID Number')
        # always redraw at the end
        self.pop_analysis_repeats_graph.canvas.draw()
        QtCore.QCoreApplication.processEvents()
        """
        
    def plot_3D_graph(self):
        for i in reversed(range(self.colormap_layout.count())): ### Clear out old widges in layout
            self.colormap_layout.itemAt(i).widget().setParent(None)
            
        self.colormap_canvas = MplCanvas(self, width=5, height=4, dpi=100) ###Initialize new Canvas
        self.colormap_layout.addWidget(self.colormap_canvas) ### Add canvas to colormap layout 
        self.colormap_figure.setLayout(self.colormap_layout) ### Add colormap layout to colormap plot widget

#        self.colormap_canvas.axes.clear()
#        try:
#            self.colormap_canvas.cbar.remove()
#        except:
#            None
#        self.colormap_canvas.figure.set_visible(True)

        rows, cols = (self.total_org_number, self.total_org_number)
        arr = [[0 for i in range(cols)] for j in range(rows)]

        self.names = list(self.org_names.keys())

        for pair in list(itertools.combinations(self.db_files, 2)):
            shared_seeds = list(self.get_shared_seeds(list(pair), False)[0])[0]

            arr[self.db_files.index(pair[0])][self.db_files.index(pair[1])] += int(shared_seeds)
            arr[self.db_files.index(pair[1])][self.db_files.index(pair[0])] += int(shared_seeds)

        for i in range(len(arr)):
            conn = sqlite3.connect(self.db_files[i])
            c = conn.cursor()
            arr[i][i] = int(list(c.execute("select count(*) from repeats;").fetchall()[0])[0])
            c.close()
            conn.close()

        labels = copy.deepcopy(arr)

        for i in range(len(arr)):
            arr[i][i] = 0

        ax = self.colormap_canvas.axes
        im = self.colormap_canvas.axes.imshow(arr, cmap='summer')
        self.colormap_canvas.cbar = self.colormap_canvas.axes.figure.colorbar(im, ax=self.colormap_canvas.axes)
        self.colormap_canvas.cbar.ax.set_ylabel("", rotation=-90, va="bottom",fontsize=8)
        cursor = mplcursors.cursor(im, hover=True)
        @cursor.connect("add")
        def on_add(sel):
            sel.annotation.arrow_patch.set(arrowstyle="simple", fc="white", alpha=.5)
            sel.annotation.set_bbox(None)
            i,j = sel.target.index
            sel.annotation.set_text(labels[i][j])

        ax.set_xticks(np.arange(len(arr)))
        ax.set_yticks(np.arange(len(arr)))

        #get labels based on org table rows
        ax.set_xticklabels(self.rows)
        ax.set_yticklabels(self.rows)
        ax.set_xlabel("Organism", fontsize = 10)
        ax.set_ylabel("Organism", fontsize = 10)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=2)
        #ax.set_frame_on(False)
        #ax.set_aspect('equal')
#        for i in range(len(arr)):
#            for j in range(len(arr)):
#                text = ax.text(j, i, arr[i][j],
#                               ha="center", va="center", color="black")

        self.colormap_canvas.draw()

    def find_locations(self):

        #error checking
        if len(self.table2.selectedItems()) == 0:
            QtWidgets.QMessageBox.question(self, "Error", "Please select at least 1 seed to find locations of.",
                                                QtWidgets.QMessageBox.Ok)
            self.loc_finder_table.setRowCount(0)
            return

        #get seeds from selected rows in table
        seeds = []
        for i in self.table2.selectionModel().selectedRows():
            item = self.table2.item(i.row(), 0)
            seeds.append(item.text())

        index = 0
        org_names = list(self.org_names.keys())
        #loop through each db table looking for the seed, if found, insert data in to locations table
        for db_file in self.db_files:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            for seed in seeds:
                data = c.execute("SELECT chromosome, location, five, three FROM repeats WHERE seed = ? ", (seed,)).fetchone()
                #make sure seed was found in table, then parse and store in table
                if data != None:
                    data = list(data)
                    chroms = data[0].split(',')
                    locs = data[1].split(',')
                    fives = data[2].split(',')
                    threes = data[3].split(',')
                    if threes[0] == '':
                        threes = []
                    if fives[0] == '':
                        fives = []
                    i = 0
                    for chrom in chroms:
                        self.loc_finder_table.setRowCount(index + 1)
                        seed_table = QtWidgets.QTableWidgetItem()
                        sequence_table = QtWidgets.QTableWidgetItem()
                        organism_table = QtWidgets.QTableWidgetItem()
                        chromsome_table = QtWidgets.QTableWidgetItem()
                        location_table = QtWidgets.QTableWidgetItem()
                        seed_table.setData(QtCore.Qt.EditRole, seed)
                        if len(fives) == 0 and len(threes) != 0:
                            sequence_table.setData(QtCore.Qt.EditRole, seed + threes[i])
                        elif len(threes) == 0 and len(fives) != 0:
                            sequence_table.setData(QtCore.Qt.EditRole, fives[i] + seed)
                        else:
                            sequence_table.setData(QtCore.Qt.EditRole, fives[i] + seed + threes[i])
                        organism_table.setData(QtCore.Qt.EditRole, org_names[self.db_files.index(db_file)])
                        chromsome_table.setData(QtCore.Qt.EditRole, chrom)
                        location_table.setData(QtCore.Qt.EditRole, abs(int(locs[i])))
                        self.loc_finder_table.setItem(index, 0, seed_table)
                        self.loc_finder_table.setItem(index, 1, sequence_table)
                        self.loc_finder_table.setItem(index, 2, organism_table)
                        self.loc_finder_table.setItem(index, 3, chromsome_table)
                        self.loc_finder_table.setItem(index, 4, location_table)
                        i += 1
                        index += 1

        self.loc_finder_table.resizeColumnsToContents()


    # this function clears the loc_finder_table
    def clear_loc_table(self):
        self.loc_finder_table.clearContents()
        self.loc_finder_table.setRowCount(0)


    def table_sorting(self, logicalIndex):
        self.switcher[logicalIndex] *= -1
        if self.switcher[logicalIndex] == -1:
            self.table2.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.table2.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)


    # sorting to table2: IE the table in top-right
    def table2_sorting(self, logicalIndex):
        self.switcher_table2[logicalIndex] *= -1
        if self.switcher_table2[logicalIndex] == -1:
            self.table2.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.table2.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)


    # sorting for location table: IE table in bottom right
    def loc_table_sorter(self, logicalIndex):
        self.switcher_loc_table[logicalIndex] *= -1
        if (self.switcher_loc_table[logicalIndex] == -1):
            self.loc_finder_table.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.loc_finder_table.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)


    def clear(self):
        self.table2.setRowCount(0)


    def go_back(self):
        GlobalSettings.mainWindow.getData()
        GlobalSettings.mainWindow.mwfg.moveCenter(GlobalSettings.mainWindow.cp)  ##Center window
        GlobalSettings.mainWindow.move(GlobalSettings.mainWindow.mwfg.topLeft())  ##Center window
        GlobalSettings.mainWindow.show()
        self.hide()


    # this function calls the close window class. Allows the user to choose what files they want to keep/delete
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()


class loading_window(QtWidgets.QWidget):
    def __init__(self):
        super(loading_window, self).__init__()
        uic.loadUi(GlobalSettings.appdir + "loading_data_form.ui", self)
        self.loading_bar.setValue(0)
        self.setWindowTitle("Loading Data")
        self.hide()

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi,tight_layout=True)
        self.axes = fig.add_subplot(111)
        self.axes.clear()
        super(MplCanvas, self).__init__(fig)