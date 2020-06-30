from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import os
from CSPRparser import CSPRparser
import Algorithms
from functools import partial
import numpy as np
from PyQt5.QtWidgets import *
from matplotlib_venn import venn3_unweighted
import show_names_ui
import show_names2_ui
import sys
import gzip
import sqlite3
from collections import Counter
import statistics

class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'populationanalysis.ui', self)
        self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
        self.goBackButton.clicked.connect(self.go_back)
        self.analyze_button.clicked.connect(self.pre_analyze)
        self.clear_Button.clicked.connect(self.clear)
        self.ncbi_search_button.clicked.connect(self.launch_ncbi_seacher)
        self.meta_genomic_cspr_checkbox.stateChanged.connect(self.get_data)
        self.parser = CSPRparser("")
        self.Endos = dict()
        self.fna_files = dict()
        self.cspr_files = {}
        self.sq = Algorithms.SeqTranslate()
        self.ref_para_list = list()
        self.mode = 0
        self.find_locs_button.clicked.connect(self.find_locations)
        self.clear_loc_button.clicked.connect(self.clear_loc_table)
        self.directory = ""
        self.names = []
        self.names_venn = []
        self.show_names.clicked.connect(self.show_names_func)
        self.show_names_2.clicked.connect(self.show_names_func2)
        self.name_form = show_names_ui.show_names_table()
        self.name_form2 = show_names2_ui.show_names_table2()


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
        self.loc_finder_table.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.loc_finder_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.loc_finder_table.resizeColumnsToContents()

        # action buttons
        self.actionMetaGenome_Parser.triggered.connect(self.launch_chrom_selector)

        self.combinerWindow = fna_and_cspr_combiner()

        self.total_org_number = 0

        self.switcher_table2 = [1, 1, 1, 1, 1, 1, 1, 1,
                                1]  # for keeping track of where we are in the sorting clicking for each column
        self.switcher_loc_table = [1, 1, 1, 1, 1]

        self.mwfg = self.frameGeometry()
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()

        self.cspr_file = ""
        self.db_file = ""


    def launch_ncbi_seacher(self):
        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(0)
        GlobalSettings.mainWindow.ncbi_search_dialog.show()


    # launches the chromesome selector window
    def launch_chrom_selector(self):
        GlobalSettings.mainWindow.cspr_selector.launch(self.cspr_files)


    # this function calls the popParser function and fills all the tables
    def pre_analyze(self):
        # if the user is wanting to go with 1 meta genomic cspr file
        if self.meta_genomic_cspr_checkbox.isChecked():
            selectedList = self.org_Table.selectedItems()
            # error check
            if len(selectedList) == 0 or len(selectedList) > 1:
                QtWidgets.QMessageBox.question(self, "Error",
                                               "Please select no more than 1 CSPR file for analysis.",
                                               QtWidgets.QMessageBox.Ok)
                return

            orgName = selectedList[0].text()
            self.cspr_file = str(self.fna_files[orgName])
            self.db_file = self.cspr_file.strip('.cspr') + '_repeats.db'
            self.fill_data()



        # if the user is wanting to go with creating a new meta genomic cspr file
        else:
            selectedList = self.org_Table.selectedItems()

            # if the table is showing only fna/fasta files
            if not self.meta_genomic_cspr_checkbox.isChecked():
                # rules for selecting FNA/Fasta files
                # check to make sure that the user selected at least 2 organisms, and 1 endonuclease
                if len(selectedList) < 1 or self.endoBox.currentText() == 'None Selected':
                    QtWidgets.QMessageBox.question(self, "Nothing Seleted",
                                                   "No items selected. Please select at least 1 organism and only 1 endonuclease.",
                                                   QtWidgets.QMessageBox.Ok)
                    return
                if len(selectedList) == 1:
                    error = QtWidgets.QMessageBox.question(self, "Only 1 Organism Selected",
                                                           "Population Analysis works with multiple organisms, or a metagenome. If the file selected is not a metagenome, the program may not function correctly. Do you wish to continue? \n\n"
                                                           "Do you wish to continue?",
                                                           QtWidgets.QMessageBox.Yes |
                                                           QtWidgets.QMessageBox.No,
                                                           QtWidgets.QMessageBox.No)
                    if (error == QtWidgets.QMessageBox.No):
                        return

                submitList = list()
                for item in selectedList:
                    submitList.append(self.fna_files[item.text()])
                self.combinerWindow.launch(submitList)
            # rules for selecting cspr files
            elif self.meta_genomic_cspr_checkbox.isChecked():
                if len(selectedList) == 0:
                    QtWidgets.QMessageBox.question(self, "Nothing Seleted",
                                                   "No items selected. Please select one metagenome for Population Analysis.",
                                                   QtWidgets.QMessageBox.Ok)
                    return
                elif len(selectedList) > 1:
                    QtWidgets.QMessageBox.question(self, "Too Many Selected",
                                                   "Only 1 metagenomic CSPR file is allowed to be selected.",
                                                   QtWidgets.QMessageBox.Ok)
                    return


    def launch(self, path):
        os.chdir(path)
        self.directory = path
        self.get_data()


    def get_data(self):
        if self.directory == '':
            return

        onlyfiles = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        self.fna_files.clear()

        # show/hide the stuff that isn't needed
        if self.meta_genomic_cspr_checkbox.isChecked():
            self.endoBox.hide()
            self.ncbi_search_button.hide()
            self.label_3.hide()
            self.label_2.setText('Select a Metagenomic CSPR File')
        elif not self.meta_genomic_cspr_checkbox.isChecked():
            self.endoBox.show()
            self.ncbi_search_button.show()
            self.label_3.show()
            self.label_2.setText('Select organism(s) and endonuclease:')

        # if the user wants the FNA/Fast files to be sown
        if not self.meta_genomic_cspr_checkbox.isChecked():
            index = 0
            for file in onlyfiles:
                if file.find('.fna') != -1 or file.find('.fasta') != -1:
                    # find the organism name
                    f = open(file, 'r')
                    hold = f.readline()
                    f.close()
                    hold = str(hold)
                    hold = hold[:len(hold) - 4]
                    # get the organism name
                    spaceIndex = hold.find(' ') + 1
                    commaIndex = hold.find(',')
                    buf = hold[spaceIndex:commaIndex]

                    # store the name in the dict of fna_files, that keys the name with the file path
                    self.fna_files[buf] = file

                    # store the data in the table
                    tabWidget = QtWidgets.QTableWidgetItem(buf)
                    self.org_Table.setRowCount(index + 1)
                    self.org_Table.setItem(index, 0, tabWidget)
                    index += 1
        # if the user wants the metagenomic cspr files to be shown
        else:
            index = 0
            for file in onlyfiles:
                if file.find('.cspr') != -1:
                    f = gzip.open(file, 'r')
                    hold = f.readline()
                    f.close()
                    hold = str(hold)
                    hold = hold.strip("'b")
                    hold = hold[:len(hold) - 4]
                    # only show the files that are metagenomic
                    if '(meta)' in hold:
                        colonIndex = hold.find(':') + 1
                        commaIndex = hold.find(',')
                        orgName = hold[colonIndex:commaIndex]

                        self.fna_files[orgName] = file

                        tabWidget = QtWidgets.QTableWidgetItem(orgName)
                        self.org_Table.setRowCount(index + 1)
                        self.org_Table.setItem(index, 0, tabWidget)
                        index += 1
            if index == 0:
                self.org_Table.clearContents()
                self.org_Table.setRowCount(0)

        #        self.org_Table.resizeColumnsToContents() ##Commenting this out allows the header to remain full sized

        self.fillEndo()
        # self.changeEndos()


    # this function opens CASPERinfo and builds the dropdown menu of selectable endonucleases
    def fillEndo(self):
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
                    default_seed_length = line_tokened[2]
                    default_tot_length = line_tokened[3]
                    self.Endos[endo + " PAM: " + p_pam] = (endo, p_pam, default_seed_length, default_tot_length)

                break
        f.close()
        self.endoBox.addItem("None Selected")
        self.endoBox.addItems(self.Endos.keys())
        # self.endoBox.currentIndexChanged.connect(self.changeEndos)


    def fill_data(self):
        #get misc line from cspr file to know the number of orgs in file
        with gzip.open(self.cspr_file,'r') as f:
            i = 0
            self.org_list = {}
            for line in f:
                line = str(line).strip(r"'b\r\n")
                if i == 2:
                    colonIndex = line.find(':') + 2
                    usefulData = line[colonIndex:]
                    usefulData = usefulData.split('|')
                    usefulData.pop()
                    j = 0
                    while j < len(usefulData):
                        temp = usefulData[j].split(',')
                        if temp[0] not in self.org_list.keys():
                            self.org_list[temp[0]] = 1
                        else:
                            self.org_list[temp[0]] += 1
                        j += 1
                    break
                i += 1
        self.total_org_number = len(self.org_list.keys())
        self.table2.setRowCount(0)
        index = 0
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        for seeds in c.execute("select * from repeats"):
            self.table2.setRowCount(index + 1)
            seeds = list(seeds)

            #seed
            seed = QtWidgets.QTableWidgetItem()
            s = seeds[0]
            seed.setData(QtCore.Qt.EditRole, str(s))
            seed.setTextAlignment(QtCore.Qt.AlignHCenter)

            #total repeats
            total_repeats = QtWidgets.QTableWidgetItem()
            count = seeds[-1]
            total_repeats.setData(QtCore.Qt.EditRole, count)
            total_repeats.setTextAlignment(QtCore.Qt.AlignHCenter)

            #score
            score = QtWidgets.QTableWidgetItem()
            sc = str(seeds[4])
            sc = sc.split(',')
            c = Counter(sc)
            score.setData(QtCore.Qt.EditRole, list(c.most_common(1)[0])[0])
            score.setTextAlignment(QtCore.Qt.AlignHCenter)

            #PAM and Strand
            tail = str(seeds[3])
            tail = tail.split(',')[0]
            if tail.find('+') != -1:
                st = '+'
                tail = tail.split('+')
                p = tail[1]
            else:
                st = '-'
                tail = tail.split('-')
                p = tail[1]
            strand = QtWidgets.QTableWidgetItem()
            pam = QtWidgets.QTableWidgetItem()
            strand.setData(QtCore.Qt.EditRole, st)
            pam.setData(QtCore.Qt.EditRole, p)
            pam.setTextAlignment(QtCore.Qt.AlignHCenter)
            strand.setTextAlignment(QtCore.Qt.AlignHCenter)


            #location counts
            counts = {}
            locs = str(seeds[1])
            locs = locs.split(',')
            if locs.count(locs[0]) == len(locs):
                ind = int(locs[0])
                cnt = 0
                org = ""
                for k in self.org_list.keys():
                    cnt += self.org_list[k]
                    if ind <= cnt:
                        org = k
                        break
                counts[org] = len(locs)
            else:
                for loc in locs:
                    ind = int(loc)
                    cnt = 0
                    org = ""
                    for k in self.org_list.keys():
                        cnt += self.org_list[k]
                        if ind <= cnt:
                            org = k
                            break
                    if org not in counts.keys():
                        counts[org] = 1
                    else:
                        counts[org] += 1

            #avg repeats per scaffold
            scaffold_count = 0
            rep_count = 0
            for key in counts.keys():
                scaffold_count += 1
                rep_count += counts[key]
            avg_rep_per_scaff = rep_count / scaffold_count
            avg_rep_per_scaff = float("%.2f" % avg_rep_per_scaff)
            avg_rep = QtWidgets.QTableWidgetItem()
            avg_rep.setData(QtCore.Qt.EditRole, str(avg_rep_per_scaff))
            avg_rep.setTextAlignment(QtCore.Qt.AlignHCenter)


            #% coverage
            coverage = (len(counts) / len(self.org_list.keys())) * 100
            coverage = float("%.2f" % coverage)
            perc_coverage = QtWidgets.QTableWidgetItem(str(coverage) + '%')
            perc_coverage.setTextAlignment(QtCore.Qt.AlignHCenter)


            #%consensus
            tails = str(seeds[3]).split(',')
            pam_counts = {}
            for tail in tails:
                if tail.find('+') != -1:
                    tail = tail.split('+')
                    p = tail[1]
                else:
                    tail = tail.split('-')
                    p = tail[1]
                if p not in pam_counts.keys():
                    pam_counts[p] = 1
                else:
                    pam_counts[p] += 1

            perc_cons = (max(pam_counts.values()) / sum(pam_counts.values())) * 100
            perc_cons = float("%.2f" % perc_cons)
            consensus = QtWidgets.QTableWidgetItem(str(perc_cons) + '%')
            consensus.setTextAlignment(QtCore.Qt.AlignHCenter)

            #consensus seq
            temp = ""
            cnt = 0
            for p in pam_counts.keys():
                if pam_counts[p] > cnt:
                    cnt = pam_counts[p]
                    temp = p
            for tail in tails:
                if tail.find('+') != -1:
                    tail = tail.split('+')
                else:
                    tail = tail.split('-')
                if temp == tail[1]:
                    temp = tail[0]
                    break


            cons_seq = temp + s
            consensus_seq = QtWidgets.QTableWidgetItem()
            consensus_seq.setData(QtCore.Qt.EditRole, cons_seq)
            consensus_seq.setTextAlignment(QtCore.Qt.AlignHCenter)

            self.table2.setItem(index, 0, seed)
            self.table2.setItem(index, 1, perc_coverage)
            self.table2.setItem(index, 2, total_repeats)
            self.table2.setItem(index, 3, avg_rep)
            self.table2.setItem(index, 4, consensus_seq)
            self.table2.setItem(index, 5, consensus)
            self.table2.setItem(index, 6, score)
            self.table2.setItem(index, 7, pam)
            self.table2.setItem(index, 8, strand)

            index += 1
            #break
        self.table2.resizeColumnsToContents()
        self.plot_repeats_vs_seeds()
        self.plot_3D_graph()
        self.plot_venn()


    def plot_repeats_vs_seeds(self):
        self.pop_analysis_repeats_graph.canvas.axes.clear()
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT seed, count from repeats").fetchall()
        c.close()
        x1 = list(range(0, len(data)))
        y1 = []
        for obj in data:
            y1.append(int(list(obj)[1]))
        y1.sort(reverse=True)

        # get stats
        self.average = statistics.mean(y1)
        self.mode = statistics.mode(y1)
        self.median = statistics.median(y1)
        self.repeat_count = len(data)

        # clear axes
        self.pop_analysis_repeats_graph.canvas.axes.clear()
        # the following are for plotting / formatting
        self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1)
        self.pop_analysis_repeats_graph.canvas.axes.set_xlabel('Seed ID Number')
        self.pop_analysis_repeats_graph.canvas.axes.set_ylabel('Number of Repeats')
        self.pop_analysis_repeats_graph.canvas.axes.set_title('Number of Repeats per Seed ID Number')
        # always redraw at the end
        self.pop_analysis_repeats_graph.canvas.draw()


    def plot_3D_graph(self):
        rows, cols = (self.total_org_number, self.total_org_number)
        arr = [[0 for i in range(cols)] for j in range(rows)]

        x3 = []
        y3 = []
        z3 = np.zeros(int((self.total_org_number * (self.total_org_number - 1)) / 2))
        dz = []
        self.names = list(self.org_list.keys())
        dx = np.ones(int((self.total_org_number * (self.total_org_number - 1)) / 2))
        dy = np.ones(int((self.total_org_number * (self.total_org_number - 1)) / 2))
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        for seeds in c.execute("select * from repeats"):
            seeds = list(seeds)
            temp_names = []
            chroms = str(seeds[1]).split(',')
            for c in chroms:
                ind = int(c)
                cnt = 0
                c_name = ""
                for k in self.org_list.keys():
                    cnt += self.org_list[k]
                    if ind <= cnt:
                        c_name = k
                        break
                if c_name not in temp_names:
                    temp_names.append(c_name)

            if len(temp_names) >= 2:
                for i in range(len(temp_names)-1):
                    j = i + 1
                    while j != len(temp_names):
                        arr[self.names.index(temp_names[i])][self.names.index(temp_names[j])] += 1
                        arr[self.names.index(temp_names[j])][self.names.index(temp_names[i])] += 1
                        j += 1

        for j in range(cols):
            i = len(self.names)-1
            while i != j:
                x3.append(i)
                y3.append(j)
                dz.append(arr[i][j])
                i -= 1

        self.pop_analysis_3dgraph.canvas.axes.clear()
        self.pop_analysis_3dgraph.canvas.axes.bar3d(x3, y3, z3, dx, dy, dz)

        new_names = []
        for n in range(len(self.names)):
            new_names.append(n)


        self.pop_analysis_3dgraph.canvas.axes.set_xlabel('x')
        self.pop_analysis_3dgraph.canvas.axes.set_ylabel('y')
        self.pop_analysis_3dgraph.canvas.axes.set_zlabel('z')
        self.pop_analysis_3dgraph.canvas.axes.set_xticks(np.arange(1, self.total_org_number + 1, 1))
        self.pop_analysis_3dgraph.canvas.axes.set_yticks(np.arange(0, self.total_org_number, 1))
        self.pop_analysis_3dgraph.canvas.axes.tick_params(labelsize=8)
        self.pop_analysis_3dgraph.canvas.axes.set_xticklabels(new_names, rotation=45)
        self.pop_analysis_3dgraph.canvas.axes.set_yticklabels(new_names, rotation=-45)
        self.pop_analysis_3dgraph.canvas.draw()


    def plot_venn(self):
        self.pop_analysis_venn_diagram.canvas.figure.clf()
        rows, cols = (self.total_org_number, self.total_org_number)
        arr = [[0 for i in range(cols)] for j in range(rows)]
        all_3 = 0
        singles = [0 for i in range(cols)]
        self.names_venn = list(self.org_list.keys())

        if len(self.names_venn) >= 3:
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            for seeds in c.execute("select * from repeats"):
                seeds = list(seeds)
                temp_names = []
                chroms = str(seeds[1]).split(',')
                for c in chroms:
                    ind = int(c)
                    cnt = 0
                    c_name = ""
                    for k in self.org_list.keys():
                        cnt += self.org_list[k]
                        if ind <= cnt:
                            c_name = k
                            break
                    if c_name not in temp_names:
                        temp_names.append(c_name)

                if len(temp_names) >= 2:
                    for i in range(len(temp_names) - 1):
                        j = i + 1
                        while j != len(temp_names):
                            arr[self.names_venn.index(temp_names[i])][self.names_venn.index(temp_names[j])] += 1
                            arr[self.names_venn.index(temp_names[j])][self.names_venn.index(temp_names[i])] += 1
                            j += 1
                else:
                    singles[self.names_venn.index(temp_names[0])] += 1

                if all(x in temp_names for x in [self.names_venn[0], self.names_venn[1], self.names_venn[2]]):
                    all_3 += 1

            venn3_unweighted(subsets=(singles[0], singles[1], arr[0][1], singles[2], arr[0][2],
                                      arr[1][2], all_3), set_labels=('0', '1', '2'))
            self.pop_analysis_venn_diagram.canvas.draw()

        else:
            self.pop_analysis_venn_diagram.canvas.figure.clf()
            self.pop_analysis_venn_diagram.canvas.draw()


    def find_locations(self):

        #error checking
        if len(self.table2.selectedItems()) == 0:
            QtWidgets.QMessageBox.question(self, "Error", "Please select at least 1 seed to find locations of.",
                                                QtWidgets.QMessageBox.Ok)
            self.loc_finder_table.setRowCount(0)
            return

        #get data from sql table on seeds selected and fill in locations table
        seeds = []
        i = 0
        for item in self.table2.selectedItems():
            if i % 9 == 0:
                seeds.append(item.text())
            i += 1
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        index = 0
        for seed in seeds:
            data = c.execute("SELECT * FROM repeats WHERE seed = ? ", (seed,)).fetchone()
            data = list(data)
            chroms = str(data[1]).split(',')
            locs = str(data[2]).split(',')
            tails = str(data[3]).split(',')
            tail_data = []
            for tail in tails:
                if tail.find('+') != -1:
                    tail = tail.split('+')
                    tail_data.append(tail[0])
                else:
                    tail = tail.split('-')
                    tail_data.append(tail[0])
            i = 0
            for chrom in chroms:
                ind = int(chrom)
                cnt = 0
                org = ""
                for k in self.org_list.keys():
                    cnt += self.org_list[k]
                    if ind <= cnt:
                        org = k
                        break

                self.loc_finder_table.setRowCount(index + 1)
                seed_table = QtWidgets.QTableWidgetItem()
                sequence_table = QtWidgets.QTableWidgetItem()
                organism_table = QtWidgets.QTableWidgetItem()
                chromsome_table = QtWidgets.QTableWidgetItem()
                location_table = QtWidgets.QTableWidgetItem()

                seed_table.setData(QtCore.Qt.EditRole, seed)
                sequence_table.setData(QtCore.Qt.EditRole, tail_data[i] + seed)
                organism_table.setData(QtCore.Qt.EditRole, org)
                chromsome_table.setData(QtCore.Qt.EditRole, chrom)
                location_table.setData(QtCore.Qt.EditRole, locs[i])
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


    def show_names_func(self):
        # print(self.names)
        self.name_form.fill_table(self.names)
        self.name_form.show()


    def show_names_func2(self):
        # print(self.names)
        if len(self.names_venn) >= 3:
            self.name_form2.fill_table(self.names_venn[0:3])
            self.name_form2.show()
        else:
            self.name_form2.name_table2.setRowCount(0)


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
        GlobalSettings.mainWindow.show()
        self.hide()


    # this function calls the close window class. Allows the user to choose what files they want to keep/delete
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()









class fna_and_cspr_combiner(QtWidgets.QDialog):
    def __init__(self):
        # Qt init stuff
        super(fna_and_cspr_combiner, self).__init__()
        uic.loadUi(GlobalSettings.appdir + "pop_analysis_fna_combiner.ui", self)
        self.sequencer_prog_bar.setValue(0)

        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.start_button.clicked.connect(self.run_analysis)

        # variables below
        self.ref_para_list = list()
        self.fna_file_names = list()
        self.sq = Algorithms.SeqTranslate()
        self.proc_running = False
        self.combined_fna_file = ''
        self.cancelled = False

    # for when the user clicks the 'x' button
    def closeEvent(self, event):
        closeWindow = self.cancel_function()

        # if the user is doing OT and does not decide to cancel it ignore the event
        if closeWindow == -2:
            event.ignore()
        else:
            event.accept()

    # this is the function that will eventually run the analysis. For now, it is just running the combine/build new cspr function
    def run_analysis(self):
        # make sure the process isn't already running
        if self.proc_running:
            return

        # make sure the user inputs whats needed
        if self.orgName_line_edit.text() == '' or self.org_code_line_edit.text() == '' or self.orgNum_lineEdit.text() == '':
            QtWidgets.QMessageBox.question(self, "Missing Information",
                                           "Please input an organism name, an organism code, and the number of organisms you are analyzing.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # make sure the user inputs an integer for the number of organisms
        if not self.orgNum_lineEdit.text().isdigit():
            QtWidgets.QMessageBox.question(self, "Error",
                                           "Organism number must be integers only!",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.combine_fna_files()
        self.build_new_cspr_file()

    # cancel function, just clears everything and closes the window
    def cancel_function(self):
        # check to see if the sequencer is running. If so ask the user if they wish to close out
        if self.proc_running:
            error = QtWidgets.QMessageBox.question(self, "Sequencer Is Running",
                                                   "Sequencer is running. Closing this window will cancel that process and return to the Population Analysis window. \n\n"
                                                   "Do you wish to continue?",
                                                   QtWidgets.QMessageBox.Yes |
                                                   QtWidgets.QMessageBox.No,
                                                   QtWidgets.QMessageBox.No)
            if (error == QtWidgets.QMessageBox.No):
                return -2
            else:
                self.cancelled = True
                self.process.kill()

        # close out and leave
        self.orgName_line_edit.setText('')
        self.org_code_line_edit.setText('')
        self.sequencer_prog_bar.setValue(0)
        self.hide()

    # open the window and get the fna files that the user wishes to use
    def launch(self, filenames):
        self.fna_file_names = filenames
        self.process = QtCore.QProcess()
        self.show()

    # this function takes a list of file paths, which are FNA or Fasta files, and combines them into one file
    # this function also builds a parallel/reference list of the chromosomes
    def combine_fna_files(self):
        self.ref_para_list.clear()

        # open the output file (currently just a test file)
        self.combined_fna_file = GlobalSettings.CSPR_DB + '/' + 'temp.fna'
        out_stream = open(self.combined_fna_file, 'w')

        # for each file in the list
        for file in self.fna_file_names:
            file_stream = open(file, 'r')

            buf = file_stream.readline()

            spaceIndex = buf.find(' ') + 1
            commaIndex = buf.find(',')
            orgName = buf[spaceIndex:commaIndex]

            # read the whole file
            while buf != "":

                # if buf is empty, break out
                if buf == '' or buf == '\n':
                    break

                # if buf starts with a '>', it's a chromosome so make sure to store it in the ref_para list
                if buf.startswith('>'):
                    self.ref_para_list.append((orgName, buf))
                # right the buf to the new file, and re-read it
                out_stream.write(buf)
                buf = file_stream.readline()
            file_stream.close()
        out_stream.close()

    # this function builds a new cspr file from the combined FNA file
    # very similar to New Genome
    def build_new_cspr_file(self):
        self.num_chromo_next = False
        self.num_chromo = 0

        # this function reads output. Just used to populate the progress bar
        def output_stdout(p):
            line = str(p.readAll())
            line = line[2:]
            line = line[:len(line) - 1]
            line = line.strip('\\n')
            line = line.strip('\n')
            for lines in filter(None, line.split(r'\n')):
                lines = lines.strip('\n')
                lines = lines.strip('\\n')
                if (lines == 'Finished reading in the genome file.'):
                    self.num_chromo_next = True
                elif (self.num_chromo_next == True):
                    lines = lines.strip('\\n')
                    self.num_chromo_next = False
                    self.num_chromo = int(lines)
                elif (lines.find('Chromosome') != -1 and lines.find('complete.') != -1):
                    temp = lines
                    temp = temp.replace('Chromosome ', '')
                    temp = temp.replace(' complete.', '')
                    temp = temp.strip('\n')
                    if (int(temp) == self.num_chromo):
                        self.sequencer_prog_bar.setValue(99)
                    else:
                        self.sequencer_prog_bar.setValue(int(temp) / self.num_chromo * 100)
                elif (lines == 'Finished Creating File.'):
                    self.sequencer_prog_bar.setValue(100)

        # this function will end up doing stuff when the process is finished.
        def finished():
            self.proc_running = False

            if self.cancelled == False:
                # get the file name
                cspr_file_name = GlobalSettings.CSPR_DB + '/' + self.org_code_line_edit.text() + '_' + \
                                 GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0] + '.cspr'
                self.process.kill()
                endoChoice = GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0]
                GlobalSettings.pop_Analysis.total_org_number, GlobalSettings.pop_Analysis.ref_para_list = GlobalSettings.pop_Analysis.parser.popParser(
                    cspr_file_name, endoChoice)
                GlobalSettings.pop_Analysis.fill_data()
            os.remove(self.combined_fna_file)
            self.cancelled = False
            self.cancel_function()

        # --------------getting the arugments---------------------------------
        # get the file path to the combined fna file
        fna_file_path = self.combined_fna_file
        path_to_fna = fna_file_path
        # get the endo_choice, hard coded for now, will eventually from with the user
        endo_choice = GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0]
        # get the pam itself, taken from the endo_choice. Also, if there's multiple, always take the first one
        pam = self.sq.endo_info[endo_choice][0].split(',')[0]
        # get the code. currently hard coded to test_code until I know what to put it to later one
        code = self.org_code_line_edit.text()
        # check to see if the seq_finder should do 5' prime or no
        if int(self.sq.endo_info[endo_choice][3]) == 3:
            pamdir = False
        else:
            pamdir = True
        # get the output folder location
        output_location = GlobalSettings.CSPR_DB
        # get the path to CASPERinfo
        path_to_info = GlobalSettings.appdir + 'CASPERinfo'
        # make org name something that will make more sense, this is just for testing right now
        orgName = self.orgName_line_edit.text() + '  , (meta), ' + self.orgNum_lineEdit.text()
        # get the seed and RNA length, based on the endo choice
        gRNA_length = self.sq.endo_info[endo_choice][2]
        seed_length = self.sq.endo_info[endo_choice][1]

        secondCode = ''
        # get the notes here
        for i in range(len(self.ref_para_list)):
            secondCode = secondCode + self.ref_para_list[i][0] + ',' + self.ref_para_list[i][1] + '|'
            secondCode = secondCode.replace('\n', '')
        #------------------done getting the arguments------------------------------

        # get the program path
        program = '"' + GlobalSettings.appdir + 'Casper_Seq_Finder' + '" '

        # compile all of the arguments into one line
        args = '"' + endo_choice + '" '
        args = args + '"' + pam + '" '
        args = args + '"' + code + '" '
        args = args + str(pamdir) + ' '
        args = args + '"' + output_location + '/" '
        args = args + '"' + path_to_info + '" '
        args = args + '"' + path_to_fna + '" '
        args = args + '"' + orgName + '" '
        args = args + gRNA_length + ' '
        args = args + seed_length + ' '
        # combine the program and arguments into 1
        program = program + args

        self.process.readyReadStandardOutput.connect(partial(output_stdout, self.process))
        self.proc_running = True
        self.process.start(program)
        self.process.finished.connect(finished)