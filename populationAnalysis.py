from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import os
from CSPRparser import CSPRparser
import Algorithms
from functools import partial

from statistics import mode

class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi('populationanalysis.ui', self)
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.goBackButton.clicked.connect(self.go_back)
        self.analyze_button.clicked.connect(self.fill_data)
        self.clear_Button.clicked.connect(self.clear)
        self.ncbi_search_button.clicked.connect(self.launch_ncbi_seacher)
        self.meta_genomic_cspr_checkbox.stateChanged.connect(self.get_data)
        self.parser = CSPRparser("")
        self.Endos = dict()
        self.fna_files = dict()
        self.cspr_files = {}
        self.sq=Algorithms.SeqTranslate()
        self.ref_para_list = list()


        #orgonaism table
        self.org_Table.setColumnCount(1)
        self.org_Table.setShowGrid(False)
        self.org_Table.setHorizontalHeaderLabels(["Organism"])
        self.org_Table.horizontalHeader().setSectionsClickable(True)
        self.org_Table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.org_Table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.org_Table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.org_Table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)

        #top right table
        self.table2.setColumnCount(9)
        self.table2.setShowGrid(False)
        self.table2.setHorizontalHeaderLabels(["Seed","% Conserved","Total Repeats","Avg. Repeats/Organism", "Consensus Sequence", "% Consensus", "Strand","PAM", "Score"])
        self.table2.horizontalHeader().setSectionsClickable(True)
        self.table2.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table2.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.table2.resizeColumnsToContents()

        #Finder table
        self.loc_finder_table.setColumnCount(4)
        self.loc_finder_table.setShowGrid(False)
        self.loc_finder_table.setHorizontalHeaderLabels(["Sequence", "Organism", "Chromosome", "Location"])
        self.loc_finder_table.horizontalHeader().setSectionsClickable(True)
        self.loc_finder_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.loc_finder_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.loc_finder_table.resizeColumnsToContents()

        # action buttons
        self.actionMetaGenome_Parser.triggered.connect(self.launch_chrom_selector)

        self.combinerWindow = fna_and_cspr_combiner()

    def launch_ncbi_seacher(self):
        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(0)
        GlobalSettings.mainWindow.ncbi_search_dialog.show()

    # launches the chromesome selector window
    def launch_chrom_selector(self):
        GlobalSettings.mainWindow.cspr_selector.launch(self.cspr_files)

    def launch(self,path):
        os.chdir(path)
        self.directory = path
        self.get_data()

    # this function builds the Select Organisms table
    def get_data(self):
        onlyfiles = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        self.fna_files.clear()

        # show/hide the stuff that isn't needed 
        if self.meta_genomic_cspr_checkbox.isChecked():
            self.endoBox.hide()
            self.ncbi_search_button.hide()
            self.label_3.hide()
            self.label_2.setText('Select 1 Meta Genomic CSPR File')
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
                    f = open(file, 'r')
                    hold = f.readline()
                    f.close()
                    # only show the files that are metagenomic
                    if '_meta' in hold:
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
        
        self.org_Table.resizeColumnsToContents()
        """
        orgsandendos = {}
        shortName = {}
        index = 0


        for file in onlyfiles:
            if file.find('.cspr')!=-1:
                newname = file[0:-4]
                s = newname.split('_')
                hold = open(file)
                buf = (hold.readline())
                species = buf[8:buf.find('\n')]
                endo = str(s[1])
                if species not in shortName:
                    shortName[species] = s[0]
                if species in orgsandendos:
                    orgsandendos[species].append(endo)
                else:
                    orgsandendos[species] =[endo]
                    self.cspr_files[str(species)] = file
                    index+=1

        # data is a dict. Key is the organism (taken from the cspr file) key is a list of endonuclease's that the organism has or that the user has
        self.data = orgsandendos
        # shortHand is another dict where the key is the organism name (Taken from the cspr file). The value the short hand for that org (bsu, sce, yli)
        self.shortHand = shortName
        """
        self.fillEndo()
        #self.changeEndos()

    # this function opens CASPERinfo and builds the dropdown menu of selectable endonucleases
    def fillEndo(self):
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            f = open(GlobalSettings.appdir + "\\CASPERinfo")
        else:
            f = open(GlobalSettings.appdir + "/CASPERinfo")
        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if(line[0]=="-"):
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
        #self.endoBox.currentIndexChanged.connect(self.changeEndos)

    # this function displays all of the organisms of which the user has that endo in their DB
    def changeEndos(self):
        endo_box = str(self.endoBox.currentText())
        endo_box = endo_box[:endo_box.find(" ")]
        self.org_Table.setRowCount(0)
        index = 0
        for keys in self.cspr_files.keys():
            filename = str(self.cspr_files[keys])
            endo = filename[filename.find("_")+1:]
            endo = endo.replace(".cspr","")
            if(endo == endo_box or endo_box == "None"):
                self.org_Table.setRowCount(index + 1)
                name = QtWidgets.QTableWidgetItem(str(keys))
                self.org_Table.setItem(index, 0, name)
                index += 1
        self.org_Table.resizeColumnsToContents()

    # this function calls the popParser function and fills all the tables
    def fill_data(self):
        selectedList = self.org_Table.selectedItems()

        # if the table is showing only fna/fasta files
        if not self.meta_genomic_cspr_checkbox.isChecked():
            # rules for selecting FNA/Fasta files
            # check to make sure that the user selected at least 2 organisms, and 1 endonuclease
            if len(selectedList) < 1 or self.endoBox.currentText() == 'None Selected':
                QtWidgets.QMessageBox.question(self, "Nothing Seleted", "No items selected. Please select at least 1 organism, and only 1 endonuclease",
                                            QtWidgets.QMessageBox.Ok)
                return
            if len(selectedList) == 1:
                error = QtWidgets.QMessageBox.question(self, "Only 1 Organism Selected",
                                                        "Population Analysis works with multiple organisms, or a meta genome. If the file selected it not a meta genome, the program may not function correctly. Do you wish to continu?.\n\n"
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
                QtWidgets.QMessageBox.question(self, "Nothing Seleted", "No items selected. Please select one meta genome for Population Analysis.",
                                            QtWidgets.QMessageBox.Ok)
                return
            elif len(selectedList) > 1:
                QtWidgets.QMessageBox.question(self, "Too many Selected", "Only 1 meta genomic CSPR file is allowed to be selected",
                                            QtWidgets.QMessageBox.Ok)
                return
    """
        endo = str(self.endoBox.currentText())
        endo = endo[:endo.find(" ")]
        selected_files = self.org_Table.selectedItems()
        error = False
        filenames = []
        for files in selected_files:
            cspr_name = self.cspr_files[files.text()]
            cspr_name = str(cspr_name)
            cspr_name = cspr_name[cspr_name.find('_')+1:]
            cspr_name = cspr_name[:cspr_name.find('.')]
            current_endo = str(self.endoBox.currentText())
            if(cspr_name !=  current_endo[:current_endo.find(" ")]):
                error = True
            filenames.append(self.cspr_files[files.text()])



        if(error != True):
            self.parser.popParser(filenames, endoChoice=endo)
            print(self.parser.popData)
            self.table2.setRowCount(0)
            index = 0
            for seeds in self.parser.popData:
                self.table2.setRowCount(index + 1)

                seed = QtWidgets.QTableWidgetItem()
                total_repeats = QtWidgets.QTableWidgetItem()
                total_repeats.setData(QtCore.Qt.EditRole, len(self.parser.popData[seeds]))
                seed.setData(QtCore.Qt.EditRole, str(seeds))

                self.table2.setItem(index, 0, seed)
                self.table2.setItem(index, 2, total_repeats)


                #loop through the tuples for each seed
                sequences = []
                for tuples in self.parser.popData[seeds]:
                    sequences.append(tuples[3])

                con_seq_temp = str(max(set(sequences), key=sequences.count))
                consensus_seq = QtWidgets.QTableWidgetItem()
                consensus_seq.setData(QtCore.Qt.EditRole, con_seq_temp)
                self.table2.setItem(index, 4, consensus_seq)

                consensus_percentage = sequences.count(con_seq_temp) / len(sequences) * 100
                consensus_percentage = round(consensus_percentage, 1)
                consensus_perc = QtWidgets.QTableWidgetItem()
                consensus_perc.setData(QtCore.Qt.EditRole, consensus_percentage)
                self.table2.setItem(index, 5, consensus_perc)


                index += 1
            self.table2.resizeColumnsToContents()


            # for files in selected_files:
            #     self.parser.fileName = self.cspr_files[files.text()]
            #     self.parser.read_chromesome()
            #     for items in self.parser.chromesomeList:
            #         first = True
            #         for data in items:
            #             if first == True:
            #                 first = False
            #             else:
            #                 self.table2.setRowCount(index + 1)
            #                 seq = QtWidgets.QTableWidgetItem(data[1])
            #                 strand = QtWidgets.QTableWidgetItem(str(data[4]))
            #                 PAM = QtWidgets.QTableWidgetItem(data[2])
            #                 num1 = int(data[3])
            #                 score = QtWidgets.QTableWidgetItem()
            #                 score.setData(QtCore.Qt.EditRole, num1)
            #
            #                 #if data[1] in self.parser.popData[files.text()]:
            #                  #   print(self.parser.popData[files.text()][data[1]])
            #
            #
            #                 self.table2.setItem(index, 0, seq)
            #                 self.table2.setItem(index, 1, strand)
            #                 self.table2.setItem(index, 2, PAM)
            #                 self.table2.setItem(index, 3, score)
            #
            #                 index += 1
            #     self.table2.resizeColumnsToContents()
            # self.plot_repeats_vs_seeds()
        else:
            msg = QtWidgets.QMessageBox()
            msg.setWindowTitle("Error")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("<font size=4>" + "Endo does not match." + "</font>")
            msg.exec()
    """
    def clear(self):
        self.table2.setRowCount(0)

    def plot_repeats_vs_seeds(self):
        print('graph')
        # selected_files = self.org_Table.selectedItems()
        # first = True
        # cspr_filenames = []
        # plots = []
        # i = 0
        # endo = str(self.endoBox.currentText())
        # endo = endo[:endo.find(" ")]
        # print(endo)
        # for files in selected_files:
        #     self.parser.fileName = self.cspr_files[files.text()]
        #     cspr_filenames.append(files.text())
        #     #self.parser.read_repeats(endoChoice=endo)
        #     self.parser.popParser(endoChoice=endo)
        #     index = 0
        #     data = {}
        #     max = 0
        #     y1 = []
        #     x1 = []
        #
        #     for seed in self.parser.popData:
        #
        #
        #
        #     # for seed in self.parser.repeats:
        #     #     number = self.parser.repeats[seed]
        #     #     if number in data:
        #     #         data[number]+=1
        #     #     else:
        #     #         data[number] =1
        #     #
        #     # for number in self.order(data):
        #     #
        #     #     if int(data[number]) >max:
        #     #         max = int(data[number])
        #     #         self.mode = number
        #     #
        #     #     hold = 0
        #     #     while hold < data[number]:
        #     #         if index == int(round(len(self.parser.repeats) / 2)):
        #     #             self.median = number
        #     #         x1.append(index)
        #     #         y1.append(number)
        #     #         index= index+1
        #     #         hold +=1
        #
        #     if(first == True):
        #         first = False
        #         #clear axes
        #         self.pop_analysis_repeats_graph.canvas.axes.clear()
        #         #the following are for plotting / formatting
        #         self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1, label=cspr_filenames[i])
        #         self.pop_analysis_repeats_graph.canvas.axes.set_xlabel('Seed Id Number')
        #         self.pop_analysis_repeats_graph.canvas.axes.set_ylabel('Number of Repeats')
        #         self.pop_analysis_repeats_graph.canvas.axes.set_title('Number of Repeats per Seed Id Number')
        #         #always redraw at the end
        #         self.pop_analysis_repeats_graph.canvas.draw()
        #         i += 1
        #     else:
        #         self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1, label=cspr_filenames[i])
        #         self.pop_analysis_repeats_graph.canvas.draw()
        #         i += 1
        #
        # self.pop_analysis_repeats_graph.canvas.axes.legend(loc=0)
        # self.pop_analysis_repeats_graph.canvas.draw()

    def order(self,data_par):
        data = dict(data_par)
        data2  = []
        while len(data)>0:
            max=0
            for item in data:
                if item>max:
                    max=item
            data2.append(max)
            if len(data2) ==1:
                self.max_repeats =max
            del data[max]
        return data2

    def go_back(self):
        GlobalSettings.mainWindow.getData()
        GlobalSettings.mainWindow.show()
        self.hide()

    # this function calls the close window class. Allows the user to choose what files they want to keep/delete
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()


###############################################################################
# class name: fna_and_cspr_combiner
# this opens a window for the user to select the Organism Name, FNA File Name, and the organism code.
# It combines the FNA files, and then runs the sequencer. This way it creates 1 cspr file.
# after the new cspr file is created, it will run the population analysis
###############################################################################
class fna_and_cspr_combiner(QtWidgets.QDialog):
    def __init__(self):
        # Qt init stuff
        super(fna_and_cspr_combiner, self).__init__()
        uic.loadUi("pop_analysis_fna_combiner.ui", self)
        self.sequencer_prog_bar.setValue(0)

        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.start_button.clicked.connect(self.run_analysis)

        # variables below
        self.ref_para_list = list()
        self.fna_file_names = list()
        self.sq = Algorithms.SeqTranslate()
        self.proc_running = False
        self.generated_files = list()
        self.combined_fna_file = ''

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

        if self.orgName_line_edit.text() == '' or self.org_code_line_edit.text() == '':
            QtWidgets.QMessageBox.question(self, "Missing Information",
                                           "Please input an Organism Name and an Organism Code.",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.combine_fna_files()
        self.build_new_cspr_file()

    # cancel function, just clears everything and closes the window
    def cancel_function(self):
        # check to see if the sequencer is running. If so ask the user if they wish to close out
        if self.proc_running:
            error = QtWidgets.QMessageBox.question(self, "Sequencer is Running",
                                                   "Sequencer is running. Closing this window will cancel that process, and return to the Population Analysis window. .\n\n"
                                                   "Do you wish to continue?",
                                                   QtWidgets.QMessageBox.Yes |
                                                   QtWidgets.QMessageBox.No,
                                                   QtWidgets.QMessageBox.No)
            if (error == QtWidgets.QMessageBox.No):
                return -2
            else:
                self.off_target_running = False
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
            for lines in filter(None, line.split(r'\r\n')):
                if (lines == 'Finished reading in the genome file.'):
                    self.num_chromo_next = True
                elif (self.num_chromo_next == True):
                    self.num_chromo_next = False
                    self.num_chromo = int(lines)
                elif (lines.find('Chromosome') != -1 and lines.find('complete.') != -1):
                    temp = lines
                    temp = temp.replace('Chromosome ', '')
                    temp = temp.replace(' complete.', '')
                    if (int(temp) == self.num_chromo):
                        self.sequencer_prog_bar.setValue(99)
                    else:
                        self.sequencer_prog_bar.setValue(int(temp) / self.num_chromo * 100)
                elif (lines == 'Finished Creating File.'):
                    self.sequencer_prog_bar.setValue(100)
        # this function will end up doing stuff when the process is finished.
        def finished():
            self.proc_running = False

            # get the file name
            cspr_file_name = GlobalSettings.CSPR_DB + '/' + self.org_code_line_edit.text() + '_' + GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0] + '.cspr'
            self.generated_files.append(cspr_file_name)
            os.remove(self.combined_fna_file)
            self.process.kill()
            self.cancel_function()

        #--------------getting the arugments---------------------------------
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
        path_to_info = GlobalSettings.appdir + '/CASPERinfo'
        # make org name something that will make more sense, this is just for testing right now
        orgName = self.orgName_line_edit.text() + '  ,_meta'
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
        program = '"' + GlobalSettings.appdir + '/Casper_Seq_Finder_Windows' + '" '

        # compile all of the arguments into one line
        args =  '"' + endo_choice + '" '
        args = args + '"' + pam + '" '
        args = args + '"' + code + '" '
        args = args + str(pamdir) + ' '
        args = args + '"' + output_location + '" '
        args = args + '"' + path_to_info + '" '
        args = args + '"' + path_to_fna + '" '
        args = args + '"' + orgName + '" '
        args = args + gRNA_length + ' '
        args = args + seed_length + ' '
        args = args + '"' + secondCode + '"'

        # combine the program and arguments into 1
        program = program + args


        self.process.readyReadStandardOutput.connect(partial(output_stdout, self.process))
        self.proc_running = True
        self.process.start(program)
        self.process.finished.connect(finished)