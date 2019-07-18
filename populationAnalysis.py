from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import os
from CSPRparser import CSPRparser
import Algorithms

class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi('populationanalysis.ui', self)
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.goBackButton.clicked.connect(self.go_back)
        self.analyze_button.clicked.connect(self.fill_data)
        self.clear_Button.clicked.connect(self.clear)
        self.parser = CSPRparser("")
        self.Endos = dict()
        self.cspr_files = {}
        self.sq=Algorithms.SeqTranslate()


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
        self.table2.setColumnCount(7)
        self.table2.setShowGrid(False)
        self.table2.setHorizontalHeaderLabels(["Sequence","Strand","PAM", "Score","% Conserved","Total Repeats","Avg. Repeats/Organism"])
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

    # launches the chromesome selector window
    def launch_chrom_selector(self):
        GlobalSettings.mainWindow.cspr_selector.launch(self.cspr_files)

    def launch(self,path):
        os.chdir(path)
        self.directory = path
        self.get_data()

    def get_data(self):
        onlyfiles = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
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
        self.data = orgsandendos
        self.shortHand = shortName
        self.fillEndo()
        self.changeEndos()

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
                    default_seed_length = line_tokened[2].split(",")[0]
                    default_tot_length = line_tokened[2].split(",")[1]
                    self.Endos[endo + " PAM: " + p_pam] = (endo, p_pam, default_seed_length, default_tot_length)

                break
        f.close()
        self.endoBox.addItem("None Selected")
        self.endoBox.addItems(self.Endos.keys())
        self.endoBox.currentIndexChanged.connect(self.changeEndos)

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

    def fill_data(self):

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
            self.parser.popParser(filenames)
            #print(self.parser.popData.keys())
            self.table2.setRowCount(0)
            index = 0
            for files in selected_files:
                self.parser.fileName = self.cspr_files[files.text()]
                self.parser.read_chromesome()
                for items in self.parser.chromesomeList:
                    first = True
                    for data in items:
                        if first == True:
                            first = False
                        else:
                            self.table2.setRowCount(index + 1)
                            seq = QtWidgets.QTableWidgetItem(data[1])
                            strand = QtWidgets.QTableWidgetItem(str(data[4]))
                            PAM = QtWidgets.QTableWidgetItem(data[2])
                            num1 = int(data[3])
                            score = QtWidgets.QTableWidgetItem()
                            score.setData(QtCore.Qt.EditRole, num1)

                            if data[1] in self.parser.popData[files.text()]:
                                print(self.parser.popData[files.text()][data[1]])


                            self.table2.setItem(index, 0, seq)
                            self.table2.setItem(index, 1, strand)
                            self.table2.setItem(index, 2, PAM)
                            self.table2.setItem(index, 3, score)

                            index += 1
                self.table2.resizeColumnsToContents()
            self.plot_repeats_vs_seeds()
        else:
            msg = QtWidgets.QMessageBox()
            msg.setWindowTitle("Error")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("<font size=4>" + "Endo does not match." + "</font>")
            msg.exec()

    def clear(self):
        self.table2.setRowCount(0)

    def plot_repeats_vs_seeds(self):
        selected_files = self.org_Table.selectedItems()
        first = True
        cspr_filenames = []
        plots = []
        i = 0
        for files in selected_files:
            self.parser.fileName = self.cspr_files[files.text()]
            cspr_filenames.append(files.text())
            self.parser.read_repeats()
            index = 0
            data = {}
            max = 0
            y1 = []
            x1 = []
            for seed in self.parser.repeats:
                number = self.parser.repeats[seed]
                if number in data:
                    data[number]+=1
                else:
                    data[number] =1

            for number in self.order(data):

                if int(data[number]) >max:
                    max = int(data[number])
                    self.mode = number

                hold = 0
                while hold < data[number]:
                    if index == int(round(len(self.parser.repeats) / 2)):
                        self.median = number
                    x1.append(index)
                    y1.append(number)
                    index= index+1
                    hold +=1

            if(first == True):
                first = False
                #clear axes
                self.pop_analysis_repeats_graph.canvas.axes.clear()
                #the following are for plotting / formatting
                self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1, label=cspr_filenames[i])
                self.pop_analysis_repeats_graph.canvas.axes.set_xlabel('Seed Id Number')
                self.pop_analysis_repeats_graph.canvas.axes.set_ylabel('Number of Repeats')
                self.pop_analysis_repeats_graph.canvas.axes.set_title('Number of Repeats per Seed Id Number')
                #always redraw at the end
                self.pop_analysis_repeats_graph.canvas.draw()
                i += 1
            else:
                self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1, label=cspr_filenames[i])
                self.pop_analysis_repeats_graph.canvas.draw()
                i += 1

        self.pop_analysis_repeats_graph.canvas.axes.legend(loc=0)
        self.pop_analysis_repeats_graph.canvas.draw()

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
        GlobalSettings.mainWindow.show()
        self.hide()

    # this function calls the close window class. Allows the user to choose what files they want to keep/delete
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()