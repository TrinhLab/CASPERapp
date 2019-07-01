from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import os

class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi('populationanalysis.ui', self)
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.goBackButton.clicked.connect(self.go_back)
        self.Endos = dict()

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
        self.table2.setColumnCount(6)
        self.table2.setShowGrid(False)
        self.table2.setHorizontalHeaderLabels(["Sequence","Strand","PAM","% Conserved","Total Repeats","Repeats/Organism"])
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
                    self.org_Table.setRowCount(index + 1)
                    orgsandendos[species] =[endo]
                    name = QtWidgets.QTableWidgetItem(str(species))
                    self.org_Table.setItem(index, 0, name)
                    index+=1
        self.org_Table.resizeColumnsToContents()
        self.data = orgsandendos
        self.shortHand = shortName
        self.fillEndo()

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
                    self.Endos[endo + "PAM: " + p_pam] = (endo, p_pam, default_seed_length, default_tot_length)

                break
        f.close()
        self.endoBox.addItems(self.Endos.keys())

    def go_back(self):
        GlobalSettings.mainWindow.show()
        self.hide()