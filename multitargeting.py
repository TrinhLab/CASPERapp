from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import GlobalSettings
from PyQt5.QtChart import QChartView
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
from matplotlib.ticker import MaxNLocator
import os
import sqlite3
import gzip
from collections import Counter
import statistics
import repeats_vs_seeds_line


class Multitargeting(QtWidgets.QMainWindow):


    def __init__(self, parent = None):
        super(Multitargeting, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'multitargetingwindow.ui', self)
        self.setWindowIcon(QtGui.QIcon(GlobalSettings.appdir + "cas9image.png"))
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
                        font: 15pt "Arial";
                        font: bold;
                        margin-top: 10px;}"""
        self.groupBox.setStyleSheet(groupbox_style)
        self.groupBox_2.setStyleSheet(groupbox_style.replace("groupBox","groupBox_2"))
        self.groupBox_3.setStyleSheet(groupbox_style.replace("groupBox","groupBox_3"))

        # Initializes layouts for the graphs
        self.global_line = QtWidgets.QVBoxLayout()
        self.global_bar = QtWidgets.QVBoxLayout()
        self.seed_bar = QtWidgets.QVBoxLayout()
        self.global_line.setContentsMargins(0,0,0,0)
        self.global_bar.setContentsMargins(0,0,0,0)
        self.seed_bar.setContentsMargins(0,0,0,0)

        self.data = ""
        self.shortHand = ""
        self.chromo_length = list()

        # Listeners for changing the seed sequence or the .cspr file
#        self.chromo_seed.currentIndexChanged.connect(self.seed_chromo_changed)
#        self.update_min_max.clicked.connect(self.update)
        self.Analyze_Button.clicked.connect(self.make_graphs)
        self.statistics_overview.clicked.connect(self.show_statistics)

        # go back to main button
        self.back_button.clicked.connect(self.go_back)

        # Tool Bar options
        self.actionCASPER.triggered.connect(self.changeto_main)

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
        self.directory = 'Cspr files'
        self.info_path = os.getcwd()

        ##################################
        self.scene = QtWidgets.QGraphicsScene()
        self.graphicsView.setScene(self.scene)
        self.scene2 = QtWidgets.QGraphicsScene()
        self.graphicsView_2.setScene(self.scene2)
        self.graphicsView.viewport().installEventFilter(self)

        self.loading_window = loading_window()
        screen = QtGui.QGuiApplication.screenAt(QtGui.QCursor().pos())

        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.mwfg.moveCenter(self.cp)  ##Center window
        self.move(self.mwfg.topLeft())  ##Center window
        self.hide()

    def show_statistics(self):
        if (self.line_bool and self.bar_bool):
            self.multitargeting_statistics.show()
        else:
            QtWidgets.QMessageBox.question(self, "No analysis run.", 'Multitargeting Analysis must be performed before viewing statistics.\n\nSelect an organism and endonuclease and click "Analyze" then try again.', QtWidgets.QMessageBox.Ok)
            return True

    def eventFilter(self, source, event):
        if (event.type() == QtCore.QEvent.MouseMove and source is self.graphicsView.viewport()):
            coord = self.graphicsView.mapToScene(event.pos())
            for i in self.bar_coords:
                x = i[1]
                y1 = i[2]
                y2 = i[3]
                dups = 0
                if ((coord.x() == x or coord.x() == x + 1 or coord.x() == x - 1) and (
                        coord.y() >= y1 and coord.y() <= y2)):

                    listtemp = []
                    for a in self.bar_coords:
                        if (x == a[1] and y1 == a[2] and y2 == a[3]):
                            listtemp.append(a)
                            dups += 1
                    self.scene2 = QtWidgets.QGraphicsScene()
                    self.graphicsView_2.setScene(self.scene2)
                    output = str()
                    i = 1
                    for item in listtemp:
                        ind = item[0]
                        temp = self.event_data[ind]
                        if len(listtemp) > 1 and i < len(listtemp):
                            output += 'Location: ' + str(temp[0]) + ' | Seq: ' + str(temp[1]) + ' | PAM: ' + str(
                                temp[2]) + ' | SCR: ' + str(temp[3]) + ' | DIRA: ' + str(temp[4]) + '\n'
                        else:
                            output += 'Location: ' + str(temp[0]) + ' | Seq: ' + str(temp[1]) + ' | PAM: ' + str(
                                temp[2]) + ' | SCR: ' + str(temp[3]) + ' | DIRA: ' + str(temp[4])
                        i += 1
                    text = self.scene2.addText(output)
                    font = QtGui.QFont()
                    font.setBold(True)
                    font.setPointSize(12)
                    text.setFont(font)

        return Qt.QWidget.eventFilter(self, source, event)


    def launch(self,path):
        os.chdir(path)
        self.directory = path
        self.get_data()


    def get_data(self):
        #disconnect index changed signal on endo dropdown if there is one
        try:
            self.endo_drop.currentIndexChanged.disconnect()
        except:
            pass

        onlyfiles = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        self.organisms_to_files = {}
        self.organisms_to_endos = {}
        # shortName = {}
        self.endo_drop.clear()
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
                    if self.organism_drop.findText(species) == -1:
                        self.organism_drop.addItem(species)

        #fill in endos dropdown based on current organism
        endos = self.organisms_to_endos[str(self.organism_drop.currentText())]
        self.endo_drop.addItems(endos)
        self.organism_drop.currentIndexChanged.connect(self.update_endos)
        self.endo_drop.currentIndexChanged.connect(self.change_endos)

        #update file names for current org/endo combo
        self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][endos[0]][0]
        self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][endos[0]][1]

    def change_endos(self):
        #update file names based on current org/endo combo
        self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][str(self.endo_drop.currentText())][0]
        self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][str(self.endo_drop.currentText())][1]

    def update_endos(self):
        #try to disconnect index changed signal on endo dropdown if there is one
        try:
            self.endo_drop.currentIndexChanged.disconnect()
        except:
            pass

        #clear endo dropdown and fill in with endos relative to the current organism
        self.endo_drop.clear()
        endos = self.organisms_to_endos[str(self.organism_drop.currentText())]
        self.endo_drop.addItems(endos)
        self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][endos[0]][0]
        self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][endos[0]][1]

        #reconnect index changed signal on endo dropdown
        self.endo_drop.currentIndexChanged.connect(self.change_endos)


    def make_graphs(self):
        # self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][0]
        # self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][0]

        self.loading_window.loading_bar.setValue(0)
        self.loading_window.show()
        QtCore.QCoreApplication.processEvents()
        self.chromo_length.clear()
        self.loading_window.loading_bar.setValue(5)
        self.plot_repeats_vs_seeds()
        self.loading_window.loading_bar.setValue(20)
        self.bar_seeds_vs_repeats()
        self.loading_window.loading_bar.setValue(40)
        #self.fill_min_max()
        self.loading_window.loading_bar.setValue(60)
        #self.fill_seed_id_chrom()
        self.loading_window.loading_bar.setValue(80)
        #self.fill_Chromo_Text(self.chromo_seed.currentText())
        self.loading_window.loading_bar.setValue(100)
        self.multitargeting_statistics.avg_rep.setText(str(round(float(self.average),1)))
        self.multitargeting_statistics.med_rep.setText(str(round(float(self.median),1)))
        self.multitargeting_statistics.mode_rep.setText(str(round(float(self.mode),1)))
        self.multitargeting_statistics.nbr_seq.setText(str(round(float(self.repeat_count),1)))
        self.loading_window.hide()


    #fill in chromo bar visualization
    def fill_Chromo_Text(self, seed):
        # self.cspr_file = self.organisms_to_files[str(self.organism_drop.currentText())][0]
        # self.db_file = self.organisms_to_files[str(self.organism_drop.currentText())][0]

        chromo_pos = {}
        # get kstats
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        kstats = []
        with gzip.open(self.cspr_file, "r") as f:
            for line in f:
                buf = str(line)
                buf = buf.strip("'b")
                buf = buf[:len(buf) - 2]
                if buf.find("KARYSTATS") != -1:
                    buf = buf.replace("KARYSTATS: ", "")
                    kstats = buf.split(',')
                    kstats = kstats[:-1]

                    break

        seed = self.chromo_seed.currentText() ###Change to selection from table
        data = c.execute("SELECT chromosome, location FROM repeats WHERE seed = ? ", (seed,)).fetchone()
        c.close()
        if data != None:
            data = list(data)
            chromo = data[0].split(',')
            pos = data[1].split(',')
            self.event_data = {}
            for i in range(len(chromo)):
                c = int(chromo[i])
                p = abs(int(pos[i]))
                k = int(kstats[c - 1])
                new_pos = int((p / k) * 485)
                if c in chromo_pos.keys():
                    chromo_pos[c].append(new_pos)
                else:
                    chromo_pos[c] = [new_pos]
            i = 0
            self.scene = QtWidgets.QGraphicsScene()
            self.graphicsView.setScene(self.scene)
            self.bar_coords.clear()  # clear bar_coords list before creating visual
            ind = 0
            for chromo in chromo_pos:
                pen_blk = QtGui.QPen(QtCore.Qt.black)
                pen_red = QtGui.QPen(QtCore.Qt.red)
                pen_blk.setWidth(3)
                pen_red.setWidth(3)
                if i == 0:
                    text = self.scene.addText(str(chromo))
                    text.setPos(0, 0)
                    font = QtGui.QFont()
                    font.setBold(True)
                    font.setPointSize(10)
                    text.setFont(font)
                    self.scene.addRect(40, (i * 25), 525, 25, pen_blk)
                else:
                    text = self.scene.addText(str(chromo))
                    font = QtGui.QFont()
                    font.setBold(True)
                    font.setPointSize(10)
                    text.setFont(font)
                    text.setPos(0, i * 25 + 10 * i)
                    self.scene.addRect(40, (i * 25) + 10 * i, 525, 25, pen_blk)
                for k in chromo_pos[chromo]:
                    line = self.scene.addLine(k + 40, (i * 25) + 3 + 10 * i, k + 40, (i * 25) + 22 + 10 * i, pen_red)
                    temp = []  # used for storing coordinates and saving them in self.bar_coords[]
                    temp.append(ind)  # index value
                    temp.append(k + 40)  # x value
                    temp.append((i * 25) + 3 + 10 * i)  # y1
                    temp.append((i * 25) + 22 + 10 * i)  # y2
                    self.bar_coords.append(temp)  # push x, y1, and y2 to this list
                    ind += 1
                i = i + 1
            self.generate_event_data()
            return True
        else:
            QtWidgets.QMessageBox.information(self, "Seed Error",
                                           "No such seed exists in the repeats section of this organism.",
                                           QtWidgets.QMessageBox.Ok)
            return False


    # get data for chromsome viewer to display
    def generate_event_data(self):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
#        seed = self.chromo_seed.currentText()
        data = list(c.execute("SELECT * FROM repeats WHERE seed = ? ", (seed,)).fetchone())
        c.close()
        chromo = data[1].split(',')
        loc = data[2].split(',')
        three = data[3].split(',')
        five = data[4].split(',')
        pam = data[5].split(',')
        score = data[6].split(',')
        five_bool = True
        three_bool = True
        if five[0] == '':
            five_bool = False
        if three[0] == '':
            three_bool = False

        self.event_data = {}
        for i in range(len(chromo)):
            if int(loc[i]) < 0:
                dira = '-'
            else:
                dira = '+'

            if five_bool and not three_bool:
                seq = five[i] + seed
            elif not five_bool and three_bool:
                seq = seed + three[i]
            else:
                seq = five[i] + seed + three[i]

            self.event_data[i] = [str(abs(int(loc[i]))), seq, pam[i], score[i], dira]


    # creates bar graph num of repeats vs. chromsome
    # this graphs is connected to the repeats_vs_chromo.py file
    # to represent the widget space in the UI file
    def chro_bar_create(self, seed):
        ###Clear out old widgets in layout
        for i in reversed(range(self.seed_bar.count())): 
            self.seed_bar.itemAt(i).widget().setParent(None)
        self.seed_canvas = MplCanvas(self, width=5, height=4, dpi=100) ###Initialize new Canvas
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
        self.seed_canvas.axes.set_xlabel('Chromosome', fontsize = 10)
        self.seed_canvas.axes.set_ylabel('Number of Repeats', fontsize=10)
        self.line_canvas.axes.set_title('Repeats per Scaffold/Chromosome',fontsize=10)
        self.line_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
        self.line_canvas.draw()


    # plots the sequences per Number Repeats bar graph
    # this graph is connected to the seeds_vs_repeats_bar.py file
    # to represent the wdiget space in the UI file
    def bar_seeds_vs_repeats(self):
        ###Clear out old widgets in layout
        for i in reversed(range(self.global_bar.count())): 
            self.global_bar.itemAt(i).widget().setParent(None)
        self.bar_canvas = MplCanvas(self, width=5, height=4, dpi=100) ###Initialize new Canvas
        self.global_bar.addWidget(self.bar_canvas) ### Add canvas to global line layout 
        self.seeds_vs_repeats_bar.setLayout(self.global_bar) ### Add global line layout to repeats vs. seeds line plot widget

        """ Get the data """ 
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT seed, count from repeats").fetchall()
        c.close()
        repeat_data = []
        for obj in data:
            cnt = list(obj)[1]
            repeat_data.append(cnt)
        repeat_data = Counter(repeat_data)
        max_rep = max(repeat_data.values())
        x_labels = []
        y = []
        for rep in sorted(repeat_data, key=repeat_data.get, reverse=True):
            #if repeat_data[rep] / max_rep > 0.01:
            x_labels.append(rep)
            y.append(repeat_data[rep])
        x = list(range(0, len(x_labels)))
        # the following are plotting / formatting for the graph
        self.bar_canvas.axes.bar(x, y)
        self.bar_canvas.axes.set_xticks(x)
        self.bar_canvas.axes.set_xticklabels(x_labels)
        self.bar_canvas.axes.set_xlabel('Number of Repeats', fontsize=10)
        self.bar_canvas.axes.set_ylabel('Number of Sequences', fontsize=10)
        self.bar_canvas.axes.set_title('Number of Sequences per Number of Repeats',fontsize=10)
        self.bar_canvas.axes.tick_params(axis='both', which='major', labelsize=8)

        # rects are all the bar objects in the graph
        rects = self.bar_canvas.axes.patches
        rect_vals = []
        # this for loop will calculate the height and create an annotation for each bar
        for rect in rects:
            height = rect.get_height()
            temp = self.bar_canvas.axes.text(rect.get_x() + rect.get_width() / 2, height,
                                                              '%d' % int(height),
                                                              ha='center', va='bottom')
            temp.set_visible(False)
            rect_vals.append(temp)

        # function used for when user cursor is hovering over the bar, if hovering over a bar, the
        # height annotatin will appear above the bar, otherwise it will be hidden
        def on_plot_hover(event):
            i = 0
            for rect in rects:
                if rect.contains(event)[0]:
                    rect_vals[i].set_visible(True)
                else:
                    rect_vals[i].set_visible(False)
                i = i + 1
            self.bar_canvas.draw()

        # statement to detect cursor hovering over the bars
        self.bar_canvas.mpl_connect('motion_notify_event', on_plot_hover)
        # must redraw after every change
        self.bar_canvas.draw()


    # plots the repeats per ID number graph as line graph
    # this graph is connected to the repeats_vs_seeds_line.py file
    
    def plot_repeats_vs_seeds(self):
        ###Clear out old widgets in layout
        for i in reversed(range(self.global_line.count())): 
            self.global_line.itemAt(i).widget().setParent(None)
        self.line_canvas = MplCanvas(self, width=5, height=4, dpi=100) ###Initialize new Canvas
        self.global_line.addWidget(self.line_canvas) ### Add canvas to global line layout 
        self.repeats_vs_seeds_line.setLayout(self.global_line) ### Add global line layout to repeats vs. seeds line plot widget

        """ Fetch all the data """
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT seed, count from repeats").fetchall()
        c.close()
        x1 = list(range(0, len(data)))
        y1 = []
        for obj in data:
            y1.append(int(list(obj)[1]))

        #y1.sort(reverse=True)

        #get stats
        self.average = statistics.mean(y1)
        self.mode = statistics.mode(y1)
        self.median = statistics.median(y1)
        self.repeat_count = len(data)

        # clear axes
        self.line_canvas.axes.clear()
        #Plotting / formatting
        self.line_canvas.axes.plot(x1, y1)
        self.line_canvas.axes.set_xlabel('Seed ID Number',fontsize=10)
        self.line_canvas.axes.set_ylabel('Number of Repeats',fontsize=10)
        self.line_canvas.axes.set_title('Number of Repeats per Seed ID Number',fontsize=10)
        self.line_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
        # always redraw at the end
        self.line_bool = True
        self.line_canvas.draw()


    #fills min and max dropdown windows ###NEED TO EDIT FOR TEXTEDIT
    def fill_min_max(self,run_seed_fill=True):
#        self.update_min_max.clicked.disconnect()
        self.max_chromo.clear()
        self.min_chromo.clear()
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT MAX(count) from repeats").fetchone()
        c.close()
        max_rep = int(data[0])
        nums = list(range(1, max_rep + 1))
        self.min_chromo.setText(str(min(nums)))
        self.max_chromo.setText(str(max(nums)))
#########        self.update_min_max.clicked.connect(self.update)


    #fill_seed_id_chrom will fill the seed ID dropdown, and create the chromosome graph
    def fill_seed_id_chrom(self):
#        self.chromo_seed.disconnect()
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        min_v = int(self.min_chromo.toPlainText())
        max_v = int(self.max_chromo.toPlainText())
        data = c.execute("SELECT seed, count from repeats where count >= ?  and count <= ? LIMIT 0,1000", (min_v, max_v,)).fetchall()
        c.close()
        model = QtGui.QStandardItemModel()
        parentItem = model.invisibleRootItem()
        i = 0
        self.seeds = []
        for obj in data:
            seed = (list(obj)[0])
            if i > 999:
                break
            else:
                if i == 0:
                    self.chro_bar_create(seed)
                item = QtGui.QStandardItem(seed)
                parentItem.appendRow(item)
                self.seeds.append(seed)
            i += 1
        self.seed_counter = 0
        self.seed_max_counter = len(self.seeds) / 1000
#        self.chromo_seed.setModel(model)
#        self.chromo_seed.currentIndexChanged.connect(self.seed_chromo_changed)


    def update(self):
        self.scene2 = QtWidgets.QGraphicsScene()
        self.graphicsView_2.setScene(self.scene2)
        if self.min_chromo.toPlainText() != '' and self.max_chromo.toPlainText() != '':
            if int(self.min_chromo.toPlainText()) > int(self.max_chromo.toPlainText()):
                QtWidgets.QMessageBox.question(self, "Maximum cant be less than Minimum",
                                               "The Minimum number of repeats cant be more than the Maximum",
                                               QtWidgets.QMessageBox.Ok)
            else:
                self.loading_window.loading_bar.setValue(5)
                self.loading_window.show()
                QtCore.QCoreApplication.processEvents()
                if self.seed.toPlainText() == '':
                    self.fill_seed_id_chrom()
                    self.loading_window.loading_bar.setValue(50)
#                    self.fill_Chromo_Text(self.chromo_seed.currentText())
                    self.loading_window.loading_bar.setValue(100)
                else:
                    result = self.fill_Chromo_Text(self.seed.toPlainText())
                    self.loading_window.loading_bar.setValue(50)
                    if result == True:
                        self.chro_bar_create(self.seed.toPlainText())
                    self.loading_window.loading_bar.setValue(100)
        self.loading_window.hide()


    def seed_chromo_changed(self):
        self.loading_window.loading_bar.setValue(5)
        self.loading_window.show()
        QtCore.QCoreApplication.processEvents()
        self.seed.setText('')
#        self.chro_bar_create(self.chromo_seed.currentText())
        self.loading_window.loading_bar.setValue(50)
#        self.fill_Chromo_Text(self.chromo_seed.currentText())
        self.loading_window.loading_bar.setValue(100)
        self.loading_window.hide()


    #connects to view->CASPER to switch back to the main CASPER window
    def changeto_main(self):
        GlobalSettings.mainWindow.mwfg.moveCenter(GlobalSettings.mainWindow.cp)  ##Center window
        GlobalSettings.mainWindow.move(GlobalSettings.mainWindow.mwfg.topLeft())  ##Center window
        GlobalSettings.mainWindow.show()
        self.hide()


    #connects to go back button in bottom left to switch back to the main CASPER window
    def go_back(self):
        GlobalSettings.mainWindow.mwfg.moveCenter(GlobalSettings.mainWindow.cp)  ##Center window
        GlobalSettings.mainWindow.move(GlobalSettings.mainWindow.mwfg.topLeft())  ##Center window
        GlobalSettings.mainWindow.show()
        self.hide()

    # this function calls the closingWindow class.
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

class Multitargeting_Statistics(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(Multitargeting_Statistics, self).__init__(parent)
        uic.loadUi(GlobalSettings.appdir + 'multitargeting_statistics.ui', self)
        self.hide()