from PyQt6 import QtWidgets, uic, QtGui, QtCore
from PyQt6.QtWidgets import QHeaderView
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import mplcursors
import numpy as np
import matplotlib.patches as patches
from utils.ui import show_error, scale_ui, center_ui

class PopulationAnalysisWindowView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super(PopulationAnalysisWindowView, self).__init__()
        self.global_settings = global_settings
        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(self.global_settings.get_ui_dir() + '/pop.ui', self)
            self.setWindowIcon(QtGui.QIcon(self.global_settings.get_assets_dir() + "/cas9image.ico"))
            self.setWindowTitle('Population Analysis')
            self.setup_tables()
            self.setup_buttons()
            self.setup_colormap()
            self.setup_styles()
            scale_ui(self, base_width=1920, base_height=1080, font_size=12, header_font_size=30)
        except Exception as e:
            show_error(self.global_settings, "Error initializing PopulationAnalysisWindowView.", str(e))

    def setup_tables(self):
        # Organism table
        self.org_Table.setColumnCount(1)
        self.org_Table.setShowGrid(False)
        self.org_Table.setHorizontalHeaderLabels(["Organism"])
        self.org_Table.horizontalHeader().setSectionsClickable(True)
        self.org_Table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.org_Table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.org_Table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.org_Table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)

        # Shared seeds table
        self.table2.setColumnCount(9)
        self.table2.setShowGrid(False)
        self.table2.setHorizontalHeaderLabels(["Seed","% Coverage","Total Repeats","Avg. Repeats/Scaffold", "Consensus Sequence", "% Consensus", "Score","PAM", "Strand"])
        self.table2.horizontalHeader().setSectionsClickable(True)
        self.table2.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table2.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)
        self.table2.setSelectionBehavior(QtWidgets.QTableView.SelectionBehavior.SelectRows)
        self.table2.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
        self.table2.resizeColumnsToContents()

        # Location finder table
        self.loc_finder_table.setColumnCount(5)
        self.loc_finder_table.setShowGrid(False)
        self.loc_finder_table.setHorizontalHeaderLabels(["Seed ID", "Sequence", "Organism", "Scaffold", "Location"])
        self.loc_finder_table.horizontalHeader().setSectionsClickable(True)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.loc_finder_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.loc_finder_table.resizeColumnsToContents()

    def setup_buttons(self):
        self.goBackButton.setText("Go Back")
        self.analyze_button.setText("Analyze")
        self.clear_Button.setText("Clear")
        self.export_button.setText("Export")
        self.find_locs_button.setText("Find Locations")
        self.clear_loc_button.setText("Clear Locations")
        self.query_seed_button.setText("Search Seeds")

    def setup_colormap(self):
        self.colormap_layout = QtWidgets.QVBoxLayout()
        self.colormap_layout.setContentsMargins(0, 0, 0, 0)
        self.colormap_canvas = MplCanvas(self)
        self.colormap_layout.addWidget(self.colormap_canvas)
        self.colormap_figure.setLayout(self.colormap_layout)

    def setup_styles(self):
        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 15px;}
        QGroupBox#groupBox{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        font: bold 14pt 'Arial';
                        margin-top: 10px;}"""
        self.groupBox.setStyleSheet(groupbox_style)
        self.groupBox_2.setStyleSheet(groupbox_style.replace("groupBox","groupBox_2"))

    def get_selected_endo(self):
        return self.endoBox.currentText()

    def get_selected_organisms(self):
        return [index.row() for index in self.org_Table.selectionModel().selectedRows()]

    def get_selected_seeds(self):
        return [self.table2.item(index.row(), 0).text() for index in self.table2.selectionModel().selectedRows()]

    def get_seed_input(self):
        return self.seed_input.text()

    def get_selected_seeds_for_export(self):
        return [item.text() for item in self.table2.selectedItems() if item.column() == 0]

    def update_org_table(self, org_data):
        self.org_Table.setRowCount(len(org_data))
        for row, (org_name, cspr_file, db_file) in enumerate(org_data):
            item = QtWidgets.QTableWidgetItem(org_name)
            item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.org_Table.setItem(row, 0, item)
        self.org_Table.resizeColumnsToContents()

    def update_shared_seeds_table(self, seed_data):
        self.table2.setRowCount(len(seed_data))
        for row, data in enumerate(seed_data):
            for col, value in enumerate(data):
                item = QtWidgets.QTableWidgetItem(str(value))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table2.setItem(row, col, item)
        self.table2.resizeColumnsToContents()

    def update_loc_finder_table(self, loc_data):
        self.loc_finder_table.setRowCount(len(loc_data))
        for row, data in enumerate(loc_data):
            for col, key in enumerate(['seed', 'sequence', 'organism', 'chromosome', 'location']):
                item = QtWidgets.QTableWidgetItem(str(data[key]))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.loc_finder_table.setItem(row, col, item)
        self.loc_finder_table.resizeColumnsToContents()

    def plot_heatmap(self, data, labels):
        self.colormap_canvas.axes.clear()
        im = self.colormap_canvas.axes.imshow(data, cmap='summer')
        self.colormap_canvas.cbar = self.colormap_canvas.axes.figure.colorbar(im, ax=self.colormap_canvas.axes)
        self.colormap_canvas.cbar.ax.set_ylabel("", rotation=-90, va="bottom", fontsize=8)

        cursor = mplcursors.cursor(im, hover=True)
        @cursor.connect("add")
        def on_add(sel):
            sel.annotation.arrow_patch.set(arrowstyle="simple", fc="white", alpha=.5)
            sel.annotation.set_bbox(None)
            i, j = sel.target.index
            sel.annotation.set_text(labels[i][j])

        ax = self.colormap_canvas.axes
        ax.set_xticks(np.arange(len(data)))
        ax.set_yticks(np.arange(len(data)))
        ax.set_xticklabels(range(1, len(data) + 1))
        ax.set_yticklabels(range(1, len(data) + 1))
        ax.set_xlabel("Organism", fontsize=10)
        ax.set_ylabel("Organism", fontsize=10)
        ax.tick_params(axis='both', which='major', labelsize=8)

        def plot_cell_grid(data, ax=None, **kwargs):
            for x in range(data[0]):
                for y in range(data[1]):
                    rect = patches.Rectangle((x - .5, y - .5), 1, 1, fill=False, **kwargs)
                    ax.add_patch(rect)

        plot_cell_grid([len(data), len(data)], ax, color="black", linewidth=1)
        self.colormap_canvas.draw()

    def clear_shared_seeds_table(self):
        self.table2.setRowCount(0)

    def clear_loc_finder_table(self):
        self.loc_finder_table.setRowCount(0)

    def show_loading_window(self, value):
        self.loading_window.loading_bar.setValue(value)
        self.loading_window.show()

    def hide_loading_window(self):
        self.loading_window.hide()

    def update_loading_window(self, value, text):
        self.loading_window.loading_bar.setValue(value)
        self.loading_window.info_label.setText(text)

    def update_endo_dropdown(self, endos):
        self.endoBox.clear()
        self.endoBox.addItems(endos)

    def sort_table2(self, column):
        self.table2.sortItems(column)

    def sort_loc_finder_table(self, column):
        self.loc_finder_table.sortItems(column)

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=400, height=250, dpi=100):
        try:
            fig = Figure(dpi=dpi, tight_layout=True)
            self.axes = fig.add_subplot(111)
            self.axes.clear()
            super(MplCanvas, self).__init__(fig)
        except Exception as e:
            show_error("Error initializing MplCanvas class in population analysis.", e)

class LoadingWindow(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super(LoadingWindow, self).__init__()
        uic.loadUi(global_settings.get_ui_dir() + "/loading_data_form.ui", self)
        self.loading_bar.setValue(0)
        self.setWindowTitle("Loading Data")
        self.setWindowIcon(QtGui.QIcon(global_settings.get_assets_dir() + "cas9image.ico"))
        scale_ui(self, base_width=1920, base_height=1080, font_size=12, header_font_size=30, custom_scale_width=450, custom_scale_height=125)
