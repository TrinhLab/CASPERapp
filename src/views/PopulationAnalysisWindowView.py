from PyQt6 import QtWidgets, uic, QtGui, QtCore
from PyQt6.QtWidgets import QHeaderView, QAbstractItemView
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import mplcursors
import numpy as np
import matplotlib.patches as patches
from utils.ui import show_error

class PopulationAnalysisWindowView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/population_analysis.ui', self)
            self._init_ui_components()
        except Exception as e:
            show_error(self.settings, "Error initializing PopulationAnalysisWindowView", str(e))

    def _init_ui_components(self):
        self._init_grpSelectOrganisms()
        self._init_grpSeedAnalysis()
        # self._init_colormap()

    def _init_grpSelectOrganisms(self):
        self.combo_box_endonuclease = self._find_widget('cmbEndonuclease', QtWidgets.QComboBox)
        print(self.combo_box_endonuclease)
        self.table_organism = self._find_widget('tblOrganism', QtWidgets.QTableWidget)
        self.push_button_analyze_organism = self._find_widget('pbtnAnalyzeOrganism', QtWidgets.QPushButton)

        self.tab_widget_shared_seeds_heatmap = self._find_widget('tabsSharedSeedHeatmap', QtWidgets.QTabWidget)
        self.tab_shared_seed_heatmap = self._find_widget('tabSharedSeedHeatmap', QtWidgets.QWidget)
        self.heatmap_seed = self._find_widget('heatmapSeed', QtWidgets.QWidget)

        self.table_organism.setColumnCount(1)
        self.table_organism.setShowGrid(False)
        self.table_organism.setHorizontalHeaderLabels(["Organism"])
        self.table_organism.horizontalHeader().setSectionsClickable(True)
        self.table_organism.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table_organism.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_organism.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_organism.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)

    def _init_grpSeedAnalysis(self):
        self.line_edit_seed = self._find_widget('ledSeed', QtWidgets.QLineEdit)
        self.push_button_query_seed = self._find_widget('pbtnQuerySeed', QtWidgets.QPushButton)
        self.push_button_clear_seeds = self._find_widget('pbtnClearSeeds', QtWidgets.QPushButton)
        self.table_seed = self._find_widget('tblSeed', QtWidgets.QTableWidget)

        self.table_seed.setColumnCount(9)
        self.table_seed.setShowGrid(False)
        self.table_seed.setHorizontalHeaderLabels([
            "Seed", "% Coverage", "Total Repeats", "Avg. Repeats/Scaffold",
            "Consensus Sequence", "% Consensus", "Score", "PAM", "Strand"
        ])
        self.table_seed.horizontalHeader().setSectionsClickable(True)
        self.table_seed.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_seed.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_seed.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)

        self.push_button_find_locations = self._find_widget('pbtnFindLocations', QtWidgets.QPushButton)
        self.push_button_clear_locations = self._find_widget('pbtnClearLocations', QtWidgets.QPushButton)

        self.table_locations = self._find_widget('tblLocation', QtWidgets.QTableWidget)

        self.table_locations.setColumnCount(5)
        self.table_locations.setShowGrid(False)
        self.table_locations.setHorizontalHeaderLabels([
                "Seed ID", "Sequence", "Organism", "Scaffold", "Location"
        ])
        self.table_locations.horizontalHeader().setSectionsClickable(True)
        self.table_locations.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table_locations.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

    # def _init_colormap(self):
    #     self.colormap_figure = self._find_widget('wgtColormap', QtWidgets.QWidget)
    #     if self.colormap_figure:
    #         self.colormap_layout = QtWidgets.QVBoxLayout()
    #         self.colormap_layout.setContentsMargins(0, 0, 0, 0)
    #         self.colormap_canvas = MplCanvas(self)
    #         self.colormap_layout.addWidget(self.colormap_canvas)
    #         self.colormap_figure.setLayout(self.colormap_layout)

    def _find_widget(self, name: str, widget_type: type) -> QtWidgets.QWidget:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def get_selected_endo(self):
        return self.combo_box_endonuclease.currentText()

    def get_selected_organisms(self):
        return [index.row() for index in self.org_Table.selectionModel().selectedRows()]

    def get_selected_seeds(self):
        return [self.table_seed.item(index.row(), 0).text() for index in self.table_seed.selectionModel().selectedRows()]

    def get_seed_input(self):
        return self.seed_input.text()

    def get_selected_seeds_for_export(self):
        return [item.text() for item in self.table_seed.selectedItems() if item.column() == 0]

    def update_org_table(self, org_data):
        self.table_organism.setRowCount(len(org_data))
        for row, (org_name, cspr_file, db_file) in enumerate(org_data):
            item = QtWidgets.QTableWidgetItem(org_name)
            item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.table_organism.setItem(row, 0, item)
        self.table_organism.resizeColumnsToContents()

    def update_shared_seeds_table(self, seed_data):
        self.table_seed.setRowCount(len(seed_data))
        for row, data in enumerate(seed_data):
            print(data)
            for col, value in enumerate(data):
                item = QtWidgets.QTableWidgetItem(str(value))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table_seed.setItem(row, col, item)
        self.table_seed.resizeColumnsToContents()

    def update_loc_finder_table(self, loc_data):
        self.table_locations.setRowCount(len(loc_data))
        for row, data in enumerate(loc_data):
            print(data)
            for col, key in enumerate(['seed', 'sequence', 'organism', 'chromosome', 'location']):
                item = QtWidgets.QTableWidgetItem(str(data[key]))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table_locations.setItem(row, col, item)
        self.table_locations.resizeColumnsToContents()

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
        self.table_seed.setRowCount(0)

    def clear_loc_finder_table(self):
        self.loc_finder_table.setRowCount(0)

    def update_endo_dropdown(self, endos):
        """Update the endonuclease dropdown with the provided options"""
        try:
            self.logger.info("Starting update_endo_dropdown")
            self.logger.debug(f"Received endos: {endos}")

            print(self.combo_box_endonuclease)
            
            # if not self.combo_box_endonuclease:
                # self.logger.error("combo_box_endonuclease is None")
                # return
                
            self.combo_box_endonuclease.clear()
            self.combo_box_endonuclease.addItems(endos)
            
            self.logger.info(f"Updated endonuclease dropdown with {len(endos)} options")
            self.logger.debug(f"Current items in dropdown: {[self.combo_box_endonuclease.itemText(i) for i in range(self.combo_box_endonuclease.count())]}")
        except Exception as e:
            self.logger.error(f"Error updating endonuclease dropdown: {str(e)}")
            self.logger.exception("Full traceback:")
            show_error(self.settings, "Error updating endonuclease dropdown", str(e))

    def sort_table2(self, column):
        self.table_seed.sortItems(column)

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
