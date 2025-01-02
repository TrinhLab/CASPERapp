from PyQt6 import QtWidgets, uic, QtGui, QtCore
from PyQt6.QtWidgets import QHeaderView, QAbstractItemView
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import mplcursors
import numpy as np
import matplotlib.patches as patches
from utils.ui import show_error
import copy

class PopulationAnalysisWindowView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        
        # Connect to theme change signal
        self.settings.theme_changed.connect(self._on_theme_changed)
        
        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/population_analysis.ui', self)
            self._init_ui_components()
            self._set_styles()
        except Exception as e:
            show_error(self.settings, "Error initializing PopulationAnalysisWindowView", str(e))

    def _init_ui_components(self):
        self._init_grpSelectOrganisms()
        self._init_grpSeedAnalysis()

        self.push_button_export_selected_gRNAs = self._find_widget('pbtnExportSelectedgRNAs', QtWidgets.QPushButton)

    def _init_grpSelectOrganisms(self):
        try:
            self.combo_box_endonuclease = self._find_widget('cmbEndonuclease', QtWidgets.QComboBox)
            self.table_organism = self._find_widget('tblOrganism', QtWidgets.QTableWidget)
            self.push_button_analyze_organism = self._find_widget('pbtnAnalyzeOrganism', QtWidgets.QPushButton)

            # Find the tab widget and heatmap widget
            self.tab_widget_shared_seeds_heatmap = self._find_widget('tabsSharedSeedHeatmap', QtWidgets.QTabWidget)
            self.tab_shared_seed_heatmap = self._find_widget('tabSharedSeedHeatmap', QtWidgets.QWidget)
            self.heatmap_seed = self._find_widget('heatmapSeed', QtWidgets.QWidget)
            
            # Create layout for heatmap
            self.colormap_layout = QtWidgets.QVBoxLayout(self.heatmap_seed)
            self.colormap_layout.setContentsMargins(0, 0, 0, 0)
            
            # Create the matplotlib canvas
            self.colormap_canvas = MplCanvas(self)
            self.colormap_layout.addWidget(self.colormap_canvas)

            # Set up the organism table
            self.table_organism.setColumnCount(1)
            self.table_organism.setShowGrid(False)
            self.table_organism.setHorizontalHeaderLabels(["Organism"])
            self.table_organism.horizontalHeader().setSectionsClickable(True)
            self.table_organism.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.table_organism.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
            self.table_organism.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            self.table_organism.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
            
        except Exception as e:
            self.logger.error(f"Error in _init_grpSelectOrganisms: {str(e)}")
            self.logger.exception("Full traceback:")
            show_error(self.settings, "Error initializing select organisms group", str(e))

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
            for col, value in enumerate(data):
                item = QtWidgets.QTableWidgetItem(str(value))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table_seed.setItem(row, col, item)
        self.table_seed.resizeColumnsToContents()

    def update_loc_finder_table(self, loc_data):
        self.table_locations.setRowCount(len(loc_data))
        for row, data in enumerate(loc_data):
            for col, key in enumerate(['seed', 'sequence', 'organism', 'chromosome', 'location']):
                item = QtWidgets.QTableWidgetItem(str(data[key]))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table_locations.setItem(row, col, item)
        self.table_locations.resizeColumnsToContents()

    def plot_heatmap(self, data, labels):
        """Plot heatmap of shared seeds between organisms"""
        try:
            self.logger.debug("Starting plot_heatmap")
            self.logger.debug(f"Data shape: {np.array(data).shape}")
            self.logger.debug(f"Data: {data}")
            self.logger.debug(f"Labels: {labels}")
            
            # Clear the previous plot and colorbar safely
            self.colormap_canvas.axes.clear()
            if hasattr(self.colormap_canvas, 'cbar'):
                try:
                    self.colormap_canvas.cbar.remove()
                except:
                    self.logger.debug("Could not remove old colorbar, creating new figure")
                    # If colorbar removal fails, create new figure and canvas
                    self.colormap_canvas.fig.clear()
                    self.colormap_canvas.axes = self.colormap_canvas.fig.add_subplot(111)
            
            # Create a copy of data for labels
            labels_data = copy.deepcopy(data)
            
            # Set diagonal elements to 0 for visualization
            for i in range(len(data)):
                data[i][i] = 0
            
            # Create the heatmap
            self.logger.debug("Creating heatmap")
            im = self.colormap_canvas.axes.imshow(data, cmap='summer')
            
            # Add colorbar
            self.logger.debug("Adding colorbar")
            self.colormap_canvas.cbar = self.colormap_canvas.fig.colorbar(im, ax=self.colormap_canvas.axes)
            
            # Determine text color based on theme
            text_color = 'white' if self.settings.get_theme() == 'dark' else 'black'
            
            self.colormap_canvas.cbar.ax.set_ylabel("", rotation=-90, va="bottom", fontsize=8, color=text_color)
            self.colormap_canvas.cbar.ax.tick_params(colors=text_color)

            # Add hover annotations
            self.logger.debug("Setting up hover annotations")
            cursor = mplcursors.cursor(im, hover=True)
            @cursor.connect("add")
            def on_add(sel):
                sel.annotation.arrow_patch.set(arrowstyle="simple", fc="white", alpha=.5)
                sel.annotation.set_bbox(None)
                i, j = sel.target.index
                # Show the actual number of shared seeds
                sel.annotation.set_text(str(labels_data[i][j]))

            # Set up axes
            self.logger.debug("Setting up axes")
            ax = self.colormap_canvas.axes
            ax.set_xticks(np.arange(len(data)))
            ax.set_yticks(np.arange(len(data)))
            
            # Use numbers for both x-axis and y-axis
            x_labels = [str(i+1) for i in range(len(data))]
            y_labels = [str(i+1) for i in range(len(data))]  # Removed "Organism" prefix
            ax.set_xticklabels(x_labels, color=text_color)
            ax.set_yticklabels(y_labels, color=text_color)
            
            # Rotate labels
            self.logger.debug("Rotating labels")
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
            
            # Add grid
            self.logger.debug("Adding grid")
            for i in range(len(data)):
                for j in range(len(data)):
                    ax.add_patch(patches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, 
                        fill=False, color="black", linewidth=1
                    ))
            
            ax.set_xlabel("Organism", fontsize=10, color=text_color)
            ax.set_ylabel("Organism", fontsize=10, color=text_color)
            ax.tick_params(axis='both', which='major', labelsize=8, colors=text_color)
            
            # Adjust layout and draw
            self.logger.debug("Adjusting layout")
            self.colormap_canvas.fig.tight_layout()
            
            self.logger.debug("Drawing canvas")
            self.colormap_canvas.draw()
            
            self.logger.debug("Completed plot_heatmap")
            
        except Exception as e:
            self.logger.error(f"Error plotting heatmap: {str(e)}")
            self.logger.exception("Full traceback:")
            show_error(self.settings, "Error plotting heatmap", str(e))

    def clear_shared_seeds_table(self):
        self.table_seed.setRowCount(0)

    def clear_loc_finder_table(self):
        self.table_locations.setRowCount(0)

    def update_endo_dropdown(self, endos):
        """Update the endonuclease dropdown with the provided options"""
        try:
            self.combo_box_endonuclease.clear()
            self.combo_box_endonuclease.addItems(endos)
        except Exception as e:
            self.logger.error(f"Error updating endonuclease dropdown: {str(e)}")
            show_error(self.settings, "Error updating endonuclease dropdown", str(e))

    def sort_table2(self, column):
        self.table_seed.sortItems(column)

    def sort_loc_finder_table(self, column):
        self.table_locations.sortItems(column)

    def _on_theme_changed(self, theme):
        """Handle theme changes by updating the plot"""
        try:
            if hasattr(self, 'colormap_canvas') and hasattr(self.colormap_canvas, 'axes'):
                text_color = 'white' if theme == 'dark' else 'black'
                
                # Update axis labels
                self.colormap_canvas.axes.xaxis.label.set_color(text_color)
                self.colormap_canvas.axes.yaxis.label.set_color(text_color)
                
                # Update tick labels
                self.colormap_canvas.axes.tick_params(colors=text_color)
                
                # Update colorbar if it exists
                if hasattr(self.colormap_canvas, 'cbar'):
                    self.colormap_canvas.cbar.ax.set_ylabel("", rotation=-90, va="bottom", fontsize=8, color=text_color)
                    self.colormap_canvas.cbar.ax.tick_params(colors=text_color)
                
                # Redraw the canvas
                self.colormap_canvas.draw()
                
        except Exception as e:
            self.logger.error(f"Error updating plot theme: {str(e)}")

    def _set_styles(self):
        """Apply the global groupbox style"""
        try:
            style = self.settings.get_groupbox_style()
            for groupbox in self.findChildren(QtWidgets.QGroupBox):
                groupbox.setStyleSheet(style)
        except Exception as e:
            self.logger.error(f"Error setting styles: {str(e)}")

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super().__init__(self.fig)
        
        # Set background colors based on current theme
        self.update_colors(parent.settings.get_theme() if parent else 'light')
        
        # Enable tight layout
        self.fig.tight_layout()
        
    def update_colors(self, theme):
        """Update figure and axes colors based on theme"""
        if theme == 'dark':
            self.fig.patch.set_facecolor('none')
            self.axes.set_facecolor('#2d2d2d')  # Dark background for plot area
        else:
            self.fig.patch.set_facecolor('none')
            self.axes.set_facecolor('white')
