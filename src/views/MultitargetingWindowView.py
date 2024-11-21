from typing import Optional
from PyQt6 import QtWidgets, uic, QtGui
from PyQt6.QtWidgets import QTableWidgetItem, QAbstractItemView
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from utils.ui import show_error

class MultitargetingWindowView(QtWidgets.QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        
        self.init_ui()

    def init_ui(self):
        try:
            uic.loadUi(self.settings.get_ui_dir_path() + '/multitargeting_window.ui', self)
            self._init_ui_components()
        except Exception as e:
            show_error(self.settings, "Error initializing MultitargetingWindowView", str(e))

    def _init_ui_components(self):
        self._init_grpSelectOrganism()
        self._init_grpSeedAnalysis()
        self._init_grpGlobalAnalysis()

        self.push_button_export_selected_gRNAs = self._find_widget('pbtnExportSelectedgRNAs', QtWidgets.QPushButton)

    def _init_grpSelectOrganism(self):
        self.combo_box_organism = self._find_widget('cmbOrganism', QtWidgets.QComboBox)
        self.combo_box_endonuclease = self._find_widget('cmbEndonuclease', QtWidgets.QComboBox)
        self.line_edit_max_results = self._find_widget('ledMaxResults', QtWidgets.QLineEdit)
        self.push_button_analyze = self._find_widget('pbtnAnalyze', QtWidgets.QPushButton)
        self.check_box_select_all = self._find_widget('chkSelectAll', QtWidgets.QCheckBox)
        self.table_seeds = self._find_widget('tblSeeds', QtWidgets.QTableWidget)

        # Set up table columns
        self.table_seeds.setColumnCount(8)
        self.table_seeds.setHorizontalHeaderLabels([
            "Seed", "Total Repeats", "Avg. Repeats/Scaffold", 
            "Consensus Sequence", "% Consensus", "Score", "PAM", "Strand"
        ])
        
        # Set table properties
        self.table_seeds.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_seeds.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_seeds.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        # Set minimum width for the table
        self.table_seeds.setMinimumWidth(650)

        # Add validation for max results line edit
        self.line_edit_max_results.setValidator(QtGui.QIntValidator())
        # Set default value
        self.line_edit_max_results.setText("1000")

    def _init_grpSeedAnalysis(self):
        # Get the tab widget
        self.tab_widget_seed_analysis = self._find_widget('tabsSeedAnalysis', QtWidgets.QTabWidget)

        # Get screen height for scaling
        screen = self.screen()
        height = screen.geometry().height()
        # Set tab widget to take up less vertical space (e.g., 30% of screen height)
        tab_height = int((height * 30) / 100)
        self.tab_widget_seed_analysis.setMinimumHeight(tab_height)
        self.tab_widget_seed_analysis.setMaximumHeight(tab_height)

        # Initialize other widgets
        self.tab_chromosome_viewer = self._find_widget('tabChromosomeViewer', QtWidgets.QWidget)
        self.graphical_view_chromosome = self._find_widget('graphviewChromosome', QtWidgets.QGraphicsView)
        self.scroll_chromosome = self._find_widget('scrollChromosome', QtWidgets.QScrollArea)
        
        # Create a widget to hold the chromosome visualizations
        self.chromosome_content_widget = QtWidgets.QWidget()
        self.chromosome_layout = QtWidgets.QVBoxLayout(self.chromosome_content_widget)
        self.scroll_chromosome.setWidget(self.chromosome_content_widget)
        self.scroll_chromosome.setWidgetResizable(True)

        self.tab_seed_distribution = self._find_widget('tabSeedDistribution', QtWidgets.QWidget)
        self.plot_repeat_vs_chromosome = self._find_widget('plotRepeatVsChromosome', QtWidgets.QWidget)

        # Initialize scene for chromosome details
        self.scene = QtWidgets.QGraphicsScene()
        self.scene2 = QtWidgets.QGraphicsScene()
        self.graphical_view_chromosome.setScene(self.scene2)
        
        # Set up event filters
        self.scroll_chromosome.viewport().installEventFilter(self)
        self.graphical_view_chromosome.viewport().installEventFilter(self)

        self.canvas_chromosome_map = {}

    def _init_grpGlobalAnalysis(self):
        self.tab_statistics_overview = self._find_widget('tabStatisticsOverview', QtWidgets.QWidget)

        self.tab_repeats_vs_seed = self._find_widget('tabRepeatsVsSeed', QtWidgets.QWidget)
        self.plot_repeats_vs_seed = self._find_widget('plotRepeatsVsSeed', QtWidgets.QWidget)

        self.tab_sequences_vs_repeats = self._find_widget('tabSequencesVsRepeats', QtWidgets.QWidget)
        self.plot_sequences_vs_repeats = self._find_widget('plotSequencesVsRepeats', QtWidgets.QWidget)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget

    def update_seeds_table(self, data):
        self.table_seeds.setRowCount(len(data))
        for row, row_data in enumerate(data):
            # Unpack the data
            (seed, total_count, avg_per_scaffold, sequences, 
             consensus_percent, score, pam, strand) = row_data

            # Set data in table cells
            self.table_seeds.setItem(row, 0, QTableWidgetItem(str(seed)))
            self.table_seeds.setItem(row, 1, QTableWidgetItem(str(total_count)))
            self.table_seeds.setItem(row, 2, QTableWidgetItem(f"{avg_per_scaffold:.2f}"))
            self.table_seeds.setItem(row, 3, QTableWidgetItem(sequences.split(',')[0]))  # First sequence
            self.table_seeds.setItem(row, 4, QTableWidgetItem(f"{consensus_percent:.1f}"))
            self.table_seeds.setItem(row, 5, QTableWidgetItem(str(score)))
            self.table_seeds.setItem(row, 6, QTableWidgetItem(str(pam)))
            self.table_seeds.setItem(row, 7, QTableWidgetItem(str(strand)))

            # Set alignment for all cells
            for col in range(8):
                item = self.table_seeds.item(row, col)
                if item:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.table_seeds.resizeColumnsToContents()

    def setup_plots(self):
        """Initialize the matplotlib plots"""
        self.repeats_vs_seed_canvas = MplCanvas(self, width=8, height=6)
        self.sequences_vs_repeats_canvas = MplCanvas(self, width=8, height=6)
        self.repeat_vs_chromosome_canvas = MplCanvas(self, width=8, height=6)

        # Add canvases to their respective layouts without toolbars
        for plot_widget, canvas in [
            (self.plot_repeats_vs_seed, self.repeats_vs_seed_canvas),
            (self.plot_sequences_vs_repeats, self.sequences_vs_repeats_canvas),
            (self.plot_repeat_vs_chromosome, self.repeat_vs_chromosome_canvas)
        ]:
            layout = QtWidgets.QVBoxLayout(plot_widget)
            layout.setContentsMargins(0, 0, 0, 0)  # Reduce margins
            layout.addWidget(canvas)

    def update_plots(self, repeats_data, sequences_data, chromosome_data):
        """Update all plots with new data"""
        self._update_repeats_vs_seed_plot(repeats_data)
        self._update_sequences_vs_repeats_plot(sequences_data)
        self._update_repeat_vs_chromosome_plot(chromosome_data)

    def _update_repeats_vs_seed_plot(self, data):
        """Update the repeats vs seed line plot"""
        try:
            self.logger.debug("Starting repeats vs seed plot update")
            print(f"data: {data}")
            
            self.repeats_vs_seed_canvas.axes.clear()
            
            if data and 'counts' in data:
                y1 = data['counts']
                x = range(len(y1))
                
                self.repeats_vs_seed_canvas.axes.plot(x, y1, linewidth=1.5, marker='.', markersize=3)
                self.repeats_vs_seed_canvas.axes.set_xlabel('Seed ID Number', fontsize=10)
                self.repeats_vs_seed_canvas.axes.set_ylabel('Number of Repeats', fontsize=10)
                self.repeats_vs_seed_canvas.axes.set_title('Number of Repeats per Seed ID Number', fontsize=10)
                self.repeats_vs_seed_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
                self.repeats_vs_seed_canvas.axes.grid(True, linestyle='--', alpha=0.7)
                
                # Store statistics if needed
                if 'stats' in data:
                    self.average = data['stats']['average']
                    self.mode = data['stats']['mode']
                    self.median = data['stats']['median']
                    self.repeat_count = data['stats']['repeat_count']
                
                # Force draw
                self.repeats_vs_seed_canvas.draw()
                
            else:
                self.logger.warning("No valid data for plotting")
                self.repeats_vs_seed_canvas.draw()  # Still need to draw even when clearing
                
        except Exception as e:
            self.logger.error(f"Error updating repeats vs seed plot: {str(e)}")

    def _update_sequences_vs_repeats_plot(self, data):
        """Update the sequences vs repeats plot"""
        try:
            self.sequences_vs_repeats_canvas.axes.clear()
            
            if data and 'x_vals' in data and 'y_vals' in data:
                x = data['x_vals']  # Number of repeats
                y = data['y_vals']  # Number of sequences
                
                self.sequences_vs_repeats_canvas.axes.scatter(x, y, s=10)
                self.sequences_vs_repeats_canvas.axes.set_yscale('log')
                self.sequences_vs_repeats_canvas.axes.set_xlabel('Number of Repeats', fontsize=10)
                self.sequences_vs_repeats_canvas.axes.set_ylabel('Number of Sequences', fontsize=10)
                self.sequences_vs_repeats_canvas.axes.set_title('Number of Sequences per Number of Repeats', fontsize=10)
                self.sequences_vs_repeats_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
                self.sequences_vs_repeats_canvas.axes.grid(True, linestyle='--', alpha=0.7)
                
                if x:  
                    self.sequences_vs_repeats_canvas.axes.set_xlim(x[0] - 0.5, x[-1] + 0.5)
                
            self.sequences_vs_repeats_canvas.draw()
            
        except Exception as e:
            self.logger.error(f"Error updating sequences vs repeats plot: {str(e)}")

    def _update_repeat_vs_chromosome_plot(self, data):
        """Update the chromosome bar plot"""
        try:
            self.repeat_vs_chromosome_canvas.axes.clear()
            
            if data:
                y = []
                x_labels = []
                
                # Get sorted chromosome numbers and their counts
                for chromo in sorted(data.keys()):
                    x_labels.append(chromo)
                    y.append(data[chromo])
                    
                x = list(range(0, len(x_labels)))

                # Create the bar plot
                self.repeat_vs_chromosome_canvas.axes.bar(x, y, align='center')
                
                # Set integer y-axis
                self.repeat_vs_chromosome_canvas.axes.yaxis.set_major_locator(MaxNLocator(integer=True))
                
                # Set y-axis limits
                self.repeat_vs_chromosome_canvas.axes.set_ylim(0, max(y) + 1)
                
                # Set x-axis ticks and labels
                self.repeat_vs_chromosome_canvas.axes.set_xticks(x)
                self.repeat_vs_chromosome_canvas.axes.set_xticklabels(x_labels)
                
                # If many chromosomes, show only some labels
                if len(x_labels) > 10:
                    tick_spacing = round(len(x_labels)/10)
                    for i, t in enumerate(self.repeat_vs_chromosome_canvas.axes.get_xticklabels()):
                        if (i % tick_spacing) != 0:
                            t.set_visible(False)
                
                # Set labels and title
                self.repeat_vs_chromosome_canvas.axes.set_xlabel('Chromosome', fontsize=10)
                self.repeat_vs_chromosome_canvas.axes.set_ylabel('Number of Repeats', fontsize=10)
                self.repeat_vs_chromosome_canvas.axes.set_title('Repeats per Chromosome', fontsize=10)
                
                # Set tick label size
                self.repeat_vs_chromosome_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
                
            self.repeat_vs_chromosome_canvas.draw()
            
        except Exception as e:
            self.logger.error(f"Error updating repeat vs chromosome plot: {str(e)}")

    def fill_chromosome_viewer(self, seed_data, event_data):
        try:
            # Clear out old widgets in layout
            for i in reversed(range(self.chromosome_layout.count())):
                self.chromosome_layout.itemAt(i).widget().setParent(None)

            # Get sorted list of chromosomes
            chromo_keys = sorted(list(seed_data.keys()))

            # Get screen height for scaling
            screen = self.screen()
            height = screen.geometry().height()
            groupbox_height = int((height * 100) / 1080)

            # Create visualization for each chromosome
            for chromo in chromo_keys:
                group_box = QtWidgets.QGroupBox()
                group_box.setTitle(f"Chromosome {chromo}")
                group_box.setMinimumHeight(groupbox_height)
                group_box.setMaximumHeight(groupbox_height)
                layout = QtWidgets.QVBoxLayout(group_box)

                # Create canvas for this chromosome
                canvas = MplCanvas()
                canvas.axes.eventplot(seed_data[chromo])
                canvas.mpl_connect("motion_notify_event", self._chromosome_event_handler)

                # Add border lines
                canvas.axes.hlines(1.5, -0.01, 1.01, colors="Black", linewidth=1.5)
                canvas.axes.hlines(0.5, -0.01, 1.01, colors="Black", linewidth=1.5)
                canvas.axes.vlines(-0.01, 0.5, 1.5, colors="Black", linewidth=1.5)
                canvas.axes.vlines(1.01, 0.5, 1.5, colors="Black", linewidth=1.5)
                
                # Set axis limits
                canvas.axes.set_ylim(0.45, 1.55)
                canvas.axes.set_xlim(-0.05, 1.05)
                canvas.axes.axis('off')
                canvas.draw()

                # Store canvas mapping
                self.canvas_chromosome_map[canvas] = chromo
                self.event_data = event_data  # Store event data for hover details

                layout.addWidget(canvas)
                self.chromosome_layout.addWidget(group_box)

        except Exception as e:
            show_error(self.settings, "Error in fill_chromosome_viewer", str(e))

    def _chromosome_event_handler(self, event):
        """Handle mouse events on chromosome visualization"""
        try:
            # Get current mouse location
            x = event.xdata
            y = event.y
            if x is None:  # Mouse outside the plot
                return

            # Get event data relative to the canvas
            curr_chromosome = self.canvas_chromosome_map[event.canvas]
            chromosome_seed_data = self.event_data[curr_chromosome]

            # Get targets within range of mouse location
            local_targets = []
            for entry in chromosome_seed_data:
                try:
                    if abs(x - entry[0]) <= 0.001:
                        local_targets.append(entry)
                except:
                    pass

            # Update viewer with target details if found
            if local_targets:
                self.scene2 = QtWidgets.QGraphicsScene()
                self.graphical_view_chromosome.setScene(self.scene2)
                
                output = ""
                for target in local_targets:
                    output += f"Location: {target[1]} | Seq: {target[2]} | PAM: {target[3]} | SCR: {target[4]} | DIRA: {target[5]}\n"

                text = self.scene2.addText(output)
                font = QtGui.QFont()
                font.setPointSize(self.settings.fontSize if hasattr(self.settings, 'fontSize') else 12)
                text.setFont(font)

        except Exception as e:
            self.logger.error(f"Error in chromosome event handler: {str(e)}")

    def update_statistics_labels(self, total_repeats, avg_repeats, median_repeats, mode_repeats):
        """Update the statistics overview labels with new values"""
        try:
            # Find and update the statistics labels
            total_label = self._find_widget('lblTotalRepeatsValue', QtWidgets.QLabel)
            avg_label = self._find_widget('lblAverageRepeatsValue', QtWidgets.QLabel)
            median_label = self._find_widget('lblMedianRepeatsValue', QtWidgets.QLabel)
            mode_label = self._find_widget('lblModeRepeatsValue', QtWidgets.QLabel)
            
            if total_label:
                total_label.setText(str(round(float(total_repeats), 1)))
            if avg_label:
                avg_label.setText(str(round(float(avg_repeats), 1)))
            if median_label:
                median_label.setText(str(round(float(median_repeats), 1)))
            if mode_label:
                mode_label.setText(str(round(float(mode_repeats), 1)))
            
        except Exception as e:
            self.logger.error(f"Error updating statistics labels: {str(e)}")

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
