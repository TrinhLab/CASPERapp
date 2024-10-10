from PyQt6.QtWidgets import QMainWindow
from views.MultitargetingWindowView import MultitargetingWindowView
from models.MultitargetingWindowModel import MultitargetingWindowModel

class MultitargetingWindowController(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self._view = MultitargetingWindowView(global_settings)
        self._model = MultitargetingWindowModel(global_settings)
        self.setCentralWidget(self._view)
        self.setup_connections()

    @property
    def view(self):
        return self._view

    @property
    def model(self):
        return self._model

    def show(self):
        super().show()

    def setup_connections(self):
        # Connect signals from view to controller methods
        self.view.comboBox_organism.currentIndexChanged.connect(self.update_endo_list)
        self.view.comboBox_endo.currentIndexChanged.connect(self.load_data)
        self.view.pushButton_load.clicked.connect(self.load_data)
        # Add more connections as needed

    def update_endo_list(self):
        organism = self.view.comboBox_organism.currentText()
        endos = self.model.get_endos_for_organism(organism)
        self.view.comboBox_endo.clear()
        self.view.comboBox_endo.addItems(endos)

    def load_data(self):
        organism = self.view.comboBox_organism.currentText()
        endo = self.view.comboBox_endo.currentText()
        self.model.set_files(organism, endo)
        
        # Load and process data
        repeats_data = self.model.get_repeats_data()
        kstats = self.model.get_kstats()
        
        # Update view with data
        self.view.update_table(repeats_data)
        self.view.update_global_stats(self.model)
        
        # Load and update graphs
        self.load_graphs()

    def load_graphs(self):
        # Load data for graphs
        seeds_vs_repeats = self.model.get_seeds_vs_repeats_data()
        repeats_vs_seeds = self.model.get_repeats_vs_seeds_data()
        
        # Update graphs in view
        self.view.update_global_line_graph(seeds_vs_repeats)
        self.view.update_global_bar_graph(repeats_vs_seeds)

    def initialize(self):
        organisms = self.model.get_organisms()
        self.view.comboBox_organism.addItems(organisms)
        self.view.scaleUI()
        self.view.centerUI()
