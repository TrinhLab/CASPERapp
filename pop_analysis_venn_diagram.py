# Each widget used for matlab potting needs its own file to be promoted to it within
# the Qt editor. This file is for the line graph of repeats vs. seeds. This just initializes
# the widget to be a matplotlib plot
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib_venn import venn3_unweighted
import matplotlib.pyplot as plt

class pop_analysis_venn_diagram(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(plt.figure(figsize=(7.5,7.5)))

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)

        self.setLayout(vertical_layout)
        self.canvas.figure.tight_layout()



