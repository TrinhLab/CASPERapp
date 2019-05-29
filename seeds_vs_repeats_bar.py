# Each widget used for matlab potting needs its own file to be promoted to it within
# the Qt editor. This file is for the bar graph of repeats vs. seeds. This just initializes
# the widget to be a matplotlib plot
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure


class seeds_vs_repeats_bar(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(Figure())

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)

        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)