# Each widget used for matlab potting needs its own file to be promoted to it within
# the Qt editor. This file is for the bar graph of repeats vs chromosome. This just initializes
# the widget to be a matplotlib plot
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure


class repeats_vs_chromo(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(Figure())

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        #special alignment to make the graph fit better in the window given to it in the UI file
        self.canvas.figure.subplots_adjust(left=0.10, right=0.95, bottom=0.15, top=0.95)
        self.setLayout(vertical_layout)