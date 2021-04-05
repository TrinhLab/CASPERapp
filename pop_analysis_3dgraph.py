from PyQt5.QtWidgets import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from mpl_toolkits.mplot3d import Axes3D



class pop_analysis_3dgraph(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvasQTAgg(Figure())
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.canvas.axes = self.canvas.figure.add_subplot()

        #self.canvas.figure.subplots_adjust(left=0.1, right=0.9, bottom=0.05, top=0.8)

        self.canvas.figure.set_visible(False)
        self.setLayout(vertical_layout)
