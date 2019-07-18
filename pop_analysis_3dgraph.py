from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvas
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.figure import Figure

class pop_analysis_3dgraph(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(Figure(figsize=(2, 2)))
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.canvas.axes = self.canvas.figure.add_subplot(111, projection='3d')
        self.setLayout(vertical_layout)
        self.canvas.figure.tight_layout()

        for c, z in zip(['r', 'g', 'b', 'y'], [30, 20, 10, 0]):
            xs = np.arange(20)
            ys = np.random.rand(20)

            # You can provide either a single color or an array. To demonstrate this,
            # the first bar of each set will be colored cyan.
            cs = [c] * len(xs)
            cs[0] = 'c'
            self.canvas.axes.bar(xs, ys, zs=z, zdir='y', color=cs, alpha=0.8)

        self.canvas.axes.set_xlabel('X')
        self.canvas.axes.set_ylabel('Y')
        self.canvas.axes.set_zlabel('Z')