# Each widget used for matlab potting needs its own file to be promoted to it within
# the Qt editor. This file is for the line graph of repeats vs. seeds. This just initializes
# the widget to be a matplotlib plot
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib_venn import venn2, venn3
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class pop_analysis_venn_diagram(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        set1 = set(['A', 'B', 'C', 'D'])
        set2 = set(['B', 'C', 'D', 'E'])
        set3 = set(['C', 'D', ' E', 'F', 'G'])


        self.canvas = FigureCanvas(plt.figure(figsize=(7.5,7.5), dpi=100))
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('A', 'B', 'C'))
        self.setLayout(vertical_layout)
        #self.canvas.figure.subplots_adjust(left=0.10, right=0.95, bottom=0.15, top=0.95)
        self.canvas.figure.tight_layout()

        #self.canvas.show()
        #self.canvas.axes = venn3([set1, set2, set3], ('set1', 'set2', 'set3'))


