from PyQt6 import QtWidgets, QtGui, QtCore
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon
from PyQt6 import uic
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import models.GlobalSettings as GlobalSettings
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from utils.ui import show_message, show_error, scale_ui, center_ui
from PyQt6.QtWidgets import QWidget, QComboBox, QVBoxLayout, QPushButton
import os

class MultitargetingWindowView(QWidget):
    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        self._init_ui()

    def _init_ui(self):
        try:
            self._load_ui_file()
            # self._init_ui_elements()
            # self._scale_ui()
        except Exception as e:
            self._handle_init_error(e)

    def _handle_init_error(self, e: Exception):
        error_msg = f"Error initializing MultitargetingWindowView: {str(e)}"
        self.settings.logger.error(error_msg, exc_info=True)
        show_error(self.settings, "Initialization Error", error_msg)
        raise

    def _load_ui_file(self):
        ui_file = os.path.join(self.settings.get_ui_dir(), "multitargeting_window.ui")
        uic.loadUi(ui_file, self)

    def scaleUI(self):
        # Implement scaling logic if needed
        pass

    def centerUI(self):
        # Implement centering logic if needed
        pass

    def update_table(self, data):
        # Implement table update logic
        pass

    def update_global_stats(self, model):
        # Implement global stats update logic
        pass

    def update_global_line_graph(self, data):
        # Implement line graph update logic
        pass

    def update_global_bar_graph(self, data):
        # Implement bar graph update logic
        pass

    # Add other methods as needed

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        try:
            fig = Figure(dpi=dpi, tight_layout=True)
            self.axes = fig.add_subplot(111)
            self.axes.clear()
            super(MplCanvas, self).__init__(fig)
        except Exception as e:
            show_error(self.settings, "Error initializing MplCanvas class in MultitargetingWindowView.", e)

class LoadingWindow(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(LoadingWindow, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "ui/loading_data_form.ui", self)
            self.setWindowIcon(QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.loading_bar.setValue(0)
            self.setWindowTitle("Loading Data")
            scale_ui(self, custom_scale_width=450, custom_scale_height=125)
        except Exception as e:
            show_error(self.settings, "Error initializing LoadingWindow class in MultitargetingWindowView.", e)

class MultitargetingStatisticsView(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        try:
            super(MultitargetingStatisticsView, self).__init__(parent)
            uic.loadUi(GlobalSettings.appdir + 'ui/multitargeting_stats.ui', self)
            self.setWindowIcon(QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.setWindowTitle("Statistics")
            scale_ui(self, custom_scale_width=275, custom_scale_height=185)
        except Exception as e:
            show_error(self.settings, "Error initializing MultitargetingStatisticsView class in MultitargetingWindowView.", e)

class SqlQuerySettingsView(QtWidgets.QMainWindow):
    def __init__(self):
        try:
            super(SqlQuerySettingsView, self).__init__()
            uic.loadUi(GlobalSettings.appdir + "ui/multitargeting_sql_settings.ui", self)
            self.setWindowTitle("SQL Settings")
            self.setWindowIcon(QIcon(GlobalSettings.appdir + "cas9image.ico"))
            self.row_count.setValidator(QtGui.QIntValidator())
            scale_ui(self, custom_scale_width=375, custom_scale_height=140)
        except Exception as e:
            show_error(self.settings, "Error initializing SqlQuerySettingsView class in MultitargetingWindowView.", e)
