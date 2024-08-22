import models.GlobalSettings as GlobalSettings
from utils.ui import show_error, scale_ui, center_ui
from PyQt5 import QtWidgets, QtCore, uic, Qt
from views.annotation_functions import *
import traceback

logger = GlobalSettings.logger

class AnnotationWindow(QtWidgets.QMainWindow):
    def __init__(self, info_path):
        super(AnnotationWindow, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'ui/annotation_details.ui', self)
        self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.ico"))
        self.Submit_button.clicked.connect(self.submit)
        self.Go_Back_Button.clicked.connect(self.go_Back)
        self.select_all_checkbox.stateChanged.connect(self.select_all_genes)
        self.fontSize = 12
        self.mainWindow = ""
        self.type = ""
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.switcher_table = [1, 1, 1, 1, 1, 1, 1, 1]
        self.tableWidget.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.tableWidget.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.tableWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.tableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tableWidget.setAutoScroll(False)
        self.tableWidget.horizontalHeader().sectionClicked.connect(self.table_sorting)

        scale_ui(self, custom_scale_width=1150, custom_scale_height=650)

    def table_sorting(self, logicalIndex):
        try:
            self.switcher_table[logicalIndex] *= -1
            order = QtCore.Qt.DescendingOrder if self.switcher_table[logicalIndex] == -1 else QtCore.Qt.AscendingOrder
            self.tableWidget.sortItems(logicalIndex, order)
        except Exception as e:
            show_error("table_sorting() in Annotation Window", e)

      #submit selected rows for results to process
    def submit(self):
        try:
            self.mainWindow.collect_table_data_nonkegg()
            self.mainWindow.show()  # Open main window back up
        except Exception as e:
            logger.critical("Error in submit() in AnnotationWindow.")
            logger.critical(e)
            logger.critical(traceback.format_exc())
            show_error("submit() in AnnotationWindow", e)
        finally:
            self.hide()  # Close annotation window regardless of success or failure

    #go back to main
    def go_Back(self):
        try:
            self.tableWidget.clear()
            self.mainWindow.searches.clear()
            self.tableWidget.setColumnCount(0)
            self.mainWindow.show()
            self.mainWindow.progressBar.setValue(0)
            self.hide()
        except Exception as e:
            show_error("go_Back() in AnnotationWindow", e)
            self.mainWindow.checkBoxes.clear()
    
    # this is the connection for the select all checkbox - selects/deselects all the genes in the table
    # this function is very similar to the other fill_table, it just works with the other types of annotation files
    def fill_table_nonKegg(self, mainWindow, results_list):
        try:
            self.tableWidget.clearContents()
            self.mainWindow = mainWindow
            self.tableWidget.setColumnCount(5)
            self.mainWindow.progressBar.setValue(85)
            self.tableWidget.setHorizontalHeaderLabels(["Feature Type", "Chromosome/Scaffold #", "Feature ID/Locus Tag", "Feature Name", "Feature Description"])

            mainWindow.checkBoxes = []
            self.type = "nonkegg"
           
            for index, result in enumerate(results_list[:1000]): # Limit to first 1000 results
                self.tableWidget.setRowCount(index + 1) # Increment table row count
                chrom, feature = result

                # Create and set table items
                items = [
                    QtWidgets.QTableWidgetItem(feature.type),
                    QtWidgets.QTableWidgetItem(str(chrom)),
                    QtWidgets.QTableWidgetItem(get_id(feature)),
                    QtWidgets.QTableWidgetItem(get_name(feature)),
                    QtWidgets.QTableWidgetItem(get_description(feature))
                ]

                for col, item in enumerate(items):
                    self.tableWidget.setItem(index, col, item)

                mainWindow.checkBoxes.append((chrom, feature, index))

            self.tableWidget.resizeColumnsToContents()
            mainWindow.hide()

            #center on current screen
            center_ui(self)
            self.show()
            self.activateWindow()

            return 0
        except Exception as e:
            show_error("fill_table_nonKegg() in AnnotationWindow", e)

    # this is the connection for the select all checkbox - selects/deselects all the genes in the table
    def select_all_genes(self):
        try:
            if self.select_all_checkbox.isChecked():
                self.tableWidget.selectAll()
            else:
                self.tableWidget.clearSelection()
        except Exception as e:
            show_error("select_all_genes() in AnnotationWindow", e)

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        try:
            GlobalSettings.mainWindow.closeFunction()
        except Exception as e:
            show_error("closeEvent() in AnnotationWindow", e)
            event.accept()