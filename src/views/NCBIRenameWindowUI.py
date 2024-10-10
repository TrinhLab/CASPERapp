from PyQt6 import QtWidgets, QtGui, QtCore, uic
import os

class NCBIRenameWindowUI(QtWidgets.QDialog):  # Changed to QDialog
    def __init__(self, settings, files):
        super(NCBIRenameWindowUI, self).__init__()
        self.settings = settings
        self.files = files
        self.setup_ui()

    def setup_ui(self):
        # Load the UI file
        uic.loadUi(os.path.join(self.settings.get_ui_dir(), 'ncbi_rename_window.ui'), self)
        self.setWindowTitle("Rename Files")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.settings.get_assets_dir(), "cas9image.ico")))

        # Initialize UI components
        self.init_ui_components()

        # Populate the table with files
        self.populate_table()

        # Set up styles
        self.set_styles()

        self.setModal(True)  # Make the dialog modal

    def init_ui_components(self):
        self.rename_table = self.findChild(QtWidgets.QTableWidget, 'rename_table')
        self.submit_button = self.findChild(QtWidgets.QPushButton, 'submit_button')
        self.go_back = self.findChild(QtWidgets.QPushButton, 'go_back')

        # Set up the table
        self.rename_table.setColumnCount(2)
        self.rename_table.setHorizontalHeaderLabels(['Original Filename', 'New Filename'])
        self.rename_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents)

    def set_styles(self):
        header = self.rename_table.horizontalHeader()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)

    def populate_table(self):
        self.rename_table.setRowCount(len(self.files))
        for row, file in enumerate(self.files):
            item = QtWidgets.QTableWidgetItem(file)
            item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
            self.rename_table.setItem(row, 0, item)
            self.rename_table.setCellWidget(row, 1, QtWidgets.QLineEdit())

    def get_new_names(self):
        new_names = []
        for row in range(self.rename_table.rowCount()):
            original = self.rename_table.item(row, 0).text()
            new_name = self.rename_table.cellWidget(row, 1).text()
            if new_name:
                if original.endswith('.gbff'):
                    new_name = new_name if new_name.endswith('.gbff') else f"{new_name}.gbff"
                elif original.endswith('.fna'):
                    new_name = new_name if new_name.endswith('.fna') else f"{new_name}.fna"
            new_names.append(new_name)
        return new_names

    def clear_table(self):
        self.rename_table.setRowCount(0)
