from PyQt5 import QtWidgets, QtCore, uic, QtGui
from PyQt5.QtGui import QStandardItemModel, QStandardItem
import os
from models.GlobalSettings import GlobalSettings1
from utils.ui import scale_ui
from Bio import SeqIO
import io
import csv

class LevelComboDelegate(QtWidgets.QStyledItemDelegate):
    level_changed = QtCore.pyqtSignal()

    def createEditor(self, parent, option, index):
        editor = QtWidgets.QComboBox(parent)
        editor.addItems(["plasmid", "genome", "chromosome"])
        editor.currentTextChanged.connect(self.level_changed.emit)
        return editor

    def setEditorData(self, editor, index):
        value = index.model().data(index, QtCore.Qt.EditRole)
        if value:
            editor.setCurrentText(value)

    def setModelData(self, editor, model, index):
        value = editor.currentText()
        model.setData(index, value, QtCore.Qt.EditRole)
        self.level_changed.emit()

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)

    def paint(self, painter, option, index):
        text = index.data(QtCore.Qt.DisplayRole)
        is_placeholder = index.data(QtCore.Qt.UserRole) == text

        if is_placeholder:
            painter.setPen(QtGui.QColor(128, 128, 128))  # Gray color for placeholder
        else:
            painter.setPen(QtGui.QColor(0, 0, 0))  # Black color for user-selected text

        painter.drawText(option.rect, QtCore.Qt.AlignCenter, text)

        if isinstance(self.parent(), PopulateFNAFiles):
            self.parent().adjust_level_column_width()

class NewNameDelegate(QtWidgets.QStyledItemDelegate):
    def createEditor(self, parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        editor.installEventFilter(self)
        return editor

    def setEditorData(self, editor, index):
        value = index.model().data(index, QtCore.Qt.EditRole) or index.model().data(index, QtCore.Qt.UserRole)
        editor.setText(value)

    def setModelData(self, editor, model, index):
        value = editor.text()
        placeholder = model.data(index, QtCore.Qt.UserRole)
        if value and value != placeholder:
            model.setData(index, value, QtCore.Qt.EditRole)
            model.setData(index, value, QtCore.Qt.DisplayRole)
        else:
            model.setData(index, "", QtCore.Qt.EditRole)
            model.setData(index, placeholder, QtCore.Qt.DisplayRole)

    def paint(self, painter, option, index):
        text = index.data(QtCore.Qt.DisplayRole) or index.data(QtCore.Qt.UserRole)
        if text == index.data(QtCore.Qt.UserRole):
            painter.setPen(QtGui.QColor(128, 128, 128))  # Gray color for placeholder
        else:
            painter.setPen(QtGui.QColor(0, 0, 0))  # Black color for entered text
        painter.drawText(option.rect, QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter, text)

    def eventFilter(self, obj, event):
        if isinstance(obj, QtWidgets.QLineEdit) and event.type() == QtCore.QEvent.KeyPress:
            if event.key() in (QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete):
                if not obj.text():
                    index = self.parent().currentIndex()
                    placeholder = index.data(QtCore.Qt.UserRole)
                    self.parent().model().setData(index, "", QtCore.Qt.EditRole)
                    self.parent().model().setData(index, placeholder, QtCore.Qt.DisplayRole)
                    return True
        return super().eventFilter(obj, event)

# window size 750 x 520
class PopulateFNAFiles(QtWidgets.QMainWindow):
    def __init__(self, settings):
        super(PopulateFNAFiles, self).__init__()
        self.settings = settings
        uic.loadUi(os.path.join(self.settings.get_ui_dir(), 'custom_combine_file_name_selector1.ui'), self)
        self.file_watcher = None
        
        self.available_files_model = QStandardItemModel()
        self.available_files_model.setHorizontalHeaderLabels(['Available FNA Files'])
        self.fna_files_unselected_table.setModel(self.available_files_model)
        
        self.selected_files_model = QStandardItemModel()
        self.selected_files_model.setHorizontalHeaderLabels(['Selected FNA Files', 'New Name', 'Level'])
        self.fna_files_selected_table.setModel(self.selected_files_model)
        self.fna_files_unselected_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.fna_files_selected_table.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked | QtWidgets.QAbstractItemView.EditKeyPressed)
        self.fna_files_selected_table.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        
        # Set the delegate for the "Level" column
        level_delegate = LevelComboDelegate(self)
        level_delegate.level_changed.connect(self.adjust_level_column_width)
        self.fna_files_selected_table.setItemDelegateForColumn(2, level_delegate)
        
        # Set the delegate for the "New Name" column
        new_name_delegate = NewNameDelegate(self)
        self.fna_files_selected_table.setItemDelegateForColumn(1, new_name_delegate)
        
        self.populate_fna_files()
        scale_ui(self, custom_scale_width=1150, custom_scale_height=600)
        
        self.setup_tables()
        self.setup_column_resize()
        
        # Connect buttons
        self.select_button.clicked.connect(self.select_files)
        self.deselect_button.clicked.connect(self.deselect_files)
        self.merge_button.clicked.connect(self.merge_selected_files)
        self.resize_columns()
        
        # Connect item changed signal
        self.selected_files_model.itemChanged.connect(self.on_item_changed)
        self.fna_files_selected_table.installEventFilter(self)
        self.fna_files_unselected_table.installEventFilter(self)

    def populate_fna_files(self):
        try:
            fna_folder = os.path.join(self.settings.CSPR_DB, "FNA")
            if not os.path.exists(fna_folder):
                raise FileNotFoundError(f"FNA folder not found at {fna_folder}")

            fna_files = [f for f in os.listdir(fna_folder) if f.endswith('.fna') and os.path.isfile(os.path.join(fna_folder, f))]
            fna_files.sort(key=str.lower)

            self.available_files_model.clear()
            self.available_files_model.setHorizontalHeaderLabels(['Available FNA Files'])

            for file in fna_files:
                item = QStandardItem(file)
                item.setEditable(False)
                self.available_files_model.appendRow(item)

            self.fna_files_unselected_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.fna_files_unselected_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.fna_files_unselected_table.horizontalHeader().setStretchLastSection(True)
            self.fna_files_unselected_table.verticalHeader().setVisible(False)

            # Set up a file system watcher to monitor changes in the FNA folder
            if self.file_watcher is None:
                self.file_watcher = QtCore.QFileSystemWatcher([fna_folder])
                self.file_watcher.directoryChanged.connect(self.populate_fna_files)

            self.setup_tables()

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"An error occurred while populating FNA files: {e}")

    def select_files(self):
        selected_indexes = self.fna_files_unselected_table.selectionModel().selectedRows()
        selected_indexes = sorted(selected_indexes, key=lambda x: x.row(), reverse=True)
        for index in selected_indexes:
            items = self.available_files_model.takeRow(index.row())
            name_item = items[0]
            file_name = name_item.text()
            name_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)

            # Extract full organism name and level from the file
            file_path = os.path.join(self.settings.CSPR_DB, "FNA", file_name)
            organism_name = ""
            level = "chromosome"  # Default level
            try:
                with open(file_path, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        organism_name = record.description.split(' ', 1)[1]
                        if 'chromosome' in record.description.lower():
                            level = "chromosome"
                        elif 'plasmid' in record.description.lower():
                            level = "plasmid"
                        else:
                            level = "genome"
                        break  # We only need the first record
            except Exception as e:
                print(f"Error extracting organism name from {file_name}: {e}")

            new_name_item = QStandardItem()
            new_name_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable)
            new_name_item.setData(organism_name if organism_name else "Enter new name", QtCore.Qt.UserRole)

            level_item = QStandardItem(level)
            level_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable)
            level_item.setData(level, QtCore.Qt.UserRole)  # Set the placeholder text as UserRole
            self.selected_files_model.appendRow([name_item, new_name_item, level_item])
        self.selected_files_model.sort(0)
        self.resize_columns()
        self.fna_files_selected_table.setFocus()
        self.adjust_level_column_width()
        self.fna_files_selected_table.update()  # Force update to show placeholder text
        self.adjust_column_widths()

    def deselect_files(self):
        selection_model = self.fna_files_selected_table.selectionModel()
        selected_indexes = selection_model.selectedIndexes()
        
        if not selected_indexes:
            return

        # Group indexes by row
        rows_to_remove = {}
        for index in selected_indexes:
            if index.row() not in rows_to_remove:
                rows_to_remove[index.row()] = []
            rows_to_remove[index.row()].append(index.column())

        # Sort rows in descending order to avoid index shifting when removing
        for row in sorted(rows_to_remove.keys(), reverse=True):
            file_name = self.selected_files_model.item(row, 0).text()
            
            self.selected_files_model.removeRow(row)
            
            new_item = QStandardItem(file_name)
            new_item.setEditable(False)
            self.available_files_model.appendRow(new_item)

        self.available_files_model.sort(0)
        self.resize_columns()
        self.adjust_level_column_width()

        if self.selected_files_model.rowCount() == 0:
            # Reset column widths to default values
            self.fna_files_selected_table.setColumnWidth(0, 200)  # "Selected FNA Files"
            self.fna_files_selected_table.setColumnWidth(1, 470)  # "New Name"
            self.fna_files_selected_table.setColumnWidth(2, 30)  # "Level"

    def merge_selected_files(self):
        if self.selected_files_model.rowCount() == 0:
            QtWidgets.QMessageBox.warning(self, "No Selection", "Please select files to merge.")
            return

        selected_files = []
        for row in range(self.selected_files_model.rowCount()):
            file_name = self.selected_files_model.item(row, 0).text()
            new_name = self.selected_files_model.item(row, 1).text()
            level = self.selected_files_model.item(row, 2).text()
            selected_files.append((file_name, new_name, level))
        
        output_file_name, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Merged File", "", "FNA Files (*.fna)")
        if not output_file_name:
            return

        with open(output_file_name, 'w') as outfile:
            for file_index, (fname, new_name, level) in enumerate(selected_files, 1):
                file_path = os.path.join(self.settings.CSPR_DB, "FNA", fname)
                for record_index, record in enumerate(SeqIO.parse(file_path, "fasta"), 1):
                    organism = new_name if new_name else record.description.split(' ', 1)[1]
                    entry_type = level if level else ('chromosome' if 'chromosome' in record.description.lower() else 'plasmid' if 'plasmid' in record.description.lower() else 'genome')
                    
                    new_id = f"{organism}_{file_index}_{entry_type}_{record_index}"
                    outfile.write(f">{new_id}\n{record.seq}\n")

        # Remove merged files from the "File Selected" table
        self.selected_files_model.removeRows(0, self.selected_files_model.rowCount())

        # Update the UI
        self.resize_columns()
        self.adjust_level_column_width()
        self.adjust_column_widths()

        QtWidgets.QMessageBox.information(self, "Merge Complete", f"Files merged and saved to {output_file_name}")

    def resize_columns(self):
        self.fna_files_unselected_table.resizeColumnsToContents()
        self.fna_files_selected_table.resizeColumnsToContents()
        
        # Set the stretch factor for the first column of both tables
        self.fna_files_unselected_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.fna_files_selected_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        
        # Set other columns to be interactive (user can resize them)
        for i in range(1, self.fna_files_selected_table.model().columnCount()):
            self.fna_files_selected_table.horizontalHeader().setSectionResizeMode(i, QtWidgets.QHeaderView.Interactive)
        
        # Adjust the "Level" column width
        self.adjust_level_column_width()

    def on_item_changed(self, item):
        if item.column() == 1:  # Only "New Name" column
            self.selected_files_model.blockSignals(True)  # Block signals to prevent recursion
            clipboard = QtWidgets.QApplication.clipboard()
            selected_indexes = self.fna_files_selected_table.selectionModel().selectedIndexes()
            new_text = item.text()
            
            # Update all selected items in the "New Name" column
            for index in selected_indexes:
                if index.column() == 1:
                    self.selected_files_model.setItem(index.row(), 1, QStandardItem(new_text))
            
            clipboard.setText(new_text)
            self.selected_files_model.blockSignals(False)  # Unblock signals

        # Handle the "Level" column
        elif item.column() == 2:  # "Level" column
            self.adjust_level_column_width()

    def eventFilter(self, obj, event):
        if obj in (self.fna_files_selected_table, self.fna_files_unselected_table):
            if event.type() == QtCore.QEvent.KeyPress:
                key = event.key()
                modifiers = event.modifiers()
                
                if key == QtCore.Qt.Key_C and modifiers == QtCore.Qt.ControlModifier:
                    self.copy_selected_cells()
                    return True
                elif key == QtCore.Qt.Key_V and modifiers == QtCore.Qt.ControlModifier:
                    self.paste_to_selected_cells()
                    return True
                elif key in (QtCore.Qt.Key_Up, QtCore.Qt.Key_Down, QtCore.Qt.Key_Left, QtCore.Qt.Key_Right):
                    return self.handle_arrow_keys(obj, key)
                elif key in (QtCore.Qt.Key_Delete, QtCore.Qt.Key_Backspace):
                    if obj == self.fna_files_selected_table:
                        self.clear_new_name_entries()
                        return True
                elif key == QtCore.Qt.Key_Return or key == QtCore.Qt.Key_Enter:
                    if obj == self.fna_files_unselected_table:
                        self.select_files()
                        return True

        return super().eventFilter(obj, event)

    def copy_selected_cells(self):
        selection = self.fna_files_selected_table.selectedIndexes()
        if selection:
            for index in selection:
                if index.column() == 1:  # "New Name" column
                    text = index.data()
                    QtWidgets.QApplication.clipboard().setText(text)
                    break  # Copy only the first selected "New Name" cell

    def paste_to_selected_cells(self):
        selection = self.fna_files_selected_table.selectedIndexes()
        if selection:
            clipboard = QtWidgets.QApplication.clipboard()
            text = clipboard.text()
            
            for index in selection:
                if index.column() == 1:  # "New Name" column
                    self.selected_files_model.setItem(index.row(), 1, QStandardItem(text))

    def handle_arrow_keys(self, obj, key):
        current = obj.currentIndex()
        model = obj.model()
        if key in (QtCore.Qt.Key_Up, QtCore.Qt.Key_Down):
            new_row = (current.row() - 1) % model.rowCount() if key == QtCore.Qt.Key_Up else (current.row() + 1) % model.rowCount()
            new_index = model.index(new_row, current.column())
        elif key in (QtCore.Qt.Key_Left, QtCore.Qt.Key_Right) and obj == self.fna_files_selected_table:
            new_col = (current.column() - 1) % model.columnCount() if key == QtCore.Qt.Key_Left else (current.column() + 1) % model.columnCount()
            new_index = model.index(current.row(), new_col)
        else:
            return False
        obj.setCurrentIndex(new_index)
        return True

    def clear_new_name_entries(self):
        selection = self.fna_files_selected_table.selectionModel()
        if selection.hasSelection():
            for index in selection.selectedIndexes():
                if index.column() == 1:  # "New Name" column
                    placeholder = self.selected_files_model.data(index, QtCore.Qt.UserRole)
                    self.selected_files_model.setData(index, "", QtCore.Qt.EditRole)
                    self.selected_files_model.setData(index, placeholder, QtCore.Qt.DisplayRole)
        self.fna_files_selected_table.update()

    def adjust_level_column_width(self):
        level_column = 2  # Assuming "Level" is the third column (index 2)
        min_width = 60  # Minimum width for the "Level" column
        max_width = 100  # Maximum width for the "Level" column
        padding = 20  # Padding for the content

        # Find the maximum width needed for the items in the "Level" column
        max_content_width = 0
        for row in range(self.selected_files_model.rowCount()):
            index = self.selected_files_model.index(row, level_column)
            item_width = self.fna_files_selected_table.fontMetrics().width(index.data()) + padding
            max_content_width = max(max_content_width, item_width)

        # Consider the header width as well
        header_width = self.fna_files_selected_table.fontMetrics().width(self.selected_files_model.headerData(level_column, QtCore.Qt.Horizontal)) + padding
        max_content_width = max(max_content_width, header_width)

        # Set the column width, but allow user to adjust between min and max
        new_width = max(min(max_content_width, max_width), min_width)
        current_width = self.fna_files_selected_table.columnWidth(level_column)

        if new_width != current_width:
            self.fna_files_selected_table.setColumnWidth(level_column, new_width)

    def setup_column_resize(self):
        header = self.fna_files_selected_table.horizontalHeader()
        header.sectionResized.connect(self.on_column_resized)
        header.setStretchLastSection(False)
        
        # Set all columns to Interactive mode
        for i in range(3):
            self.fna_files_selected_table.horizontalHeader().setSectionResizeMode(i, QtWidgets.QHeaderView.Interactive)
        
        self.adjust_column_widths()

    def on_column_resized(self, logical_index, old_size, new_size):
        if logical_index == 2:  # "Level" column
            min_width = 60
            max_width = 100
            if min_width <= new_size <= max_width:
                return  # No need to adjust, size is within acceptable range
            elif new_size < min_width:
                new_size = min_width
            elif new_size > max_width:
                new_size = max_width
            
            self.fna_files_selected_table.horizontalHeader().blockSignals(True)
            self.fna_files_selected_table.setColumnWidth(logical_index, new_size)
            self.fna_files_selected_table.horizontalHeader().blockSignals(False)

    def setup_tables(self):
        # Set up the unselected files table
        self.fna_files_unselected_table.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.fna_files_unselected_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.fna_files_unselected_table.horizontalHeader().setStretchLastSection(True)
        self.fna_files_unselected_table.verticalHeader().setVisible(False)

        # Set up the selected files table
        self.fna_files_selected_table.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.fna_files_selected_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.fna_files_selected_table.horizontalHeader().setStretchLastSection(True)
        self.fna_files_selected_table.verticalHeader().setVisible(False)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.adjust_column_widths()

    def adjust_column_widths(self):
        total_width = self.fna_files_selected_table.width()
        
        # Set the "Selected FNA Files" column to 20% of the total width
        self.fna_files_selected_table.setColumnWidth(0, int(total_width * 0.20))
        
        # Set the "New Name" column to 60% of the total width
        self.fna_files_selected_table.setColumnWidth(1, int(total_width * 0.70))
        
        # Set the "Level" column to 20% of the total width
        self.fna_files_selected_table.setColumnWidth(2, int(total_width * 0.10))
        
        # Adjust the "Level" column width
        self.adjust_level_column_width()

