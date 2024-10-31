from PyQt6 import QtWidgets, QtGui, QtCore, uic
from PyQt6.QtWidgets import QWidget, QVBoxLayout
import os
from typing import Optional

class NCBIWindowView(QtWidgets.QMainWindow):
    initialization_complete = QtCore.pyqtSignal()  # New signal

    def __init__(self, settings):
        super(NCBIWindowView, self).__init__()
        self.settings = settings
        self.logger = settings.get_logger()
        self.progress_bars = {}
        self.progress_labels = {}
        self._is_initialized = False  # Track initialization state
        self._setup_basic_ui()  # Only do basic initialization first

    def _setup_basic_ui(self):
        """Initial minimal setup to show the window quickly"""
        try:
            uic.loadUi(os.path.join(self.settings.get_ui_dir_path(), "ncbi_window_v2.ui"), self)
            
            QtCore.QTimer.singleShot(100, self._complete_initialization)
            
        except Exception as e:
            self.logger.error(f"Error in basic UI setup: {str(e)}")
            raise

    def _complete_initialization(self):
        """Complete the full initialization of UI components"""
        try:
            if self._is_initialized:
                return

            # Initialize all UI components
            self._init_ui_components()
            
            self._is_initialized = True
            self.logger.debug("NCBI Window initialization completed")
            
            # Emit signal after everything is initialized
            self.initialization_complete.emit()
            
        except Exception as e:
            self.logger.error(f"Error in complete initialization: {str(e)}")
            raise

    def _init_ui_components(self) -> None:
        """Initialize all UI components at once instead of using timers"""
        try:
            self._init_grpStep1()
            self._init_grpStep2()
            self._init_grpStep3()
        except Exception as e:
            self.logger.error(f"Error in _init_ui_components: {str(e)}")
            raise

    def _init_grpStep1(self) -> None:
        try:
            self.line_edit_organism = self._find_widget("ledOrganism", QtWidgets.QLineEdit)
            self.line_edit_strain = self._find_widget("ledStrain", QtWidgets.QLineEdit)
            self.line_edit_max_results = self._find_widget("ledMaxResults", QtWidgets.QLineEdit)
            self.check_box_complete_genomes_only = self._find_widget("chkCompleteGenomesOnly", QtWidgets.QCheckBox)
            
            # Set default values
            self.line_edit_max_results.setText("100")
            
        except Exception as e:
            self.logger.error(f"Error initializing Step 1: {str(e)}")

    def _init_grpStep2(self) -> None:
        try:
            self.push_button_search = self._find_widget("pbtnSearch", QtWidgets.QPushButton)
            self.check_box_select_all_rows = self._find_widget("chkSelectAllRows", QtWidgets.QCheckBox)
            self.table_ncbi_results = self._find_widget("tblNCBIResults", QtWidgets.QTableView)
        except Exception as e:
            self.logger.error(f"Error initializing Step 2: {str(e)}")

    def _init_grpStep3(self) -> None:
        try:
            self.radio_button_collections_refseq = self._find_widget("rbtnCollectionsRefSeq", QtWidgets.QRadioButton)
            self.radio_button_collections_genbank = self._find_widget("rbtnCollectionsGenBank", QtWidgets.QRadioButton)
            self.check_box_file_types_fna = self._find_widget("chkFileTypesFNA", QtWidgets.QCheckBox)
            self.check_box_file_types_gbff = self._find_widget("chkFileTypesGBFF", QtWidgets.QCheckBox)
            self.push_button_download_files = self._find_widget("pbtnDownloadFiles", QtWidgets.QPushButton)
            self.progress_bar_download_files = self._find_widget("pbDownloadFiles", QtWidgets.QProgressBar)
            self.label_download_files_status = self._find_widget("lblDownloadFilesStatus", QtWidgets.QLabel)
            
            # Set initial states
            self.progress_bar_download_files.setValue(0)
            self.radio_button_collections_refseq.setChecked(True)
            self.check_box_file_types_fna.setChecked(True)
            
        except Exception as e:
            self.logger.error(f"Error initializing Step 3: {str(e)}")

    def populate_ncbi_table(self, model):
        """Populate the table with data and set up proper row selection"""
        self.table_ncbi_results.setModel(model)
        
        # Set selection behavior to select entire rows
        self.table_ncbi_results.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_ncbi_results.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
        
        # Disable cell editing
        self.table_ncbi_results.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        
        # Enable sorting
        self.table_ncbi_results.setSortingEnabled(True)
        
        # Set the horizontal header to resize mode
        header = self.table_ncbi_results.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        
        # Enable horizontal scrolling
        self.table_ncbi_results.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table_ncbi_results.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Adjust row height to fit content
        self.table_ncbi_results.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        
        # Center-align the header text
        self.table_ncbi_results.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        
        # Set word wrap for all cells
        self.table_ncbi_results.setWordWrap(True)
        
        # Adjust the size of the table to fit its contents
        self.table_ncbi_results.resizeColumnsToContents()
        self.table_ncbi_results.resizeRowsToContents()
        
        # Set a minimum width for each column based on header text
        font_metrics = QtGui.QFontMetrics(self.table_ncbi_results.font())
        for i in range(self.table_ncbi_results.model().columnCount()):
            header_text = self.table_ncbi_results.model().headerData(i, QtCore.Qt.Orientation.Horizontal)
            width = font_metrics.horizontalAdvance(header_text) + 20  # Add some padding
            self.table_ncbi_results.setColumnWidth(i, max(width, self.table_ncbi_results.columnWidth(i)))
        
        # Ensure the last column doesn't stretch
        header.setStretchLastSection(False)
        
        # Set focus policy to enable keyboard selection
        self.table_ncbi_results.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

    def get_search_parameters(self):
        return {
            'organism': self.line_edit_organism.text(),
            'strain': self.line_edit_strain.text(),
            'max_results': self.line_edit_max_results.text(),
            'complete_genomes_only': self.check_box_complete_genomes_only.isChecked(),
            'refseq': self.radio_button_collections_refseq.isChecked(),
            'genbank': self.radio_button_collections_genbank.isChecked(),
            'fna': self.check_box_file_types_fna.isChecked(),
            'gbff': self.check_box_file_types_gbff.isChecked()
        }

    def get_selected_rows(self):
        return self.table_ncbi_results.selectionModel().selectedRows()

    def reset_progress(self):
        self.progress_bar_download_files.setValue(0)
        self.set_download_files_status_label("Ready to download")

    def set_progress(self, value):
        self.progress_bar_download_files.setValue(value)

    def set_download_files_status_label(self, text):
        self.label_download_files_status.setText(text)

    # def clear_form(self):
    #     self.organism_line_edit.clear()
    #     self.infra_name_line_edit.clear()
    #     self.ret_max_line_edit.setText("100")
    #     self.yes_box.setChecked(False)
    #     self.refseq_checkbox.setChecked(True)
    #     self.genbank_checkbox.setChecked(False)
    #     self.fna_checkbox.setChecked(True)
    #     self.gbff_checkbox.setChecked(True)

    def clear_table(self):
        source_model = self.table_ncbi_results.model().sourceModel()
        if source_model:
            source_model.clear()
        self.table_ncbi_results.reset()
    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.settings.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget
