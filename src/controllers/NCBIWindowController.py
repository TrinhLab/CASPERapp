from PyQt6 import QtWidgets, QtCore, QtGui
from utils.ui import show_error, show_message, position_window
from models.NCBIWindowModel import NCBIWindowModel, PandasModel, CustomProxyModel
from views.NCBIWindowView import NCBIWindowView
import os
import time
from controllers.NCBIRenameWindowController import NCBIRenameWindowController
from ftplib import FTP
from urllib.parse import urlparse
import socket
import pandas as pd

class NCBIWindowController:
    def __init__(self, settings):
        self.settings = settings
        try:
            self.logger = self.settings.get_logger()
            self.model = NCBIWindowModel(settings)
            self.view = NCBIWindowView(settings)
            
            self._init_ui()
            self.setup_connections()
        except Exception as e:
            show_error(self.settings, "Error initializing NCBIWindowController", str(e))

    def setup_connections(self):
        try:
            self.view.push_button_search.clicked.connect(self.search_ncbi)
            self.view.push_button_download_files.clicked.connect(self.download_files_wrapper)
            self.view.check_box_select_all_rows.clicked.connect(self.select_all_rows_in_table)
            self.view.radio_button_collections_genbank.toggled.connect(self.is_checked_GenBank_radio_button)

        except Exception as e:
            self.logger.error(f"Error setting up connections: {str(e)}", exc_info=True)
            show_error(self.settings, "Error setting up connections", str(e))

    def _init_ui(self):
        # self.view.ret_max_line_edit.setText("100")
        # self.view.progressBar.setValue(0)
        pass

    def search_ncbi(self):
        try:
            self.view.reset_progress()
            search_params = self.view.get_search_parameters()
            
            # Set default values if fields are empty
            if not search_params['organism'].strip():
                search_params['organism'] = "Escherichia coli"
            if not search_params['strain'].strip():
                search_params['strain'] = "K-12"

            self.logger.info(f"Querying NCBI with parameters: {search_params}")
            
            self.df = self.model.search_ncbi(search_params)
            
            if self.df.empty:
                print("No results found")
                self.logger.warning("No results found for the given search parameters.")
                show_message(12, QtWidgets.QMessageBox.Icon.Warning, "No Results", "No results found for the given search parameters.")
            else:
                print(f"Query returned {len(self.df)} results")
                self.logger.info(f"Query returned {len(self.df)} results")
                self.populate_ncbi_table(self.df)
            
            self.view.activateWindow()
        except Exception as e:
            print(f"Error in query_db: {str(e)}")
            self.logger.error(f"Error in query_db: {str(e)}", exc_info=True)
            show_error(self.settings, "Error in query_db", e)

    def populate_ncbi_table(self, data):
        try:
            self.logger.info(f"Populating table with {len(data)} rows")
            self.pandas_model = PandasModel(data)
            self.proxy_model = CustomProxyModel(self.view)
            self.proxy_model.setSourceModel(self.pandas_model)
            self.view.populate_ncbi_table(self.proxy_model)
            self.logger.info("Table populated successfully")
        except Exception as e:
            self.logger.error(f"Error in populate_table: {str(e)}", exc_info=True)
            show_error(self.settings, "Error in populate_table", e)

    def download_files_wrapper(self):
        try:
            self.view.reset_progress()

            if self.model.df.empty:
                show_message(12, QtWidgets.QMessageBox.Icon.Warning, "No Query Results", "Please run an NCBI query to fill the table with results to choose from!")
                return

            indices = self.view.get_selected_rows()
            if not indices:
                show_message("No Rows Selected", "Please select rows from the table!")
                return

            thread_count = QtCore.QThreadPool.globalInstance().maxThreadCount()
            if len(indices) > thread_count:
                show_message("Too Many Selections!", 
                             f"You only have {thread_count} threads available to download with.\n\nPlease select {thread_count} or fewer rows.")
                return

            if not self.view.radio_button_collections_refseq.isChecked() and not self.view.radio_button_collections_genbank.isChecked():
                show_message("No File Type Selected", "No file type selected. Please select the file types you want to download!")
                return

            self.download_files(indices)
        except Exception as e:
            self.logger.error(f"Error in download_files_wrapper: {str(e)}", exc_info=True)
            show_error(self.settings, "Error in download_files_wrapper", e)

    def download_files(self, selected_rows):
        try:
            self.view.reset_progress()
            self.total_files = len(selected_rows)
            self.view.progress_bar_download_files.setMaximum(self.total_files)
            self.view.set_download_files_status_label("Preparing to download...")

            self.model.clear_downloaded_files()
            self.download_threads = []
            self.completed_threads = 0
            self.unavailable_files = []

            for index in selected_rows:
                id = self.proxy_model.data(self.proxy_model.index(index.row(), 0))
                species_name = self.proxy_model.data(self.proxy_model.index(index.row(), 1))
                strain = self.proxy_model.data(self.proxy_model.index(index.row(), 2))
                self.logger.info(f"Processing ID: {id}")
                
                url = self.model.get_download_url(id, self.view.radio_button_collections_genbank.isChecked())
                self.logger.info(f"Download URL for ID {id}: {url}")
                if not url:
                    self.logger.warning(f"No download URL found for ID: {id}")
                    self.unavailable_files.append((species_name, strain))
                    self.on_thread_completed()  # Increment completed threads for unavailable files
                    continue

                downloader = DownloadThread(self, url, id, species_name, strain, self.view.check_box_file_types_fna.isChecked(), self.view.check_box_file_types_gbff.isChecked())
                downloader.finished.connect(self.on_download_finished)
                downloader.progress_updated.connect(self.update_progress)
                downloader.status_updated.connect(self.update_status)
                downloader.all_completed.connect(self.on_thread_completed)
                self.download_threads.append(downloader)
                downloader.start()

            if not self.download_threads and not self.unavailable_files:
                show_message(12, QtWidgets.QMessageBox.Icon.Information, "No Downloads", "No valid files selected for download.")
                return

        except Exception as e:
            self.logger.error(f"Error in download_files: {str(e)}", exc_info=True)
            show_error(self.settings, "Error in download_files", e)

    def on_thread_completed(self):
        self.completed_threads += 1
        self.view.set_progress(self.completed_threads)
        
        if self.completed_threads == self.total_files:
            self.view.set_download_files_status_label("Download Complete.<br>Press Back to go back to the Main window.")

            if self.unavailable_files:
                self.show_unavailable_files_warning()
            
            if self.model.files:
                self.rename_files(self.model.files)
            elif not self.unavailable_files:
                show_message(12, QtWidgets.QMessageBox.Icon.Warning, "No Files Downloaded", "No files were downloaded. Please check your selection and try again.")
            
    def show_unavailable_files_warning(self):
        warning_text = "The following files were not available for download:\n\n"
        for species_name, strain in self.unavailable_files:
            warning_text += f"{species_name} ({strain})\n"
        
        show_message(
            12,
            QtWidgets.QMessageBox.Icon.Warning,
            "Unavailable Files",
            warning_text,
            QtWidgets.QMessageBox.StandardButton.Close
        )

    def update_progress(self, id, value, max_value):
        current_value = self.view.progress_bar_download_files.value()
        self.view.set_progress(current_value + 1)

    def update_status(self, status):
        self.view.set_download_files_status_label(status)

    def on_download_finished(self, success):
        if success:
            current_value = self.view.progress_bar_download_files.value()
            self.view.set_progress(current_value + 1)
        else:
            self.view.set_download_files_status_label("Download Failed")
            show_message(12, QtWidgets.QMessageBox.Icon.Warning, "Download Failed", "Failed to download one or more files. Please check your internet connection and try again.")

    def rename_files(self, files):
        reply = QtWidgets.QMessageBox.question(self.view, "Rename Files", 
                                               "Would you like to rename the downloaded files?",
                                               QtWidgets.QMessageBox.StandardButton.Yes | 
                                               QtWidgets.QMessageBox.StandardButton.No)
        
        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            rename_controller = NCBIRenameWindowController(self.settings, files, self)
            if rename_controller.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.on_rename_complete()
        else:
            self.on_rename_complete()

    def on_rename_complete(self):
        try:
            # self.update_main_window()
            self.view.set_download_files_status_label("Download and renaming complete.<br>Press Back to go back to the Main window.")
        except Exception as e:
            self.logger.error(f"Error in on_rename_complete: {str(e)}", exc_info=True)
            show_error(self.settings, "Error after renaming", str(e))

    def update_main_window(self):
        try:
            if hasattr(self.main_window, 'fill_annotation_dropdown'):
                self.main_window.fill_annotation_dropdown()
            else:
                self.logger.warning("MainWindowController does not have fill_annotation_dropdown method")
        except Exception as e:
            self.logger.error(f"Error updating main window: {str(e)}", exc_info=True)
            show_error(self.settings, "Error updating main window", str(e))

    def select_all_rows_in_table(self):
        if self.view.check_box_select_all_rows.isChecked():
            self.view.table_ncbi_results.selectAll()
        else:
            self.view.table_ncbi_results.clearSelection()
        
        # Force the view to update
        self.view.table_ncbi_results.viewport().update()

    def is_checked_GenBank_radio_button(self, checked):
        if checked:
            show_message("Warning!", 
                         "The GenBank collection may contain poorly or partially annotated annotation files. "
                         "We highly recommend using the RefSeq collection if it is available.")

class DownloadThread(QtCore.QThread):
    finished = QtCore.pyqtSignal(bool)
    progress_updated = QtCore.pyqtSignal(int, int, int)
    status_updated = QtCore.pyqtSignal(str)
    all_completed = QtCore.pyqtSignal()  # New signal to indicate all operations are complete

    def __init__(self, controller, url, id, species_name, strain, download_fna, download_gbff):
        super().__init__()
        self.controller = controller
        self.url = url
        self.id = id
        self.species_name = species_name
        self.strain = strain
        self.download_fna = download_fna
        self.download_gbff = download_gbff

    def run(self):
        try:
            parsed_url = urlparse(self.url)
            ftp_host = parsed_url.netloc
            ftp_path = parsed_url.path

            self.controller.logger.info(f"Attempting to connect to FTP server: {ftp_host}")
            
            # Resolve the IP address
            try:
                ip_address = socket.gethostbyname(ftp_host)
                self.controller.logger.info(f"Resolved IP address: {ip_address}")
            except socket.gaierror as e:
                self.controller.logger.error(f"Failed to resolve hostname: {ftp_host}. Error: {str(e)}")
                self.finished.emit(False)
                return

            ftp = FTP(ftp_host)
            ftp.login()
            ftp.cwd(ftp_path)
            ftp.set_pasv(True)  # Use passive mode
            ftp.voidcmd('TYPE I')  # Set binary mode
            
            self.controller.logger.info(f"Successfully connected to FTP server: {ftp_host}")
            
            files_to_download = []
            if self.download_fna:
                fna_files = [f for f in ftp.nlst() if f.endswith('_genomic.fna.gz')]
                self.controller.logger.info(f"Found {len(fna_files)} FNA files")
                files_to_download.extend(fna_files)
            if self.download_gbff:
                gbff_files = [f for f in ftp.nlst() if f.endswith('_genomic.gbff.gz')]
                self.controller.logger.info(f"Found {len(gbff_files)} GBFF files")
                files_to_download.extend(gbff_files)

            self.controller.logger.info(f"Total files to download: {len(files_to_download)}")

            total_size = 0
            for file in files_to_download:
                try:
                    total_size += ftp.size(file)
                except Exception as e:
                    self.controller.logger.warning(f"Could not get size for file {file}: {str(e)}")

            downloaded_size = 0

            for file in files_to_download:
                self.status_updated.emit(f"Downloading: {file}")
                file_type = 'FNA' if file.endswith('.fna.gz') else 'GBFF'
                output_dir = os.path.join(self.controller.settings.CSPR_DB, file_type)
                
                # Create the output directory if it doesn't exist
                os.makedirs(output_dir, exist_ok=True)
                
                local_filename = os.path.join(output_dir, file)
                
                self.controller.logger.info(f"Downloading file: {file} to {local_filename}")
                
                with open(local_filename, 'wb') as local_file:
                    def callback(data):
                        local_file.write(data)
                        nonlocal downloaded_size
                        downloaded_size += len(data)
                        if total_size > 0:
                            self.progress_updated.emit(self.id, downloaded_size, total_size)

                    ftp.retrbinary(f"RETR {file}", callback)

                self.controller.logger.info(f"Download complete: {file}")
                self.status_updated.emit(f"Decompressing: {file}")

                # Decompress the file
                self.controller.model.decompress_file(local_filename)
                
                # Add the decompressed file to the model's files list
                decompressed_filename = local_filename[:-3]  # Remove .gz extension
                self.controller.model.add_downloaded_file(decompressed_filename)

            ftp.quit()
            self.controller.logger.info(f"All files downloaded and decompressed successfully for ID: {self.id}")
            self.all_completed.emit()  # Emit the new signal when everything is done
            self.finished.emit(True)
        except Exception as e:
            self.controller.logger.error(f"Download error for ID {self.id}: {str(e)}", exc_info=True)
            self.finished.emit(False)
