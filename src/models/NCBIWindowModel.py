import pandas as pd
from PyQt6 import QtCore, QtGui
from Bio import Entrez
from bs4 import BeautifulSoup, XMLParsedAsHTMLWarning
from ftplib import FTP
import gzip
import os
import ssl
import platform
import warnings
import xml.etree.ElementTree as ET
from PyQt6.QtCore import Qt
import socket
from urllib.parse import urlparse

class NCBIWindowModel:
    def __init__(self, settings):
        self.settings = settings
        self.logger = settings.get_logger()
        self.df = pd.DataFrame()
        self.genbank_ftp_dict = {}
        self.refseq_ftp_dict = {}
        self.files = []  # Initialize the files list
        Entrez.email = "your_email@example.com"  # Replace with a valid email
        
        # Make DownloadThread accessible through the model
        self.DownloadThread = DownloadThread

        # Suppress the XMLParsedAsHTMLWarning
        warnings.filterwarnings("ignore", category=XMLParsedAsHTMLWarning)

    def search_ncbi(self, search_params):
        try:
            retmax = int(search_params['max_results']) if search_params['max_results'] else 100
            term = f'"{search_params["organism"]}"[Organism]'
            
            if search_params['complete_genomes_only']:
                term += ' AND "Complete Genome"[Assembly Level]'
            if search_params['strain']:
                term += f' AND "{search_params["strain"]}"[Infraspecific name]'

            self.logger.info(f"Searching NCBI with term: {term}")
            handle = Entrez.esearch(db="assembly", retmax=retmax, term=term)
            search_results = Entrez.read(handle)
            handle.close()

            if not search_results['IdList']:
                self.logger.warning(f"No IDs found for the given search parameters. Term: {term}")
                return pd.DataFrame()

            ids = search_results['IdList']
            self.logger.info(f"Found {len(ids)} IDs")

            handle = Entrez.esummary(db="assembly", id=','.join(ids))
            summary_results = Entrez.read(handle)
            handle.close()

            self.logger.info(f"Retrieved summary results for {len(summary_results['DocumentSummarySet']['DocumentSummary'])} assemblies")
            
            # Log the keys of the first assembly for debugging
            if summary_results['DocumentSummarySet']['DocumentSummary']:
                first_assembly = summary_results['DocumentSummarySet']['DocumentSummary'][0]
                self.logger.debug(f"Keys in first assembly: {first_assembly.keys()}")
                self.logger.debug(f"Attributes of first assembly: {first_assembly.attributes}")

            self._process_ncbi_data(summary_results)

            self.logger.info(f"Processed data into DataFrame with {len(self.df)} rows")
            return self.df

        except Exception as e:
            self.logger.error(f"Error in query_ncbi: {str(e)}", exc_info=True)
            raise

    def _process_ncbi_data(self, summary_results):
        data = []
        for assembly in summary_results['DocumentSummarySet']['DocumentSummary']:
            gb_ftp = assembly.get('FtpPath_GenBank', '')
            rs_ftp = assembly.get('FtpPath_RefSeq', '')
            
            entry = {
                'ID': assembly.attributes['uid'],
                'Species Name': assembly.get('SpeciesName', 'N/A'),
                'Strain': assembly.get('Biosource', {}).get('InfraspeciesList', [{}])[0].get('Sub_value', 'N/A'),
                'Assembly Name': assembly.get('AssemblyName', 'N/A'),
                'RefSeq assembly accession': 'Not Available' if not rs_ftp else assembly.get('AssemblyAccession', 'N/A'),
                'GenBank assembly accession': 'Not Available' if not gb_ftp else assembly.get('AssemblyAccession', 'N/A'),
                'Assembly Status': assembly.get('AssemblyStatus', 'N/A')
            }
            
            if gb_ftp:
                self.genbank_ftp_dict[entry['ID']] = gb_ftp + '/'
            if rs_ftp:
                self.refseq_ftp_dict[entry['ID']] = rs_ftp + '/'

            data.append(entry)

        self.df = pd.DataFrame(data)
        self.df.replace('', 'N/A', inplace=True)

        if self.df.empty:
            self.logger.warning("No data found in NCBI response")
        else:
            self.logger.info(f"Processed {len(self.df)} entries from NCBI response")

        # Log a sample of the processed data
        if not self.df.empty:
            self.logger.debug(f"Sample of processed data:\n{self.df.head().to_string()}")

    def _get_element_text(self, element, tag):
        el = element.find(tag)
        return el.text if el is not None else 'N/A'

    def get_download_url(self, id, use_genbank):
        url = self.genbank_ftp_dict.get(id, '') if use_genbank else self.refseq_ftp_dict.get(id, '')
        if url and not url.startswith('ftp://'):
            url = 'ftp://' + url
        return url.rstrip('/')  # Remove trailing slash if present

    def decompress_file(self, filename):
        try:
            block_size = 65536
            with gzip.open(filename, 'rb') as f_in:
                with open(filename.replace('.gz', ''), 'wb') as f_out:
                    while True:
                        block = f_in.read(block_size)
                        if not block:
                            break
                        else:
                            f_out.write(block)
            os.remove(filename)
        except Exception as e:
            self.logger.error(f"Error in decompress_file: {str(e)}")
            raise

    def get_output_path(self, file_type):
        if platform.system() == "Windows":
            return os.path.join(self.settings.CSPR_DB, file_type)
        else:
            return os.path.join(self.settings.CSPR_DB, file_type)

    def rename_file(self, old_name, new_name, file_type):
        old_path = os.path.join(self.get_output_path(file_type), old_name)
        new_path = os.path.join(self.get_output_path(file_type), new_name)
        
        if not os.path.exists(new_path):
            os.rename(old_path, new_path)
            return True
        else:
            return False

    def add_downloaded_file(self, file_path):
        self.files.append(file_path)

    def clear_downloaded_files(self):
        self.files.clear()

class PandasModel(QtCore.QAbstractTableModel):
    def __init__(self, df=pd.DataFrame(), parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent=parent)
        self._df = df.copy()

    def toDataFrame(self):
        return self._df.copy()

    def headerData(self, section, orientation, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if role != QtCore.Qt.ItemDataRole.DisplayRole:
            return QtCore.QVariant()

        if orientation == QtCore.Qt.Orientation.Horizontal:
            try:
                return self._df.columns.tolist()[section]
            except (IndexError, ):
                return QtCore.QVariant()
        elif orientation == QtCore.Qt.Orientation.Vertical:
            try:
                return self._df.index.tolist()[section]
            except (IndexError, ):
                return QtCore.QVariant()

    def data(self, index, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return QtCore.QVariant()

        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            value = str(self._df.iloc[index.row(), index.column()])
            return QtCore.QVariant(value)
        elif role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
            return QtCore.Qt.AlignmentFlag.AlignCenter
        elif role == QtCore.Qt.ItemDataRole.ForegroundRole:
            value = str(self._df.iloc[index.row(), index.column()])
            if value == "Not Available":
                return QtCore.QVariant(QtGui.QColor(Qt.GlobalColor.red))

        return QtCore.QVariant()

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self._df.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self._df.columns)

    def sort(self, column, order):
        colname = self._df.columns.tolist()[column]
        self.layoutAboutToBeChanged.emit()
        self._df.sort_values(colname, ascending=order == QtCore.Qt.SortOrder.AscendingOrder, inplace=True)
        self._df.reset_index(inplace=True, drop=True)
        self.layoutChanged.emit()

    def clear(self):
        self.beginResetModel()
        self._df = pd.DataFrame()
        self.endResetModel()

class CustomProxyModel(QtCore.QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._filters = {}

    @property
    def filters(self):
        return self._filters

    def setFilter(self, expression, column):
        if expression:
            self.filters[column] = expression
        elif column in self.filters:
            del self.filters[column]
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        for column, expression in self.filters.items():
            text = self.sourceModel().index(source_row, column, source_parent).data()
            regex = QtCore.QRegularExpression(expression, QtCore.QRegularExpression.CaseInsensitiveOption)
            if regex.match(text).hasMatch():
                return False
        return True

class DownloadThread(QtCore.QThread):
    finished = QtCore.pyqtSignal(bool)
    progress_updated = QtCore.pyqtSignal(int, int, int)
    status_updated = QtCore.pyqtSignal(str)
    all_completed = QtCore.pyqtSignal()

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
            ftp.set_pasv(True)
            
            # Set binary mode before any operations
            ftp.voidcmd('TYPE I')
            
            files_to_download = []
            
            # Get list of all files
            all_files = ftp.nlst()
            
            # Process FNA files if requested
            if self.download_fna:
                # Find the main genomic FNA file (should be exactly one)
                genomic_fna = [f for f in all_files 
                             if f.endswith('_genomic.fna.gz') 
                             and not any(x in f for x in ['cds_from', 'rna_from'])]
                
                if genomic_fna:
                    files_to_download.append(genomic_fna[0])
                    self.controller.logger.info(f"Found main genomic FNA file: {genomic_fna[0]}")
            
            # Process GBFF files if requested
            if self.download_gbff:
                gbff_files = [f for f in all_files if f.endswith('_genomic.gbff.gz')]
                files_to_download.extend(gbff_files)
                self.controller.logger.info(f"Found GBFF files: {gbff_files}")

            # Calculate total size with error handling
            total_size = 0
            for file in files_to_download:
                try:
                    size = ftp.size(file)
                    if size is not None:
                        total_size += size
                except Exception as e:
                    self.controller.logger.warning(f"Could not get size for file {file}: {str(e)}")

            downloaded_size = 0

            # Download files
            for file in files_to_download:
                try:
                    self.status_updated.emit(f"Downloading: {file}")
                    file_type = 'FNA' if file.endswith('.fna.gz') else 'GBFF'
                    output_dir = os.path.join(self.controller.settings.CSPR_DB, file_type)
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

                    self.controller.model.decompress_file(local_filename)
                    decompressed_filename = local_filename[:-3]
                    self.controller.model.add_downloaded_file(decompressed_filename)
                    
                except Exception as e:
                    self.controller.logger.error(f"Error downloading file {file}: {str(e)}")
                    continue

            ftp.quit()
            self.controller.logger.info(f"All files downloaded and decompressed successfully for ID: {self.id}")
            self.all_completed.emit()
            self.finished.emit(True)
            
        except Exception as e:
            self.controller.logger.error(f"Download error for ID {self.id}: {str(e)}", exc_info=True)
            self.finished.emit(False)
