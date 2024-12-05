import pandas as pd
from PyQt6 import QtCore, QtGui
from Bio import Entrez
from ftplib import FTP
from PyQt6.QtCore import Qt
import gzip
import os
import platform
import requests
from urllib.parse import urlparse

class NCBIWindowModel:
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
            self.logger = controller.settings.get_logger()
            self.db_path = controller.settings.get_db_path()

        def run(self):
            try:
                parsed_url = urlparse(self.url)
                if parsed_url.scheme == 'ftp':
                    self._download_ftp()
                else:
                    self._download_http()
                self.all_completed.emit()
                self.finished.emit(True)
            except Exception as e:
                self.logger.error(f"Error in download thread: {str(e)}")
                self.finished.emit(False)

        def _download_http(self):
            try:
                self.status_updated.emit(f"Downloading {self.species_name} ({self.strain})")
                response = requests.get(self.url, stream=True)
                response.raise_for_status()
                
                total_size = int(response.headers.get('content-length', 0))
                is_gzipped = self.url.endswith('.gz') or 'gzip=true' in self.url
                is_fna = '.fna.' in self.url.lower() or 'fasta' in self.url.lower()
                
                file_type = 'FNA' if is_fna else 'GBFF'
                extension = '.gz' if is_gzipped else ''
                
                # Create output directory using the current database path
                output_dir = os.path.join(self.db_path, file_type)
                os.makedirs(output_dir, exist_ok=True)
                
                # Determine correct file extension and name
                if is_fna:
                    local_filename = os.path.join(output_dir, f"{self.id}.fna{extension}")
                else:
                    local_filename = os.path.join(output_dir, f"{self.id}.gbff{extension}")
                
                downloaded_size = 0
                with open(local_filename, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded_size += len(chunk)
                            if total_size:
                                progress = (downloaded_size / total_size) * 100
                                self.progress_updated.emit(self.id, progress, 100)
                
                if not os.path.exists(local_filename) or os.path.getsize(local_filename) == 0:
                    raise Exception(f"Downloaded file {local_filename} is empty or does not exist")
                
                if is_gzipped:
                    try:
                        # Read the gzipped content and write to uncompressed file
                        uncompressed_filename = local_filename[:-3]  # Remove .gz
                        with gzip.open(local_filename, 'rb') as f_in:
                            content = f_in.read()
                            if not content:
                                raise Exception("Decompressed content is empty")
                            with open(uncompressed_filename, 'wb') as f_out:
                                f_out.write(content)
                        
                        # Verify uncompressed file
                        if not os.path.exists(uncompressed_filename) or os.path.getsize(uncompressed_filename) == 0:
                            raise Exception(f"Decompressed file is empty")
                        
                        # Remove the gzipped file only if decompression was successful
                        os.remove(local_filename)
                        self.controller.model.add_downloaded_file(uncompressed_filename)
                        
                    except Exception as e:
                        self.logger.error(f"Error decompressing file: {str(e)}")
                        # Try to handle the file as non-gzipped if decompression fails
                        if os.path.exists(local_filename) and os.path.getsize(local_filename) > 0:
                            uncompressed_filename = local_filename[:-3]
                            os.rename(local_filename, uncompressed_filename)
                            self.controller.model.add_downloaded_file(uncompressed_filename)
                        else:
                            raise
                else:
                    self.controller.model.add_downloaded_file(local_filename)
                
                self.status_updated.emit(f"Download complete: {self.species_name}")
                
            except Exception as e:
                self.logger.error(f"Error in HTTP download: {str(e)}")
                raise

        def _download_ftp(self):
            try:
                self.status_updated.emit(f"Downloading {self.species_name} ({self.strain})")
                
                parsed_url = urlparse(self.url)
                ftp = FTP(parsed_url.netloc)
                ftp.login()
                
                file_size = ftp.size(parsed_url.path[1:])
                is_gzipped = self.url.endswith('.gz')
                is_fna = '.fna.' in self.url.lower() or '.fasta.' in self.url.lower()
                
                file_type = 'FNA' if is_fna else 'GBFF'
                extension = '.gz' if is_gzipped else ''
                
                local_filename = os.path.join(
                    self.db_path,
                    file_type,
                    f"{self.id}.{file_type.lower()}{extension}"
                )
                
                os.makedirs(os.path.dirname(local_filename), exist_ok=True)
                
                downloaded_size = 0
                with open(local_filename, 'wb') as f:
                    def callback(data):
                        nonlocal downloaded_size
                        f.write(data)
                        downloaded_size += len(data)
                        if file_size:
                            progress = (downloaded_size / file_size) * 100
                            self.progress_updated.emit(self.id, progress, 100)
                    
                    ftp.retrbinary(f'RETR {parsed_url.path[1:]}', callback)
                
                ftp.quit()
                
                if not os.path.exists(local_filename) or os.path.getsize(local_filename) == 0:
                    raise Exception(f"Downloaded file {local_filename} is empty or does not exist")
                
                if is_gzipped:
                    self._decompress_file(local_filename)
                else:
                    self.controller.model.add_downloaded_file(local_filename)
                    
                self.status_updated.emit(f"Download complete: {self.species_name}")
                
            except Exception as e:
                self.logger.error(f"Error in FTP download: {str(e)}")
                raise

        def _decompress_file(self, filename):
            try:
                with gzip.open(filename, 'rb') as f_in:
                    decompressed_filename = filename[:-3]
                    with open(decompressed_filename, 'wb') as f_out:
                        f_out.write(f_in.read())
                
                if not os.path.exists(decompressed_filename) or os.path.getsize(decompressed_filename) == 0:
                    raise Exception(f"Decompressed file {decompressed_filename} is empty or does not exist")
                
                os.remove(filename)
                self.controller.model.add_downloaded_file(decompressed_filename)
                
            except gzip.BadGzipFile:
                self.logger.warning(f"File {filename} is not actually gzipped, renaming without .gz extension")
                decompressed_filename = filename[:-3]
                os.rename(filename, decompressed_filename)
                self.controller.model.add_downloaded_file(decompressed_filename)

    def __init__(self, settings):
        self.settings = settings
        self.logger = settings.get_logger()
        self.df = pd.DataFrame()
        self.genbank_ftp_dict = {}
        self.refseq_ftp_dict = {}
        self.files = []
        self.current_database = None
        self.current_search_params = {}
        self.download_fna = False
        self.download_gbff = False
        Entrez.email = "your_email@example.com"
        
        self.database_apis = {
            "NCBI GenBank": self._search_ncbi,
            "ENA (European Nucleotide Archive)": self._search_ena,
            "UCSC Genome Browser": self._search_ucsc
        }

    def search_ncbi(self, search_params):
        """Search selected database with given parameters"""
        database = search_params.get('database', "NCBI GenBank")
        search_function = self.database_apis.get(database)
        
        if search_function:
            # Store current database and search parameters
            self.current_database = database
            self.current_search_params = search_params
            return search_function(search_params)
        else:
            self.logger.error(f"Unsupported database: {database}")
            raise ValueError(f"Unsupported database: {database}")

    def _search_ncbi(self, search_params):
        """Original NCBI search implementation"""
        try:
            retmax = int(search_params['max_results']) if search_params['max_results'] else 100
            term = f'"{search_params["organism"]}"[Organism]'
            
            if search_params['complete_genomes_only']:
                term += ' AND "Complete Genome"[Assembly Level]'
            if search_params['strain'] and search_params['strain'].strip():
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

    def _search_ena(self, search_params):
        """Search ENA database"""
        try:
            organism = search_params['organism']
            max_results = int(search_params['max_results'])
            
            # ENA API endpoint for text search
            base_url = "https://www.ebi.ac.uk/ena/portal/api/search"
            
            # Construct query
            query_parts = []
            query_parts.append(f'scientific_name="{organism}"')
            
            if search_params['strain']:
                query_parts.append(f'strain="{search_params["strain"]}"')
            
            if search_params['complete_genomes_only']:
                query_parts.append('assembly_level="complete genome"')
            
            query = " AND ".join(query_parts)
            
            params = {
                'result': 'assembly',
                'query': query,
                'limit': max_results,
                'offset': 0,
                'format': 'json',
                'fields': 'accession,scientific_name,strain,assembly_level,assembly_name,study_accession'
            }
            
            self.logger.info(f"Searching ENA with parameters: {params}")
            response = requests.get(base_url, params=params)
            
            # Log the actual URL being called for debugging
            self.logger.debug(f"ENA API URL: {response.url}")
            
            # Check if the response is successful
            if response.status_code != 200:
                self.logger.error(f"ENA API error: {response.status_code} - {response.text}")
                return pd.DataFrame()
            
            # Parse JSON response
            results = response.json()
            
            if not results:
                self.logger.info("No results returned from ENA")
                return pd.DataFrame()
            
            self.logger.info(f"Raw ENA response: {results[:2]}") 
            
            # Transform ENA results to match NCBI format
            data = []
            for result in results:
                entry = {
                    'ID': result.get('accession', 'N/A'),
                    'Species Name': result.get('scientific_name', 'N/A'),
                    'Strain': result.get('strain', 'N/A'),
                    'Assembly Name': result.get('assembly_name', 'N/A'),
                    'RefSeq assembly accession': 'Not Available',
                    'GenBank assembly accession': result.get('accession', 'N/A'),
                    'Assembly Status': result.get('assembly_level', 'N/A')
                }
                data.append(entry)
            
            self.df = pd.DataFrame(data)
            self.logger.info(f"Found {len(self.df)} results from ENA")
            return self.df
            
        except Exception as e:
            self.logger.error(f"Error searching ENA: {str(e)}")
            raise

    def _search_ucsc(self, search_params):
        """Search UCSC Genome Browser database"""
        try:
            organism = search_params['organism']
            max_results = int(search_params['max_results'])
            
            # Get list of genomes
            base_url = "https://api.genome.ucsc.edu"
            response = requests.get(f"{base_url}/list/ucscGenomes")
            response.raise_for_status()
            
            # Log raw response for debugging
            self.logger.debug(f"UCSC API response: {response.text[:500]}")
            
            data = response.json()
            if 'ucscGenomes' not in data:
                self.logger.warning("No ucscGenomes field in UCSC response")
                return pd.DataFrame()
            
            # Filter genomes by organism name
            matching_genomes = []
            for genome_id, genome in data['ucscGenomes'].items():
                # Search in both organism and scientificName fields
                if (organism.lower() in genome.get('organism', '').lower() or 
                    organism.lower() in genome.get('scientificName', '').lower()):
                    
                    # Get detailed track info
                    track_response = requests.get(f"{base_url}/list/tracks", params={'genome': genome_id})
                    if track_response.status_code == 200:
                        tracks = track_response.json()
                        genome['tracks'] = tracks
                        matching_genomes.append(genome)
                        
                    if len(matching_genomes) >= max_results:
                        break
            
            if not matching_genomes:
                self.logger.info("No results returned from UCSC")
                return pd.DataFrame()
            
            # Transform UCSC results to match NCBI format
            data = []
            for genome in matching_genomes:
                entry = {
                    'ID': genome.get('genome', 'N/A'),
                    'Species Name': genome.get('scientificName', genome.get('organism', 'N/A')),
                    'Strain': genome.get('description', 'N/A'),
                    'Assembly Name': genome.get('sourceName', 'N/A'),
                    'RefSeq assembly accession': 'Not Available',  # UCSC doesn't provide RefSeq accessions
                    'GenBank assembly accession': genome.get('sourceName', 'Not Available'),
                    'Assembly Status': 'Complete' if genome.get('active') else 'N/A'
                }
                data.append(entry)
            
            self.df = pd.DataFrame(data)
            self.logger.info(f"Found {len(self.df)} results from UCSC")
            return self.df
            
        except Exception as e:
            self.logger.error(f"Error searching UCSC: {str(e)}")
            raise

    def _process_ncbi_data(self, summary_results):
        data = []
        for assembly in summary_results['DocumentSummarySet']['DocumentSummary']:
            gb_ftp = assembly.get('FtpPath_GenBank', '')
            rs_ftp = assembly.get('FtpPath_RefSeq', '')
            
            # Get assembly accession and uid
            assembly_accession = assembly.get('AssemblyAccession', 'N/A')
            uid = assembly.attributes['uid']
            
            # Safely get the strain information
            biosource = assembly.get('Biosource', {})
            infraspecies_list = biosource.get('InfraspeciesList', [])
            strain = (infraspecies_list[0].get('Sub_value', 'N/A') 
                     if infraspecies_list 
                     else 'N/A')
            
            entry = {
                'ID': uid,
                'Species Name': assembly.get('SpeciesName', 'N/A'),
                'Strain': strain,
                'Assembly Name': assembly.get('AssemblyName', 'N/A'),
                'RefSeq assembly accession': assembly_accession if rs_ftp else 'Not Available',
                'GenBank assembly accession': assembly_accession if gb_ftp else 'Not Available',
                'Assembly Status': assembly.get('AssemblyStatus', 'N/A')
            }
            
            # Store FTP paths using both UID and assembly accession as keys
            if gb_ftp:
                self.genbank_ftp_dict[uid] = gb_ftp
                self.genbank_ftp_dict[assembly_accession] = gb_ftp
                
            if rs_ftp:
                self.refseq_ftp_dict[uid] = rs_ftp
                self.refseq_ftp_dict[assembly_accession] = rs_ftp
                
            data.append(entry)

        self.df = pd.DataFrame(data)
        self.df.replace('', 'N/A', inplace=True)

        if self.df.empty:
            self.logger.warning("No data found in NCBI response")
        else:
            self.logger.info(f"Processed {len(self.df)} entries from NCBI response")
            self.logger.debug(f"GenBank FTP dict has {len(self.genbank_ftp_dict)} entries")
            self.logger.debug(f"RefSeq FTP dict has {len(self.refseq_ftp_dict)} entries")

    def _get_element_text(self, element, tag):
        el = element.find(tag)
        return el.text if el is not None else 'N/A'

    def get_download_url(self, id, use_genbank):
        """Get download URL based on selected database and parameters"""
        try:
            database = self.current_database
            
            # Set download flags based on search parameters
            self.download_fna = self.current_search_params.get('fna', False)
            self.download_gbff = self.current_search_params.get('gbff', False)
            
            if database == "NCBI GenBank":
                return self._get_ncbi_download_url(id, use_genbank)
            elif database == "ENA (European Nucleotide Archive)":
                return self._get_ena_download_url(id)
            elif database == "UCSC Genome Browser":
                self.logger.warning("Downloads not supported for UCSC Genome Browser")
                return None
            
            return None
        except Exception as e:
            self.logger.error(f"Error getting download URL: {str(e)}")
            return None

    def _get_ena_download_url(self, id):
        """Get download URL for ENA database"""
        try:
            urls = []
            
            # Try to get FASTA (FNA) file
            fasta_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{id}?download=true&gzip=true"
            test_response = requests.head(fasta_url)
            if test_response.status_code == 200:
                self.logger.info(f"Found direct FASTA download URL for {id}")
                urls.append(fasta_url)
                
            # Try to get EMBL (GBFF) file
            embl_url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{id}?download=true&gzip=true"
            test_response = requests.head(embl_url)
            if test_response.status_code == 200:
                self.logger.info(f"Found direct EMBL download URL for {id}")
                urls.append(embl_url)
                
            # If no direct downloads found, try alternative sources
            if not urls:
                # Try WGS format if it's a WGS accession
                if len(id) >= 6 and id.startswith('GCA_'):
                    wgs_id = id.split('_')[1]  # Get the numeric part
                    wgs_fasta_url = f"https://ftp.ebi.ac.uk/pub/databases/ena/wgs/{wgs_id[:6].lower()}/{wgs_id}.fasta.gz"
                    wgs_embl_url = f"https://ftp.ebi.ac.uk/pub/databases/ena/wgs/{wgs_id[:6].lower()}/{wgs_id}.embl.gz"
                    
                    test_response = requests.head(wgs_fasta_url)
                    if test_response.status_code == 200:
                        self.logger.info(f"Found WGS FASTA download URL for {id}")
                        urls.append(wgs_fasta_url)
                        
                    test_response = requests.head(wgs_embl_url)
                    if test_response.status_code == 200:
                        self.logger.info(f"Found WGS EMBL download URL for {id}")
                        urls.append(wgs_embl_url)
                
                # Try assembly FTP path
                assembly_fasta_url = f"https://ftp.ebi.ac.uk/pub/databases/ena/assembly/{id}/{id}.fasta.gz"
                assembly_embl_url = f"https://ftp.ebi.ac.uk/pub/databases/ena/assembly/{id}/{id}.embl.gz"
                
                test_response = requests.head(assembly_fasta_url)
                if test_response.status_code == 200:
                    self.logger.info(f"Found assembly FASTA URL for {id}")
                    urls.append(assembly_fasta_url)
                    
                test_response = requests.head(assembly_embl_url)
                if test_response.status_code == 200:
                    self.logger.info(f"Found assembly EMBL URL for {id}")
                    urls.append(assembly_embl_url)
            
            if urls:
                self.logger.info(f"Found {len(urls)} download URLs for {id}: {urls}")
                return urls
                
            self.logger.warning(f"No suitable download URLs found for accession {id}")
            return None
                
        except Exception as e:
            self.logger.error(f"Error getting ENA download URLs for {id}: {str(e)}")
            return None

    def _get_ncbi_download_url(self, id, use_genbank):
        try:
            # First try with the ID directly
            base_url = None
            if use_genbank:
                base_url = self.genbank_ftp_dict.get(id)
            else:
                base_url = self.refseq_ftp_dict.get(id)
                
            # If not found, try to find the corresponding assembly accession in our DataFrame
            if not base_url and id in self.df.index:
                row = self.df.loc[id]
                accession = row['GenBank assembly accession'] if use_genbank else row['RefSeq assembly accession']
                if accession != 'Not Available':
                    if use_genbank:
                        base_url = self.genbank_ftp_dict.get(accession)
                    else:
                        base_url = self.refseq_ftp_dict.get(accession)

            if not base_url:
                self.logger.warning(f"No FTP path found for ID {id}")
                return None

            # Remove trailing slash if present
            base_url = base_url.rstrip('/')
            
            # Get the assembly accession from the base URL
            assembly_dir = os.path.basename(base_url)
            
            # Return both FNA and GBFF URLs if they exist
            urls = []
            
            # Check FNA file
            fna_url = f"{base_url}/{assembly_dir}_genomic.fna.gz"
            gbff_url = f"{base_url}/{assembly_dir}_genomic.gbff.gz"
            
            # Test if URLs exist
            try:
                parsed_url = urlparse(fna_url)
                ftp = FTP(parsed_url.netloc)
                ftp.login()
                
                # Check FNA file
                try:
                    ftp.size(parsed_url.path[1:])
                    urls.append(fna_url)
                except Exception:
                    # Try alternative FNA name format
                    alt_fna_url = f"{base_url}/{assembly_dir}.fna.gz"
                    try:
                        ftp.size(urlparse(alt_fna_url).path[1:])
                        urls.append(alt_fna_url)
                    except:
                        pass
                
                # Check GBFF file
                try:
                    ftp.size(urlparse(gbff_url).path[1:])
                    urls.append(gbff_url)
                except Exception:
                    # Try alternative GBFF name format
                    alt_gbff_url = f"{base_url}/{assembly_dir}.gbff.gz"
                    try:
                        ftp.size(urlparse(alt_gbff_url).path[1:])
                        urls.append(alt_gbff_url)
                    except:
                        pass
                
                ftp.quit()
                
                return urls if urls else None
                
            except Exception as e:
                self.logger.error(f"Error checking FTP URLs for {id}: {str(e)}")
                return None
                
        except Exception as e:
            self.logger.error(f"Error getting NCBI download URL for {id}: {str(e)}")
            return None

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
        db_path = self.settings.get_db_path()
        self.logger.debug(f"Using database path for downloads: {db_path}")
        return os.path.join(db_path, file_type)

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
