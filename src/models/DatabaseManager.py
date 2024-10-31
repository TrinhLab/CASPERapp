import os
from PyQt6.QtCore import QObject, pyqtSignal, QFileSystemWatcher
import sqlite3
from collections import Counter
import statistics

class DatabaseManager(QObject):
    db_state_updated = pyqtSignal(bool, str, list)  # Combined signal

    def __init__(self, logger, config_manager):
        super().__init__()
        self.logger = logger
        self.config_manager = config_manager
        self.db_path = None
        self.file_watcher = QFileSystemWatcher()  # Initialize file_watcher here
        self.file_watcher.directoryChanged.connect(self._on_directory_changed)
        self.load_database_path()
        self._update_watched_directory()

    def load_database_path(self):
        """Load the database path from .env file or set default if empty."""
        db_path = self.config_manager.get_env_value('CSPR_DB', '')
        # Remove both single and double quotes if present
        db_path = db_path.strip("'\"")
        if not db_path:
            db_path = self.get_default_database_path()
            self.save_db_path(db_path)
        self.db_path = db_path
        return self.db_path

    def validate_db_path(self, path):
        """Validate that the given path exists and contains CSPR files."""
        self.logger.debug(f"Validating DB path: {path}")
        if not os.path.isdir(path):
            self.logger.debug(f"Path is not a directory: {path}")
            return False, "The selected path is not a directory."
        has_cspr_files = any(file.endswith(".cspr") for file in os.listdir(path))
        if not has_cspr_files:
            self.logger.debug(f"Path {path} does not contain CSPR files")
            return False, "The selected directory does not contain any CSPR files."
        self.logger.debug(f"Path {path} is valid and contains CSPR files")
        return True, "Valid database path selected."

    def save_db_path(self, path):
        """Set and save the database path."""
        if not path:
            self.logger.warning("Attempting to save an empty database path")
            return False, "Empty database path is not allowed."

        # Ensure the path is a string and properly quoted
        path = str(path).strip("'\"")

        # Validate the database path
        is_valid, message = self.validate_db_path(path)
        if not is_valid:
            self.logger.warning(f"Invalid database path: {path}")
            self.db_state_updated.emit(False, message, [])
            self.db_path = path
            self.config_manager.set_env_value('CSPR_DB', path)
            self._update_watched_directory()
            return False, message

        # Set the db_path attribute
        self.db_path = path

        try:
            self.config_manager.set_env_value('CSPR_DB', path)
            self.logger.info(f"Database path set and saved: {path}")
            self.db_state_updated.emit(True, "Database path saved successfully.", [])
            self._update_watched_directory()
            return True, "Database path saved successfully."
        except Exception as e:
            error_message = f"Error saving database path: {str(e)}"
            self.logger.error(error_message)
            self.db_state_updated.emit(False, error_message, [])
            return False, error_message

    def get_db_path(self):
        return self.db_path

    def ensure_db_path_exists(self):
        """Ensure that the database path exists, creating it if necessary."""
        if not os.path.exists(self.db_path):
            try:
                os.makedirs(self.db_path)
                self.logger.info(f"Created database directory: {self.db_path}")
            except Exception as e:
                self.logger.error(f"Failed to create database directory: {self.db_path}. Error: {str(e)}")
                raise

    def get_default_database_path(self):
        """Get the default database path."""
        documents_path = os.path.expanduser("~/Documents")
        default_path = os.path.join(documents_path, 'CASPERdb')
        return self.adjust_path_for_os(default_path)

    def adjust_path_for_os(self, path):
        """Adjust the file path based on the operating system."""
        if os.name == "nt":  # Windows
            adjusted_path = path.replace("/", "\\")
        else:  # Unix-like systems
            adjusted_path = path.replace("\\", "/")
        
        # Add trailing slash for directories on non-Windows systems
        # if os.name != "nt" and os.path.isdir(adjusted_path):
        #     adjusted_path = os.path.join(adjusted_path, '')
        
        return adjusted_path

    def _update_watched_directory(self):
        """Update the directory being watched by QFileSystemWatcher."""
        if self.file_watcher.directories():
            self.file_watcher.removePaths(self.file_watcher.directories())
        if self.db_path and os.path.isdir(self.db_path):
            self.file_watcher.addPath(self.db_path)
            self.logger.debug(f"Now watching directory: {self.db_path}")

    def _on_directory_changed(self, path):
        """Handle changes in the watched directory."""
        self.logger.debug(f"Detected change in directory: {path}")
        
        # Get current state
        is_valid, message = self.validate_db_path(path)
        
        # Get list of files
        cspr_files = self._get_cspr_files()
        gbff_files = self._get_gbff_files()  # Add method to get GBFF files
        
        # Emit the signal with updated state
        self.db_state_updated.emit(is_valid, message, cspr_files)
        
        # Log the change
        self.logger.info(f"Database state updated - Valid: {is_valid}, Files: {len(cspr_files)} CSPR, {len(gbff_files)} GBFF")

    def _get_cspr_files(self):
        """Get a list of CSPR files in the current database directory."""
        if not self.db_path or not os.path.isdir(self.db_path):
            return []
        return [f for f in os.listdir(self.db_path) if f.endswith('.cspr')]

    def _get_gbff_files(self):
        """Get a list of GBFF files in the database directory."""
        if not self.db_path or not os.path.isdir(self.db_path):
            return []
        gbff_path = os.path.join(self.db_path, 'GBFF')
        if not os.path.exists(gbff_path):
            return []
        return [f for f in os.listdir(gbff_path) if f.endswith('.gbff')]

    def check_db_state(self):
        """Check the current state of the database and emit signals if changed."""
        self.logger.debug("Checking database state")
        if not self.db_path:
            self.load_database_path()

        is_valid, message = self.validate_db_path(self.db_path)
        self.logger.debug(f"Database state: valid={is_valid}, message={message}")
        
        cspr_files = self._get_cspr_files()
        gbff_files = self._get_gbff_files()
        
        message = f"Database is valid. Contains {len(cspr_files)} CSPR files and {len(gbff_files)} GBFF files."
        
        self.db_state_updated.emit(is_valid, message, cspr_files)
        self.logger.info(f"Database state checked - Valid: {is_valid}, Files: {len(cspr_files)} CSPR, {len(gbff_files)} GBFF")

    def get_organisms_and_endos(self):
        """Get mapping of organisms to their endonucleases and files"""
        try:
            if not self.db_path or not os.path.exists(self.db_path):
                self.logger.error(f"Invalid database path: {self.db_path}")
                return {}, {}
                
            onlyfiles = [f for f in os.listdir(self.db_path) if os.path.isfile(os.path.join(self.db_path, f))]
            cspr_files = [f for f in onlyfiles if f.endswith('.cspr')]
            
            organisms_to_files = {}
            organisms_to_endos = {}
            
            for file in cspr_files:
                try:
                    # Parse filename
                    newname = file[0:-5]  # Remove .cspr - changed from -4 to -5
                    endo = newname[newname.rfind("_") + 1:]  # Get endonuclease name
                    
                    # Read organism name from first line of CSPR file
                    file_path = os.path.join(self.db_path, file)
                    with open(file_path, 'r') as hold:
                        buf = hold.readline().strip()
                        species = buf.replace("GENOME: ", "")
                    
                    # Store file mappings
                    if species in organisms_to_files:
                        organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]
                    else:
                        organisms_to_files[species] = {}
                        organisms_to_files[species][endo] = [file, file.replace(".cspr", "_repeats.db")]
                    
                    # Store endonuclease mappings
                    if species in organisms_to_endos:
                        if endo not in organisms_to_endos[species]:
                            organisms_to_endos[species].append(endo)
                    else:
                        organisms_to_endos[species] = [endo]
                        
                except Exception as e:
                    self.logger.error(f"Error processing file {file}: {str(e)}")
                    continue
                    
            return organisms_to_files, organisms_to_endos
            
        except Exception as e:
            self.logger.error(f"Error getting organisms and endonucleases: {str(e)}")
            return {}, {}

    def get_repeats_data(self, db_file, row_limit=-1):
        """Get repeats data for the seeds table"""
        try:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            
            if row_limit == -1:
                sql_query = "SELECT * FROM repeats ORDER BY count DESC;"
            else:
                sql_query = f"SELECT * FROM repeats ORDER BY count DESC LIMIT 0, {row_limit};"

            repeats = c.execute(sql_query).fetchall()
            processed_data = []
            
            for repeat in repeats:
                # Extract repeat info
                seed = repeat[0]
                chroms = repeat[1].split(",")
                locs = repeat[2].split(",")
                threes = repeat[3].split(",")
                fives = repeat[4].split(",")
                pams = repeat[5].split(",")
                scores = repeat[6].split(",")
                count = repeat[7]

                # Handle missing data in threes/fives
                if len(threes) < len(fives):
                    threes.extend([''] * (len(fives) - len(threes)))
                elif len(fives) < len(threes):
                    fives.extend([''] * (len(threes) - len(fives)))

                # Find majority sequence
                majority_index = 0
                three_prime, five_prime, both_prime = False, False, False
                if threes[0] == '':
                    majority = max(set(fives), key=fives.count)
                    majority_index = fives.index(majority)
                    five_prime = True
                elif fives[0] == '':
                    majority = max(set(threes), key=threes.count)
                    majority_index = threes.index(majority)
                    three_prime = True
                else:
                    # account for both 3 and 5 present
                    threes_and_fives = []
                    for i in range(len(threes)):
                        threes_and_fives.append(threes[i] + fives[i])
                    majority = max(set(threes_and_fives), key=threes_and_fives.count)
                    majority_index = threes_and_fives.index(majority)
                    both_prime = True

                # Calculate average repeats per scaffold
                location_repeat_counts = Counter(chroms)
                avg_rep_per_scaff = sum(location_repeat_counts.values()) / len(location_repeat_counts.values())
                avg_rep_per_scaff = float("%.2f" % avg_rep_per_scaff)

                # Calculate consensus percentage
                if five_prime:
                    percent_consensus = (fives.count(fives[majority_index]) / len(fives)) * 100
                elif three_prime:
                    percent_consensus = (threes.count(threes[majority_index]) / len(threes)) * 100
                elif both_prime:
                    percent_consensus = (threes_and_fives.count(threes_and_fives[majority_index]) / len(threes_and_fives)) * 100
                percent_consensus = float("%.2f" % percent_consensus)

                # Determine strand
                strand = "+" if int(locs[majority_index]) >= 0 else "-"

                # Create processed row
                processed_row = (
                    seed,                                           # Seed
                    count,                                         # Total Repeats
                    avg_rep_per_scaff,                            # Avg. Repeats/Scaffold
                    fives[majority_index] + seed + threes[majority_index],  # Consensus Sequence
                    percent_consensus,                             # % Consensus
                    scores[majority_index],                        # Score
                    pams[majority_index],                         # PAM
                    strand                                        # Strand
                )
                processed_data.append(processed_row)

            conn.close()
            return processed_data
            
        except Exception as e:
            self.logger.error(f"Error getting repeats data: {str(e)}")
            raise

    def get_seed_data(self, db_file, seed):
        """Get detailed data for a specific seed"""
        try:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            data = c.execute("""
                SELECT chromosome, location, pam, score, 
                       five, three
                FROM repeats 
                WHERE seed = ?
            """, (seed,)).fetchall()
            conn.close()
            return data
        except Exception as e:
            self.logger.error(f"Error getting seed data: {str(e)}")
            raise

    def get_chro_bar_data(self, db_file, seed):
        """Get chromosome distribution data for a seed"""
        try:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            data = c.execute("SELECT chromosome FROM repeats WHERE seed = ?", (seed,)).fetchone()
            conn.close()
            
            if not data:
                return Counter()
                
            # Split the chromosome string and convert to integers
            chromosomes = [int(x) for x in data[0].split(',')]
            counts = Counter(chromosomes)
            return counts
            
        except Exception as e:
            self.logger.error(f"Error getting chromosome bar data: {str(e)}")
            raise

    def get_seeds_vs_repeats_data(self, db_file):
        """Get data for seeds vs repeats plot"""
        try:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            
            # Count how many sequences have each repeat count
            data = c.execute("""
                SELECT count, COUNT(*) as num_sequences
                FROM repeats 
                GROUP BY count 
                ORDER BY count ASC
            """).fetchall()
            
            conn.close()
            
            if not data:
                return None
                
            # Separate into x and y values
            x_vals = [row[0] for row in data]  # number of repeats
            y_vals = [row[1] for row in data]  # number of sequences with that count
            
            return {
                'x_vals': x_vals,
                'y_vals': y_vals
            }
            
        except Exception as e:
            self.logger.error(f"Error getting seeds vs repeats data: {str(e)}")
            raise

    def get_repeats_vs_seeds_data(self, db_file):
        """Get data for repeats vs seeds plot"""
        try:
            self.logger.debug(f"Getting repeats vs seeds data from {db_file}")
            
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            
            # Get all count values ordered by rowid to maintain order
            self.logger.debug("Executing SQL query")
            data = c.execute("SELECT count FROM repeats ORDER BY rowid;").fetchall()
            self.logger.debug(f"Raw data from database: {data[:5]}...")  # Log first 5 entries
            
            conn.close()
            
            # Extract counts from tuples
            counts = [row[0] for row in data]
            self.logger.debug(f"Processed counts (first 5): {counts[:5]}...")
            
            if not counts:
                self.logger.warning("No count data found")
                return None
                
            # Calculate statistics
            stats = {
                'average': statistics.mean(counts),
                'median': statistics.median(counts),
                'mode': statistics.mode(counts),
                'repeat_count': len(counts)
            }
            self.logger.debug(f"Calculated statistics: {stats}")
            
            result = {
                'counts': counts,
                'stats': stats
            }
            self.logger.debug("Successfully prepared repeats vs seeds data")
            return result
            
        except Exception as e:
            self.logger.error(f"Error getting repeats vs seeds data: {str(e)}")
            raise

    def calculate_statistics(self, db_file):
        """Calculate global statistics"""
        try:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            
            stats = {}
            # Add your statistics calculations here
            
            conn.close()
            return stats
        except Exception as e:
            self.logger.error(f"Error calculating statistics: {str(e)}")
            raise
