import os
from PyQt6.QtCore import QObject, pyqtSignal, QFileSystemWatcher
import sqlite3
from collections import Counter
import statistics
from enum import Enum
from typing import Set, Dict, List, Tuple
import glob

class FileChangeType(Enum):
    CSPR_ADDED = "cspr_added"
    CSPR_REMOVED = "cspr_removed"
    GBFF_ADDED = "gbff_added"
    GBFF_REMOVED = "gbff_removed"
    OTHER = "other"

class DatabaseManager(QObject):
    db_files_changed = pyqtSignal(dict)  # Emits a dict of FileChangeType: List[str]
    db_validation_changed = pyqtSignal(bool, str)  # Emits validation state and message
    db_state_changed = pyqtSignal(bool, str, dict)  # Emits (is_valid, message, changes)

    def __init__(self, logger, config_manager):
        super().__init__()
        self.logger = logger
        self.config_manager = config_manager
        self.db_path = None
        self.pending_db_path = None
        self.is_changing_directory = False
        self._is_validating = False  # Add flag to prevent recursion
        
        # Initialize last known states before file watcher
        self._last_cspr_files = set()
        self._last_gbff_files = set()
        self._last_files = {}  # Track files in each watched directory
        
        # Initialize file watcher
        self.file_watcher = QFileSystemWatcher()
        self.file_watcher.directoryChanged.connect(self._on_directory_changed)
        
        # Load database path and update states
        self.load_database_path()
        self._update_watched_directory()
        
        # Update last known states after path is loaded
        self._last_cspr_files = set(self._get_cspr_files())
        self._last_gbff_files = set(self._get_gbff_files())
        if self.db_path:
            self._last_files[self.db_path] = set(os.listdir(self.db_path))
            gbff_path = os.path.join(self.db_path, 'GBFF')
            if os.path.exists(gbff_path):
                self._last_files[gbff_path] = set(os.listdir(gbff_path))

    def load_database_path(self):
        """Load the database path from .env file."""
        try:
            db_path = self.config_manager.get_env_value('CSPR_DB', '')
            # Remove both single and double quotes if present
            db_path = db_path.strip("'\"")
            
            # Only set default if no path exists at all
            if not db_path and self.config_manager.get_env_value('FIRST_TIME_START', 'TRUE').upper() == 'TRUE':
                db_path = self.get_default_database_path()
                self.save_db_path(db_path)
            
            self.db_path = db_path
            self.logger.debug(f"Database path loaded from .env: {self.db_path}")
            return self.db_path
            
        except Exception as e:
            self.logger.error(f"Error loading database path: {str(e)}")
            return self.get_default_database_path()

    def validate_db_path(self, path):
        """
        Validate the database path without modifying it
        Returns (is_valid, message)
        """
        try:
            if not path:
                return False, "No directory selected"
            
            if not os.path.exists(path):
                return False, "The selected directory does not exist."
            
            # Check for CSPR files
            cspr_files = glob.glob(os.path.join(path, "*.cspr"))
            if not cspr_files:
                return False, "No CSPR files found"
            
            return True, "Valid database directory"
            
        except Exception as e:
            self.logger.error(f"Error validating database path: {str(e)}")
            return False, str(e)

    def save_db_path(self, path):
        """Set and save the database path."""
        if self._is_validating:  # Prevent recursive validation
            return True, "Operation in progress"
        
        try:
            self._is_validating = True
            
            if not path:
                self.logger.warning("Attempting to save an empty database path")
                return False, "Empty database path is not allowed."

            path = str(path).strip("'\"")
            is_new_genome = self._is_new_genome_context()
            
            # Validate the database path
            is_valid, message = self.validate_db_path(path)
            
            # Log the current state
            self.logger.debug(f"Saving DB path - Current state: db_path={self.db_path}, pending={self.pending_db_path}, "
                             f"is_changing={self.is_changing_directory}, is_new_genome={is_new_genome}")
            
            # If we're in new genome context or path change context, store as pending path
            if is_new_genome or is_valid:
                self.logger.debug(f"Storing pending database path: {path}")
                self.pending_db_path = path
                self.is_changing_directory = True
                
                # Create directory if needed
                if not os.path.exists(path):
                    try:
                        os.makedirs(path)
                        self.logger.info(f"Created pending directory: {path}")
                    except Exception as e:
                        self.logger.error(f"Error creating pending directory: {str(e)}")
                
                # If the path is valid, finalize the change immediately
                if is_valid:
                    success, finalize_message = self.finalize_directory_change()
                    if not success:
                        self.logger.error(f"Failed to finalize directory change: {finalize_message}")
                        return False, finalize_message
                    return True, finalize_message
                
                return True, "Path stored for new genome creation"
            
            # For invalid paths
            if not is_valid:
                self.db_validation_changed.emit(False, message)
                return False, message
            
        finally:
            self._is_validating = False

    def _is_new_genome_context(self):
        """Check if the path change is happening in new genome context"""
        try:
            import inspect
            stack = inspect.stack()
            
            # Check for NCBI context as well
            is_new_genome = any('NewGenome' in frame.filename or 'NCBI' in frame.filename for frame in stack)
            is_path_change = any('MainWindow' in frame.filename and 'change_database_directory' in frame.function 
                               for frame in stack)
            
            # Consider it a new genome context if either:
            # 1. We're in new genome/NCBI context and actively changing directory
            # 2. We're in the path change process
            is_active_change = (is_new_genome and self.is_changing_directory) or is_path_change
            
            self.logger.debug(f"Context check - New Genome: {is_new_genome}, Path Change: {is_path_change}, "
                             f"Active Change: {is_active_change}, Is Changing Directory: {self.is_changing_directory}")
            
            return is_active_change
            
        except Exception as e:
            self.logger.error(f"Error in _is_new_genome_context: {str(e)}")
            return False

    def get_active_db_path(self):
        """Get the appropriate database path based on context"""
        if self.pending_db_path and self.is_changing_directory:
            self.logger.debug(f"Using pending database path: {self.pending_db_path}")
            return self.pending_db_path
        
        self.logger.debug(f"Using current database path: {self.db_path}")
        return self.db_path

    def get_db_path(self):
        return self.db_path

    def ensure_db_path_exists(self):
        """Ensure that the database path exists, creating it if necessary."""
        path_to_check = self.get_active_db_path()  # Use active path instead of db_path
        if not os.path.exists(path_to_check):
            try:
                os.makedirs(path_to_check)
                self.logger.info(f"Created database directory: {path_to_check}")
            except Exception as e:
                self.logger.error(f"Failed to create database directory: {path_to_check}. Error: {str(e)}")
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
        """Update the watched directory and validate the new path"""
        try:
            # Remove old watched directories
            self.file_watcher.removePaths(self.file_watcher.directories())
            
            if self.db_path and os.path.isdir(self.db_path):
                # Add new directory to watch
                self.file_watcher.addPath(self.db_path)
                
                # Add GBFF subdirectory if it exists
                gbff_path = os.path.join(self.db_path, 'GBFF')
                if os.path.isdir(gbff_path):
                    self.file_watcher.addPath(gbff_path)
                
                self.logger.debug(f"Now watching directories: {self.file_watcher.directories()}")
                
                # Validate the new path and emit signals
                is_valid, message = self.validate_db_path(self.db_path)
                self.db_validation_changed.emit(is_valid, message)
                
                # Also check for any file changes
                changes = self._detect_file_changes()
                if changes:
                    self.db_files_changed.emit(changes)
                    self.db_state_changed.emit(is_valid, message, changes)
                
        except Exception as e:
            self.logger.error(f"Error updating watched directory: {str(e)}")

    def _detect_file_changes(self) -> Dict[FileChangeType, List[str]]:
        """Detect what files have changed and categorize the changes"""
        current_cspr_files = set(self._get_cspr_files())
        current_gbff_files = set(self._get_gbff_files())
        
        changes = {}
        
        # Detect CSPR changes
        cspr_added = current_cspr_files - self._last_cspr_files
        cspr_removed = self._last_cspr_files - current_cspr_files
        
        if cspr_added:
            changes[FileChangeType.CSPR_ADDED] = list(cspr_added)
        if cspr_removed:
            changes[FileChangeType.CSPR_REMOVED] = list(cspr_removed)
            
        # Detect GBFF changes
        gbff_added = current_gbff_files - self._last_gbff_files
        gbff_removed = self._last_gbff_files - current_gbff_files
        
        if gbff_added:
            changes[FileChangeType.GBFF_ADDED] = list(gbff_added)
        if gbff_removed:
            changes[FileChangeType.GBFF_REMOVED] = list(gbff_removed)
            
        # Update last known state
        self._last_cspr_files = current_cspr_files
        self._last_gbff_files = current_gbff_files
        
        return changes

    def _on_directory_changed(self, path):
        """Handle changes in the watched directory"""
        try:
            self.logger.debug(f"Detected change in directory: {path}")
            
            # Check if change is just an index file
            changed_files = set(os.listdir(path)) - set(self._last_files.get(path, []))
            if all(f.endswith('.index') for f in changed_files):
                self.logger.debug("Ignoring index file changes")
                # Update last files without triggering refresh
                self._last_files[path] = set(os.listdir(path))
                return
            
            # Get current state of CSPR files
            current_cspr_files = set(self._get_cspr_files())
            
            # Re-validate the path
            is_valid, message = self.validate_db_path(self.db_path)
            
            # Detect specific changes
            changes = self._detect_file_changes()
            
            # Update last known state
            self._last_cspr_files = current_cspr_files
            self._last_files[path] = set(os.listdir(path))
            
            # Always emit validation and state changes
            self.db_validation_changed.emit(is_valid, message)
            self.db_state_changed.emit(is_valid, message, changes)
            
            # If changes detected, emit files changed signal
            if changes:
                self.logger.debug(f"Detected file changes: {changes}")
                self.db_files_changed.emit(changes)
            
        except Exception as e:
            self.logger.error(f"Error handling directory change: {str(e)}")

    def _get_cspr_files(self):
        """Get a list of CSPR files in the current database directory"""
        if not self.db_path or not os.path.isdir(self.db_path):
            return []
        return [f for f in os.listdir(self.db_path) 
                if f.endswith('.cspr')]

    def _get_gbff_files(self):
        """Get a list of GBFF files in the database directory"""
        if not self.db_path or not os.path.isdir(self.db_path):
            return []
        gbff_path = os.path.join(self.db_path, 'GBFF')
        if not os.path.exists(gbff_path):
            return []
        return [f for f in os.listdir(gbff_path) 
                if f.endswith('.gbff')]

    def check_db_state(self):
        """
        Check database state without clearing invalid paths
        """
        try:
            current_path = self.get_db_path()
            is_valid, message = self.validate_db_path(current_path)
            
            # Get list of changes if path is valid
            changes = {}  # Initialize as dict instead of list
            if is_valid:
                changes = self._detect_file_changes()
            
            # Emit signals but don't modify the path
            self.db_validation_changed.emit(is_valid, message)
            self.db_state_changed.emit(is_valid, message, changes)
            
        except Exception as e:
            self.logger.error(f"Error checking database state: {str(e)}")

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
                    newname = file[0:-5] 
                    endonuclease = newname[newname.rfind("_") + 1:] 
                    
                    # Read organism name from first line of CSPR file
                    file_path = os.path.join(self.db_path, file)
                    with open(file_path, 'r') as hold:
                        buf = hold.readline().strip()
                        species = buf.replace("GENOME: ", "")
                    
                    # Store file mappings
                    if species in organisms_to_files:
                        organisms_to_files[species][endonuclease] = [file, file.replace(".cspr", "_repeats.db")]
                    else:
                        organisms_to_files[species] = {}
                        organisms_to_files[species][endonuclease] = [file, file.replace(".cspr", "_repeats.db")]
                    
                    # Store endonuclease mappings
                    if species in organisms_to_endos:
                        if endonuclease not in organisms_to_endos[species]:
                            organisms_to_endos[species].append(endonuclease)
                    else:
                        organisms_to_endos[species] = [endonuclease]
                        
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

    def update_db_state(self):
        """Check and update the database state"""
        self.logger.debug("Checking database state")
        if not self.db_path and not self.pending_db_path:
            self.load_database_path()

        # Use active path for validation
        path_to_check = self.get_active_db_path()
        self.logger.debug(f"Checking state for path: {path_to_check} (pending: {self.pending_db_path}, current: {self.db_path})")
        
        is_valid, message = self.validate_db_path(path_to_check)
        self.logger.debug(f"Database validation result - Path: {path_to_check}, Valid: {is_valid}, Message: {message}")
        
        # Detect any changes since last check
        changes = self._detect_file_changes()
        
        # Check if we should finalize a directory change
        if self.pending_db_path and self.is_changing_directory:
            self.logger.debug("Checking conditions for directory change finalization")
            self.logger.debug(f"Is valid: {is_valid}, Has changes: {bool(changes)}")
            
            if is_valid and not changes:  # No changes means we're not in the middle of file operations
                self.logger.debug("Attempting to finalize directory change")
                success, finalize_message = self.finalize_directory_change()
                if success:
                    self.logger.debug("Directory change finalized successfully")
                    # Emit signals
                    self.db_validation_changed.emit(True, finalize_message)
                    self.db_state_changed.emit(True, finalize_message, changes)
                    return
                else:
                    self.logger.debug(f"Directory change finalization failed: {finalize_message}")
        
        # Emit regular signals
        self.db_validation_changed.emit(is_valid, message)
        if changes:
            self.db_files_changed.emit(changes)
            
        self.logger.info(f"Database state checked - Valid: {is_valid}, Changes: {changes}")

    def cancel_directory_change(self):
        """Cancel the directory change process"""
        self.logger.debug(f"Cancelling directory change process. Previous state - Pending: {self.pending_db_path}, Changing: {self.is_changing_directory}")
        self.pending_db_path = None
        self.is_changing_directory = False
        self.logger.debug("Directory change process cancelled")

    def finalize_directory_change(self):
        """Finalize the database directory change after successful validation"""
        try:
            if self.pending_db_path and self.is_changing_directory:
                self.logger.debug(f"Finalizing directory change from {self.db_path} to {self.pending_db_path}")
                
                # Validate the pending path one final time
                is_valid, message = self.validate_db_path(self.pending_db_path)
                if not is_valid:
                    self.logger.warning(f"Cannot finalize directory change: {message}")
                    return False, message
                
                # Update the current path
                old_path = self.db_path
                self.db_path = self.pending_db_path
                
                # Update .env file - try multiple approaches to ensure it works
                try:
                    # First attempt: Direct write
                    self.config_manager.write_to_env('CSPR_DB', self.db_path)
                    
                    # Second attempt: Use set_env_value
                    self.config_manager.set_env_value('CSPR_DB', self.db_path)
                    
                    # Force reload environment variables
                    self.config_manager.load_env()
                    
                    # Verify the update
                    new_env_value = self.config_manager.get_env_value('CSPR_DB')
                    if new_env_value != self.db_path:
                        raise Exception(f"Environment variable update failed. Expected: {self.db_path}, Got: {new_env_value}")
                        
                except Exception as e:
                    self.logger.error(f"Error updating environment variable: {str(e)}")
                    return False, f"Failed to update environment variable: {str(e)}"
                
                # Clear pending state
                self.pending_db_path = None
                self.is_changing_directory = False
                
                # Update database state
                self.update_db_state()
                
                success_message = f"Successfully changed database directory to:\n{self.db_path}"
                self.logger.info(f"Successfully changed database directory from {old_path} to {self.db_path}")
                self.logger.debug("Directory change finalized successfully")
                
                return True, success_message
                
            return False, "No pending directory change to finalize"
            
        except Exception as e:
            self.logger.error(f"Error finalizing directory change: {str(e)}")
            return False, f"Error finalizing directory change: {str(e)}"
