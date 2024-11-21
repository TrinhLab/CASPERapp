from math import e
import os
import platform
from PyQt6.QtCore import QProcess, pyqtSignal, QObject
import logging
import traceback

class OffTargetModel(QObject):
    # Update signal to emit tuple of (scores, details)
    results_ready = pyqtSignal(tuple)  # Emits (scores_dict, details_dict)
    progress_updated = pyqtSignal(int, str)  # Emits (progress_value, status_message)

    def __init__(self, global_settings):
        super().__init__()  # Make sure OffTargetModel inherits from QObject
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.process = None
        self.organisms_to_files = {}
        self.organisms_to_endos = {}
        self._current_parameters = None
        self._load_organisms_and_endos()

    def _load_organisms_and_endos(self):
        """Load available organisms and endonucleases from CSPR files"""
        try:
            cspr_files = [f for f in os.listdir(self.global_settings.get_db_path()) 
                         if f.endswith('.cspr')]
            
            for file in cspr_files:
                newname = file[:-4]
                endo = newname[newname.rfind("_") + 1:-1]
                
                with open(os.path.join(self.global_settings.get_db_path(), file), 'r') as f:
                    species = f.readline().strip().replace("GENOME: ", "")
                
                # Store file mappings
                if species not in self.organisms_to_files:
                    self.organisms_to_files[species] = {}
                self.organisms_to_files[species][endo] = [
                    file,
                    file.replace(".cspr", "_repeats.db")
                ]
                
                # Store endo mappings
                if species not in self.organisms_to_endos:
                    self.organisms_to_endos[species] = []
                if endo not in self.organisms_to_endos[species]:
                    self.organisms_to_endos[species].append(endo)
                    
        except Exception as e:
            self.logger.error(f"Error loading organisms and endonucleases: {str(e)}")
            raise

    def get_organisms(self):
        """Get list of available organisms"""
        return sorted(self.organisms_to_files.keys())

    def get_endonucleases(self, organism):
        """Get available endonucleases for given organism"""
        return sorted(self.organisms_to_endos.get(organism, []))

    def get_file_paths(self, organism, endonuclease):
        """Get CSPR and DB file paths for given organism and endonuclease"""
        try:
            files = self.organisms_to_files[organism][endonuclease]
            return {
                'cspr_file': os.path.join(self.global_settings.get_db_path(), files[0]),
                'db_file': os.path.join(self.global_settings.get_db_path(), files[1])
            }
        except KeyError:
            self.logger.error(f"No files found for {organism} and {endonuclease}")
            return None
        
    def get_program_command(self):
        if platform.system() == 'Windows':
            return 'OT_Win.exe'
        else:
            return 'OT_Lin' if platform.system() == 'Linux' else 'OT_Mac'

    def start_analysis(self, parameters):
        """Start off-target analysis with given parameters"""
        try:
            self._current_parameters = parameters
            # Create QProcess if not exists
            if not self.process:
                self.process = QProcess()
                # Connect process signals for logging
                self.process.readyReadStandardOutput.connect(self._handle_process_output)
                self.process.readyReadStandardError.connect(self._handle_process_error)
                self.process.finished.connect(self._handle_process_finished)
            
            # Get file paths
            files = self.get_file_paths(parameters['organism'], parameters['endonuclease'])
            if not files:
                raise ValueError("Could not get required file paths")
                
            # Build command
            program_path, cmd = self._build_command(parameters, files)
            
            # Log analysis start and parameters
            self.logger.debug("=== Starting off-target analysis ===")
            self.logger.debug(f"Organism: {parameters['organism']}")
            self.logger.debug(f"Endonuclease: {parameters['endonuclease']}")
            self.logger.debug(f"Number of targets: {len(parameters.get('targets', []))}")
            self.logger.debug(f"Max mismatches: {parameters['max_mismatches']}")
            self.logger.debug(f"Tolerance: {parameters['tolerance']}")
            self.logger.debug(f"Average output: {parameters['average_output']}")
            
            # Set working directory
            off_target_dir = self.global_settings.get_off_target_dir_path()
            self.process.setWorkingDirectory(off_target_dir)
            
            # Set executable permissions on Unix systems
            if platform.system() != 'Windows':
                import stat
                st = os.stat(program_path)
                os.chmod(program_path, st.st_mode | stat.S_IEXEC)
                
            # Verify temp file exists and has content
            temp_path = os.path.join(off_target_dir, 'temp.txt')
            if os.path.exists(temp_path):
                with open(temp_path, 'r') as f:
                    content = f.read()
                    self.logger.debug(f"Temp file contents:\n{content}")
            else:
                self.logger.error("Temp file does not exist before starting process")

            # Start process
            self.logger.debug("Starting QProcess with command:")
            self.logger.debug(cmd)


            example_cmd = [
                "/Users/admin/Documents/proj/CASPERtest/CASPERapp/src/models/OffTarget/temp.txt",
                "spCas9",
                "/Users/admin/Documents/CASPERdb2/eck_12_spCas9.cspr",
                "/Users/admin/Documents/CASPERdb2/eck_12_spCas9_repeats.db",
                "/Users/admin/Documents/CASPERdb2/testtttt",
                "/Users/admin/Documents/proj/CASPERtest/CASPERapp/config/CASPERinfo",
                "4",
                "0.05",
                "FALSE",
                "TRUE",
                "MATRIX:HSU MATRIX-spCas9-2013"
            ]

            print(f"cmd: {cmd}")

            print(f"example_cmd: {example_cmd}")


            self.process.start(str(program_path), cmd)
            

                
            return True
            
        except Exception as e:
            self.logger.error(f"Error starting analysis: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            raise

    def _handle_process_output(self):
        """Handle standard output from process"""
        try:
            output = self.process.readAllStandardOutput().data().decode()
            self.logger.debug(f"Process output: {output}")
            
            # Update progress based on output content
            if "Reading in Data" in output:
                self.progress_updated.emit(10, "Reading data...")
            elif "Loading data for algorithm" in output:
                self.progress_updated.emit(25, "Loading algorithm data...")
            elif "Running OffTarget Analysis" in output:
                self.progress_updated.emit(50, "Running analysis...")
            elif "Writing Output" in output:
                self.progress_updated.emit(75, "Writing results...")
            elif "Analysis Complete" in output:
                self.progress_updated.emit(90, "Finalizing...")
                
        except Exception as e:
            self.logger.error(f"Error handling process output: {str(e)}")

    def _handle_process_error(self):
        """Handle standard error from process"""
        try:
            error = self.process.readAllStandardError().data().decode()
            self.logger.error(f"Process error output: {error}")
        except Exception as e:
            self.logger.error(f"Error handling process error output: {str(e)}")

    def _handle_process_finished(self, exit_code, exit_status):
        """Handle process completion and parse results"""
        try:
            self.logger.debug(f"Process finished with exit code: {exit_code}")
            
            if exit_code == 0:  # Success
                # Get output file path
                output_path = self._get_output_path(self._current_parameters)
                
                if os.path.exists(output_path):
                    # Parse results
                    results, detailed_results = self._parse_off_target_results(output_path)
                    if results is not None:
                        # Update progress to 100%
                        self.progress_updated.emit(100, "Analysis complete")
                        # Emit results signal with both scores and details
                        self.results_ready.emit((results, detailed_results))
                        self.logger.debug(f"Parsed and emitted {len(results)} off-target results with {len(detailed_results)} detailed results")
                    else:
                        self.progress_updated.emit(0, "Analysis failed - no results")
                        self.results_ready.emit(({}, {}))  # Emit empty results
                else:
                    self.logger.error(f"Output file not found: {output_path}")
                    self.progress_updated.emit(0, "Analysis failed - no output file")
                    self.results_ready.emit(({}, {}))  # Emit empty results
            else:
                self.logger.error(f"Process failed with exit code: {exit_code}")
                self.progress_updated.emit(0, f"Analysis failed with code {exit_code}")
                self.results_ready.emit(({}, {}))  # Emit empty results
                
            # Always cleanup temp file after process finishes
            self._cleanup_temp_file()
                
        except Exception as e:
            self.logger.error(f"Error in process finished handler: {str(e)}")
            self.progress_updated.emit(0, f"Error: {str(e)}")
            self.results_ready.emit(({}, {}))  # Emit empty results
            # Ensure cleanup happens even if there's an error
            self._cleanup_temp_file()

    def _cleanup_temp_file(self):
        """Clean up temporary input file"""
        try:
            # Get path to temp file
            db_path = self.global_settings.get_db_path()
            temp_path = os.path.join(db_path, 'temp.txt')
            
            # Remove temp file if it exists
            if os.path.exists(temp_path):
                os.remove(temp_path)
                self.logger.debug(f"Removed temp file: {temp_path}")
            
        except Exception as e:
            self.logger.error(f"Error cleaning up temp file: {str(e)}")

    def _parse_off_target_results(self, file_path):
        """Parse off-target analysis results file"""
        try:
            results = {}
            detailed_results = {}
            current_sequence = None
            details = []
            
            with open(file_path, 'r') as f:
                lines = f.readlines()
                
            output_type = lines[0].strip()
            
            # Process based on output type
            if output_type == "DETAILED OUTPUT":
                for line in lines[1:]:  # Skip header
                    line = line.strip()
                    if not line:
                        continue
                        
                    if ':' in line:  # This is a sequence line
                        # Save previous sequence details if they exist
                        if current_sequence and details:
                            detailed_results[current_sequence] = details
                            
                        # Start new sequence
                        parts = line.split(':')
                        current_sequence = parts[0].strip()
                        try:
                            score = float(parts[1].strip())
                            results[current_sequence] = score
                        except ValueError:
                            self.logger.warning(f"Invalid score format in line: {line}")
                        details = []
                    else:  # This is a detail line
                        details.append(line)
                        
                # Don't forget the last sequence
                if current_sequence and details:
                    detailed_results[current_sequence] = details
                    
            elif output_type == "AVG OUTPUT":
                for line in lines[1:]:  # Skip header
                    line = line.strip()
                    if not line:
                        continue
                        
                    parts = line.split(':')
                    if len(parts) == 2:
                        sequence = parts[0].strip()
                        try:
                            score = float(parts[1].strip())
                            results[sequence] = score
                        except ValueError:
                            self.logger.warning(f"Invalid score format in line: {line}")
                            
            # Return both results and detailed_results
            return results, detailed_results
            
        except Exception as e:
            self.logger.error(f"Error parsing off-target results: {str(e)}")
            return None, None

    def _build_command(self, parameters, files):
        """Build command string for off-target analysis"""
        try:
            off_target_dir = self.global_settings.get_off_target_dir_path()
            
            # Get executable path based on platform
            if platform.system() == 'Windows':
                exe_name = 'OT_Win.exe'
            elif platform.system() == 'Linux':
                exe_name = 'OT_Lin'
            else:
                exe_name = 'OT_Mac'

            db_path = self.global_settings.get_db_path()
                
            # Build paths with quotes
            program_path = f'{os.path.join(off_target_dir, exe_name)}'
            temp_path = f'{os.path.join(db_path, "temp.txt")}'
            cspr_path = f'{files["cspr_file"]}'
            db_path = f'{files["db_file"]}'
            output_path = f'{self._get_output_path(parameters)}'
            casper_info_path = f'{self.global_settings.get_casper_info_path()}'
            endo = f'{parameters["endonuclease"]}'
            
            # Build command exactly as in old version
            cmd_parts = [
                temp_path,
                endo,
                cspr_path, 
                db_path,
                output_path,
                casper_info_path,
                str(parameters['max_mismatches']),
                str(parameters['tolerance']),
                'FALSE' if parameters['average_output'] else 'TRUE',
                'TRUE' if parameters['average_output'] else 'FALSE',
                f'{self._get_hsu_value(parameters)}'
            ]
            
            return program_path, cmd_parts
            
        except Exception as e:
            self.logger.error(f"Error building command: {str(e)}")
            raise

    def _get_output_path(self, parameters):
        """Get output file path"""
        print(f"parameters.get('output_filename'): {parameters.get('output_filename')}")
        if parameters.get('output_filename'):  # If filename exists
            return os.path.join(
                self.global_settings.get_db_path(),
                parameters['output_filename']
            )
        return os.path.join(
            self.global_settings.get_off_target_dir_path(),
            'local_output.txt'
        )

    def _get_hsu_value(self, parameters):
        """Get HSU matrix name for endonuclease from CASPERinfo file"""
        try:
            casper_info_path = self.global_settings.get_casper_info_path()
            endo = parameters['endonuclease']
            
            with open(casper_info_path, 'r') as f:
                lines = f.readlines()
                
                # Find ENDONUCLEASES section
                for i, line in enumerate(lines):
                    if line.strip() == "ENDONUCLEASES":
                        # Search through endonuclease entries
                        for entry_line in lines[i+1:]:
                            if entry_line.strip() == "-----------------------------------------------------------":
                                break
                                    
                            # Parse endonuclease entry
                            parts = entry_line.strip().split(';')
                            if len(parts) >= 10:  # Make sure we have enough parts
                                endo_name = parts[1]  # The endonuclease name is the second part
                                if endo_name == endo:
                                    # Get full HSU matrix name (last part)
                                    hsu_matrix = parts[-1]
                                    # Format as "MATRIX:HSU MATRIX-{endo}-{year}"
                                    matrix_name = f"MATRIX:{hsu_matrix}"
                                    self.logger.debug(f"Found HSU matrix {matrix_name} for {endo}")
                                    return matrix_name
                                    
            # If we didn't find a match, use default
            self.logger.warning(f"Could not find HSU matrix for {endo}, using default")
            return "MATRIX:HSU MATRIX-spCas9-2013"  # Default matrix
            
        except Exception as e:
            self.logger.warning(f"Error getting HSU matrix: {str(e)}, using default")
            return "MATRIX:HSU MATRIX-spCas9-2013"  # Default matrix

    def stop_analysis(self):
        """Stop running analysis"""
        try:
            if self.process:
                self.process.kill()
                self.process = None
                
            # Clean up temp file when analysis is stopped
            self._cleanup_temp_file()
            
        except Exception as e:
            self.logger.error(f"Error stopping analysis: {str(e)}")

    def cleanup(self):
        """Clean up any temporary files"""
        try:
            # Only clean up temp file
            temp_path = os.path.join(
                self.global_settings.get_off_target_dir_path(),
                'temp.txt'
            )
            if os.path.exists(temp_path):
                os.remove(temp_path)
                self.logger.debug(f"Removed temp file: {temp_path}")
                
            # Don't remove local_output.txt here - let it be handled by the process
                
        except Exception as e:
            self.logger.error(f"Error in cleanup: {str(e)}")

    def _write_targets_to_temp(self, targets):
        """Write target sequences to temporary file for analysis"""
        try:
            db_path = self.global_settings.get_db_path()
            temp_path = os.path.join(db_path, 'temp.txt')
            
            self.logger.debug(f"Writing targets to temp file: {temp_path}")
            
            with open(temp_path, 'w') as f:
                print(targets)
                for target in targets:
                    # Format: position;sequence;pam;score;strand
                    entry = f"{target['location'].split('-')[0]};{target['sequence']};{target['pam']};{target['score']};{target['strand']}\n"
                    f.write(entry)
                        
            self.logger.debug(f"Successfully wrote {len(targets)} targets to temp file")
            
            # Store temp path for cleanup
            self._temp_path = temp_path
                        
        except Exception as e:
            self.logger.error(f"Error writing targets to temp file: {str(e)}")
            raise
