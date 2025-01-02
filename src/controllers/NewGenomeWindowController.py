from PyQt6 import QtWidgets, QtCore 
from models.NewGenomeWindowModel import NewGenomeWindowModel
from views.NewGenomeWindowView import NewGenomeWindowView
from utils.ui import show_message, show_error 
import os
from PyQt6.QtCore import pyqtSignal

class NewGenomeWindowController:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        self.directory_change_completed = pyqtSignal(str)

        try:
            self.model = NewGenomeWindowModel(self.settings)
            self.view = NewGenomeWindowView(self.settings)

            self._setup_connections()
            self._init_ui()
            self._initialize_process()
            
            self.view.endonuclease_changed.connect(self._update_endonuclease_lengths)
            
            self._load_initial_endonuclease_settings()
        except Exception as e:
            show_error(self.settings, "Error initializing NewGenomeWindowController", str(e))

        # self.model.cspr_files_created.connect(self._on_cspr_files_created)

        self.view.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding
        )

        self.settings.endonuclease_updated.connect(self._update_endonuclease_dropdown)

    def _initialize_process(self):
        self.job_process = QtCore.QProcess()
        self.job_process.setProcessChannelMode(QtCore.QProcess.ProcessChannelMode.MergedChannels)
        self.job_process.readyReadStandardOutput.connect(self._handle_process_output)
        self.job_process.finished.connect(self._handle_job_completion)

    def _handle_process_output(self):
        output = self.job_process.readAllStandardOutput().data().decode().strip()
        if output.startswith("Progress:"):
            progress = float(output.split(":")[1].strip())
            self.model.update_current_job_progress(progress)
            total_progress = self.model.update_total_progress()
            self.view.set_progress_bar_jobs(total_progress)

    def _setup_connections(self):
        # grpStep2
        self.view.push_button_ncbi_search.clicked.connect(self._open_ncbi_module)
        self.view.push_button_browse_file.clicked.connect(self._browse_fasta_file)
        self.view.push_button_add_job.clicked.connect(self._add_job_to_table)

        # grpStep3
        self.view.push_button_remove_all_jobs.clicked.connect(self._reset_table_widget_jobs)
        self.view.push_button_remove_job.clicked.connect(self._remove_selected_job)
        self.view.push_button_run_all_jobs.clicked.connect(self._run_all_jobs)

        # bottom of the window
        self.view.push_button_reset_form.clicked.connect(self._handle_reset)

    def _init_ui(self):
        self.view.update_endonuclease_dropdown(self.model.endonucleases.keys())
        # Remove this line as we'll call it in _update_endonuclease_lengths
        # self._load_endonuclease_settings()

    def _handle_reset(self):
        try:
            self.view.line_edit_organism_name.clear()
            self.view.line_edit_strain.clear()
            self.view.line_edit_organism_code.clear()

            self.model.file = ""
            self.view.line_edit_selected_file.clear()
            self.view.line_edit_selected_file.setPlaceholderText("Selected FASTA/FNA File")

            self.view.reset_table_widget_jobs()
            self.view.reset_progress_bar_jobs()

            # Reinitialize the process
            if self.job_process.state() != QtCore.QProcess.ProcessState.NotRunning:
                self.job_process.kill()
            self._initialize_process()

            # Reset the model
            self.model.reset_progress()
            self.model.jobs.clear()

            # Cancel any pending database path change
            self.settings.db_manager.pending_db_path = None
            self.logger.debug("Cancelled pending database path change")
            
            # Show confirmation to user
            show_message("Reset Complete", 
                        "Form has been reset and any pending database changes have been cancelled.")
        except Exception as e:
            self.logger.error(f"Error in handle reset: {str(e)}")
            show_error(self.settings, "Error", str(e))

    def _add_job_to_table(self):
        organism_name = self.view.get_organism_name()
        strain_name = self.view.get_strain()
        organism_code = self.view.get_organism_code()
        file_path = self.view.get_selected_file()

        warnings = self._validate_job_input(organism_name, strain_name, organism_code, file_path)
        if warnings:
            show_message(fontSize=12, icon=QtWidgets.QMessageBox.Icon.Critical,
                         title="Required Information", message="\n".join(warnings))
            return

        endonuclease = self.view.get_selected_endonuclease()
        multithreading_checked = self.view.is_multithreading_checked()
        generate_repeats_checked = self.view.is_generate_repeats_checked()

        success, job_name = self.model.add_job(organism_name, strain_name, organism_code, file_path, endonuclease, multithreading_checked, generate_repeats_checked)
        if success:
            self.view.add_job_to_table(job_name)
        else:
            show_message(fontSize=12, icon=QtWidgets.QMessageBox.Icon.Critical,
                         title="Duplicate Entry", message=job_name)

    def _validate_job_input(self, org_name, strain_name, org_code, file_path):
        warnings = []
        if not org_name:
            warnings.append("You need to include the organism's name.\n")
        if not file_path:
            warnings.append("You need to select a file.\n")
        if not strain_name:
            warnings.append("It is recommended to include the organism's subspecies/strain.\n")
        if not org_code:
            warnings.append("You must include an organism code.\n")
        return warnings

    def _browse_fasta_file(self):
        file_dialog = QtWidgets.QFileDialog()
        database_dir = self.settings.db_manager.get_active_db_path()
        self.logger.debug(f"Opening file dialog with directory: {database_dir}")
        
        file_path, _ = file_dialog.getOpenFileName(
            self.view, 
            "Choose a File", 
            database_dir, 
            "FASTA Files (*.fa *.fna *.fasta)"
        )
        
        if file_path:
            if self.model.validate_fasta_file(file_path):
                self.model.file = file_path
                self.view.set_selected_file(file_path)
                self.logger.debug(f"Selected valid FASTA file: {file_path}")
            else:
                show_message(
                    fontSize=12, 
                    icon=QtWidgets.QMessageBox.Icon.Critical,
                    title="File Selection Error",
                    message="You have selected an incorrect type of file. Please choose a FASTA/FNA file."
                )

    def _remove_selected_job(self):
        job_identifier = self.view.get_selected_job_identifier()
        if job_identifier:
            success = self.model.remove_job(job_identifier)
            if success:
                self.view.remove_selected_job_from_table()
                show_message(fontSize=12, icon=QtWidgets.QMessageBox.Icon.Information,
                             title="Job Removed", message="The selected job has been removed.")
            else:
                show_message(fontSize=12, icon=QtWidgets.QMessageBox.Icon.Critical,
                             title="Error", message="Failed to remove the selected job.")
        else:
            show_message(fontSize=12, icon=QtWidgets.QMessageBox.Icon.Warning,
                         title="No Selection", message="Please select a job to remove.")

    def _update_endonuclease_lengths(self, endonuclease):
        print(f"Updating endonuclease lengths for {endonuclease}")
        if not endonuclease:  # If endonuclease is empty, don't proceed
            return
        
        endonuclease_info = self.model.get_endonuclease_info(endonuclease)
        self.logger.debug(f"Endonuclease info: {endonuclease_info}")
        
        if endonuclease_info is not None:
            self.view.set_endonuclease_lengths(
                endonuclease_info['endonuclease_seed_length'],
                endonuclease_info['endonuclease_five_prime_length'],
                endonuclease_info['endonuclease_three_prime_length']
            )
        else:
            self.logger.warning(f"No information found for endonuclease: {endonuclease}")
            # Optionally, clear the length fields here as well
            self.view.set_endonuclease_lengths('', '', '')

    def _run_all_jobs(self):
        self.logger.debug("Starting to run all queued jobs")
        queued_indexes = self.view.get_queued_job_indexes()
        print(f"Queued job indexes: {queued_indexes}")
        
        if not queued_indexes:
            show_message(fontSize=12, icon=QtWidgets.QMessageBox.Icon.Information,
                         title="No Jobs To Run",
                         message="No queued jobs are available to run.")
            return

        self.model.reset_progress()
        self.model.set_total_jobs(len(queued_indexes))
        self.view.reset_progress_bar_jobs()
        
        self.job_indexes = queued_indexes
        if self.job_indexes:
            first_row_index = self.job_indexes[0]
            self._run_job(first_row_index)

        self.logger.debug(f"Queued {len(self.job_indexes)} jobs for running")

    def _run_job(self, row_index):
        self.logger.debug(f"Running job at row {row_index}")
        
        # Start the spinner for this job
        self.view.start_spinner(row_index)

        program = self.model.get_job_command()
        command_args = self.model.get_arguments_command_for_job(row_index)
        self.logger.debug(f"Executing command: {program} {' '.join(command_args)}")
        self.logger.debug(f"Working directory: {os.getcwd()}")

        if self.job_process.state() == QtCore.QProcess.ProcessState.NotRunning:
            # Set up process output handling
            def handle_stdout():
                output = self.job_process.readAllStandardOutput().data().decode()
                self.logger.debug(f"Process stdout: {output}")

            def handle_stderr():
                error = self.job_process.readAllStandardError().data().decode()
                self.logger.error(f"Process stderr: {error}")

            self.job_process.readyReadStandardOutput.connect(handle_stdout)
            self.job_process.readyReadStandardError.connect(handle_stderr)

            # Start the process
            self.job_process.start(program, command_args)
            
            # Check if process started successfully
            if not self.job_process.waitForStarted(3000):  # 3 second timeout
                self.logger.error("Process failed to start")
                self.logger.error(f"Process error: {self.job_process.errorString()}")
                return
            
            self.logger.debug("Job started")
        else:
            self.logger.warning("Process is still running, cannot start a new job.")
            # Kill the current process and start the new one:
            self.job_process.kill()
            self._initialize_process()
            self.job_process.start(program, command_args)

        # Ensure the table widget updates its display
        self.view.table_widget_jobs.viewport().update()

    def _handle_job_completion(self, exit_code=None, exit_status=None):
        """Handle process completion and job status updates"""
        try:
            self.logger.debug(f"Process finished with exit code: {exit_code}")
            
            # Log any remaining output
            remaining_output = self.job_process.readAllStandardOutput().data().decode()
            if remaining_output:
                self.logger.debug(f"Final process output: {remaining_output}")
            
            remaining_error = self.job_process.readAllStandardError().data().decode()
            if remaining_error:
                self.logger.error(f"Final process error output: {remaining_error}")
            
            if hasattr(self, 'job_indexes') and self.job_indexes:
                completed_row_index = self.job_indexes.pop(0)
                
                # Get the job name for the completed job
                job_name = self.model.get_job_name(completed_row_index)
                expected_cspr_file = os.path.join(self.settings.get_db_path(), f"{job_name}.cspr")
                
                if exit_code == 0 and os.path.exists(expected_cspr_file):
                    self.logger.debug(f"CSPR file created successfully: {expected_cspr_file}")
                    self.view.set_job_completed(completed_row_index)
                    
                    # Update model's completed jobs count and progress bar
                    total_progress = self.model.increment_completed_jobs()
                    if total_progress is not None:
                        self.view.set_progress_bar_jobs(total_progress)
                    else:
                        self.logger.warning("Received None for total_progress")
                    
                    # If we're in directory change mode, trigger a database state update
                    if self.settings.db_manager.is_changing_directory:
                        self.logger.debug("Directory change in progress - triggering state update")
                        self.settings.db_manager.update_db_state()
                else:
                    error_msg = f"Job failed: CSPR file not found or process error (exit code: {exit_code})"
                    self.logger.error(error_msg)
                    show_error(self.settings, f"Job Failed: {job_name}", error_msg)
                    
                # Process next job if any
                if self.job_indexes:
                    next_row_index = self.job_indexes[0]
                    self._run_job(next_row_index)
                else:
                    self.logger.info("All queued jobs completed")
                    self.view.set_progress_bar_jobs(100)
                    
                    # If we were changing directory, finalize the change
                    if self.settings.db_manager.is_changing_directory:
                        self.logger.debug("Finalizing directory change after successful genome analysis")
                        success, message = self.settings.db_manager.finalize_directory_change()
                        if success:
                            self.directory_change_completed.emit(message)
                        else:
                            self.logger.error(f"Failed to finalize directory change: {message}")
                    else:
                        # Just update the state if not changing directory
                        self.settings.update_db_state()
                    
            else:
                self.logger.warning("No job indexes found or all jobs completed")

            self.view.table_widget_jobs.viewport().update()
            
        except Exception as e:
            self.logger.error(f"Error in _handle_job_completion: {str(e)}")

    def _reset_table_widget_jobs(self):
        self.view.reset_table_widget_jobs()
        self.model.jobs.clear()  # Clear the jobs in the model as well

    def _open_ncbi_module(self):
        try:
            # Get organism and strain values from the view
            organism_name = self.view.get_organism_name()
            strain_name = self.view.get_strain()
            
            # Get NCBI controller
            ncbi_controller = self.settings.get_ncbi_window()
            
            # Connect to the initialization complete signal
            def on_init_complete():
                self.logger.debug(f"NCBI window initialized, setting organism: {organism_name}, strain: {strain_name}")
                if organism_name:
                    ncbi_controller.view.line_edit_organism.setText(organism_name)
                if strain_name:
                    ncbi_controller.view.line_edit_strain.setText(strain_name)
                # Disconnect after use to prevent multiple connections
                ncbi_controller.view.initialization_complete.disconnect(on_init_complete)
            
            # Connect the signal
            ncbi_controller.view.initialization_complete.connect(on_init_complete)
            
            # Open the NCBI tab
            self.settings.main_window.open_new_tab("NCBI Download Tool", ncbi_controller)
            
        except Exception as e:
            show_error(self.settings, "Error opening NCBI module", str(e))
            self.logger.error(f"Failed to open NCBI module: {str(e)}")
        
    def _load_initial_endonuclease_settings(self):
        initial_endonuclease = self.view.get_selected_endonuclease()
        print(f"Initial endonuclease: {initial_endonuclease}")
        self._update_endonuclease_lengths(initial_endonuclease)

    # def _on_cspr_files_created(self):
        # self.settings.update_db_state()

    def _update_endonuclease_dropdown(self):
        print(f"Updating endonuclease dropdown, {self.model.endonucleases.keys()}")
        self.view.update_endonuclease_dropdown(self.model.endonucleases.keys())
        # If there's a currently selected endonuclease, update its info
        # current_endo = self.view.get_selected_endonuclease()
        # if current_endo:
            # self._update_endonuclease_lengths(current_endo)
        pass
