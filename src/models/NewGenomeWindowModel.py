import os
import platform
from utils.sequence_utils import SeqTranslate
from collections import OrderedDict
from PyQt6.QtCore import QObject, pyqtSignal

class NewGenomeWindowModel(QObject):

    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self.jobs = []
        self.file = ""
        self.endonucleases = self.settings.get_endonucleases()
        self.total_chrom_count = 0
        self.perc_increase = 0
        self.progress = 0
        self.total_progress = 0
        self.current_job_progress = 0
        self.total_jobs = 0
        self.completed_jobs = 0
        try:
            self.seqTrans = SeqTranslate(self.settings.get_casper_info_path())
            self.load_endonucleases_data()
        except Exception as e:
            self.logger.error(f"Error initializing NewGenomeWindowModel: {e}")

    def add_job(self, organism_name, strain, organism_code, file_path, endonuclease, multithreading_checked, generate_repeats_checked):
        # Ensure the database path exists before adding a job
        self.settings.ensure_db_path_exists()

        # This method adds a new job to the queue and checks for duplicates
        
        # Get the endonuclease data from the endonucleases dictionary
        endonuclease_data = self.endonucleases[endonuclease]
        
        # Construct the job arguments using the provided information
        args = self.create_arguments_command_for_job(organism_name, strain, organism_code, file_path, endonuclease_data, multithreading_checked, generate_repeats_checked)
        
        # Create a unique job name using the organism code and endonuclease name
        job_name = f"{organism_code}_{endonuclease_data['endonuclease_abbreviation']}"

        # Create a unique identifier for the job to prevent duplicate entries
        job_description = f"{organism_name} {strain} {endonuclease_data['endonuclease_abbreviation']} {organism_code}"

        # Combine the job name as key, and job description and arguments as values
        job_entry = {job_name: {'description': job_description, 'arguments': args}}
        
        # Check if the job already exists in the queue
        if any(job_name in job for job in self.jobs):
            return False, "Duplicate entry"
        # If the job is new, add it to the queue and the job_names
        self.jobs.append(job_entry)
        
        # Return success status and the job name
        return True, job_name
    
    def remove_job(self, job_identifier):
        for index, job in enumerate(self.jobs):
            if job_identifier in job:
                del self.jobs[index]
                return True
        return False

    def create_arguments_command_for_job(self, organism_name, strain, organism_code, file_path, endonuclease_data, multithreading_checked, generate_repeats_checked):
        db_path = self.settings.get_db_path()
        
        # Ensure db_path ends with a forward slash
        if not db_path.endswith('/'):
            db_path = f"{db_path}/"
        
        self.logger.debug(f"Using database path: {db_path}")  # Add logging
        
        print(f"The endonuclease data is {endonuclease_data}")

        arguments = [
            endonuclease_data['endonuclease_abbreviation'],
            endonuclease_data['endonuclease_pam_sequence'],
            'TRUE' if multithreading_checked else 'FALSE',
            'FALSE' if endonuclease_data['endonuclease_direction'] == '3' else 'TRUE',
            'TRUE' if generate_repeats_checked else 'FALSE',
            endonuclease_data['endonuclease_five_prime_length'], 
            endonuclease_data['endonuclease_seed_length'], 
            endonuclease_data['endonuclease_three_prime_length'], 
            organism_code,
            db_path,  # This will now always end with a forward slash
            self.settings.get_casper_info_path(),
            file_path,
            f'{organism_name} {strain}',
            'notes',
            f'DATA:{endonuclease_data["endonuclease_on_target_scoring"]}' 
        ]
        
        # Add logging of the full command
        self.logger.debug(f"Generated command arguments: {arguments}")
        
        return arguments
    
    def get_arguments_command_for_job(self, job_index):
        job_entry = self.jobs[job_index]
        job_name = next(iter(job_entry))
        return job_entry[job_name]['arguments']

    def get_job_command(self):
        if platform.system() == 'Windows':
            program = f'"{os.path.join(self.settings.get_SeqFinder_dir_path(), "Casper_Seq_Finder_Win.exe")}" '
        elif platform.system() == 'Linux':
            program = f'{os.path.join(self.settings.get_SeqFinder_dir_path(), "Casper_Seq_Finder_Lin")}'
        else:
            program = f'{os.path.join(self.settings.get_SeqFinder_dir_path(), "Casper_Seq_Finder_Mac")}'
        return program 

    def update_total_progress(self):
        if self.total_jobs > 0:
            job_weight = 100 / self.total_jobs
            self.total_progress = (self.completed_jobs * job_weight) + (self.current_job_progress * job_weight / 100)
        else:
            self.total_progress = 0
        return self.total_progress

    def reset_progress(self):
        self.total_progress = 0
        self.current_job_progress = 0
        self.completed_jobs = 0
        self.total_jobs = 0
        self.total_chrom_count = 0
        self.perc_increase = 0
        self.progress = 0

    def set_total_jobs(self, total):
        self.total_jobs = total
        self.completed_jobs = 0
        self.current_job_progress = 0
        self.total_progress = 0

    def increment_completed_jobs(self):
        self.completed_jobs += 1
        self.current_job_progress = 0
        return self.update_total_progress()

    def validate_fasta_file(self, file_path):
        return file_path.lower().endswith(('.fa', '.fna', '.fasta'))

    def get_endonuclease_info(self, endonuclease):
        print(f"Getting endonuclease info for {self.endonucleases}")
        print(f"The endonuclease is {endonuclease}")
        if not endonuclease:
            return None
        return self.endonucleases.get(endonuclease, None)

    def get_job_name(self, job_index):
        """Get the name of the job at the given index"""
        if 0 <= job_index < len(self.jobs):
            job_entry = self.jobs[job_index]
            return next(iter(job_entry))  # Returns the first (and only) key
        return None
