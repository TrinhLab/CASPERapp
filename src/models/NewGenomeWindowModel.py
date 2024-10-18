import os
import platform
from utils.sequence_utils import SeqTranslate
from collections import OrderedDict
from PyQt6.QtCore import QObject, pyqtSignal

class NewGenomeWindowModel(QObject):

    def __init__(self, global_settings):
        super().__init__()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        self.jobs = []
        self.file = ""
        self.endonucleases = OrderedDict()
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

    def load_endonucleases_data(self):
        try:
            with open(self.settings.get_casper_info_path(), 'r') as casper_info_file:
                for line in casper_info_file:
                    if line.startswith('ENDONUCLEASES'):
                        self._parse_endonuclease_data(casper_info_file)
                        break
        except Exception as e:
            self.logger.error(f"Error loading endonucleases: {e}")

    def _parse_endonuclease_data(self, file):
        for line in file:
            if line.startswith('-'):
                break
            self._process_endonuclease_line(line)

    def _process_endonuclease_line(self, line):
        fields = line.strip().split(';')
        if len(fields) == 10:
            endonuclease_name = fields[0]
            pam_sequence = self._get_primary_pam(fields[1])
            five_prime_length, seed_length, three_prime_length, direction = fields[2:6]
            on_target_data = fields[8]

            endonuclease_key = f"{endonuclease_name} - PAM: {pam_sequence}"
            endonuclease_value = {
                'name': endonuclease_name,
                'pam': pam_sequence,
                'five_prime_length': five_prime_length,
                'seed_length': seed_length,
                'three_prime_length': three_prime_length,
                'direction': direction,
                'on_target_data': on_target_data
            }
            self.endonucleases[endonuclease_key] = endonuclease_value

    def _get_primary_pam(self, pam_field):
        return pam_field.split(',')[0] if ',' in pam_field else pam_field

    def add_job(self, organism_name, strain, organism_code, file_path, endonuclease, multithreading_checked, generate_repeats_checked):
        # Ensure the database path exists before adding a job
        self.settings.ensure_db_path_exists()

        # This method adds a new job to the queue and checks for duplicates
        
        # Get the endonuclease data from the endonucleases dictionary
        endonuclease_data = self.endonucleases[endonuclease]
        
        # Construct the job arguments using the provided information
        args = self.create_arguments_command_for_job(organism_name, strain, organism_code, file_path, endonuclease_data, multithreading_checked, generate_repeats_checked)
        
        # Create a unique job name using the organism code and endonuclease name
        job_name = f"{organism_code}_{endonuclease_data['name']}"

        # Create a unique identifier for the job to prevent duplicate entries
        job_description = f"{organism_name} {strain} {endonuclease_data['name']} {organism_code}"

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
        for index, job in enumerate(self.job_names):
            if job_identifier in job:
                del self.jobs[index]
                self.job_names.pop(index)
                return True
        return False

    def create_arguments_command_for_job(self, organism_name, strain, organism_code, file_path, endonuclease_data, multithreading_checked, generate_repeats_checked):
        db_path = self.settings.get_db_path()
        
        # Preserve trailing slash if present
        if db_path.endswith(os.path.sep):
            db_path = db_path.rstrip(os.path.sep) + os.path.sep
        # Add trailing slash for Darwin (macOS) machines
        elif platform.system() == 'Darwin':
            db_path = db_path.rstrip('/') + '/'

        arguments = [
            endonuclease_data['name'],  # endo
            endonuclease_data['pam'],  # pam
            'TRUE' if multithreading_checked else 'FALSE',
            'FALSE' if endonuclease_data['direction'] == '3' else 'TRUE',
            'TRUE' if generate_repeats_checked else 'FALSE',
            endonuclease_data['five_prime_length'],  # five length
            endonuclease_data['seed_length'],  # seed length
            endonuclease_data['three_prime_length'],  # three length
            organism_code,
            f'{db_path}',
            f'{self.settings.get_casper_info_path()}',
            f'{file_path}',
            f'{organism_name} {strain}',
            'notes',
            f'"DATA:{endonuclease_data["on_target_data"]}"'  # on target matrix
        ]
        return arguments
    
    def get_arguments_command_for_job(self, job_index):
        job_entry = self.jobs[job_index]
        job_name = next(iter(job_entry))
        return job_entry[job_name]['arguments']

    def get_job_command(self):
        if platform.system() == 'Windows':
            program = f'"{os.path.join(self.settings.get_SeqFinder_dir_path(), "Casper_Seq_Finder_Win.exe")}" '
        elif platform.system() == 'Linux':
            program = f'"{os.path.join(self.settings.get_SeqFinder_dir_path(), "Casper_Seq_Finder_Lin")}" '
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
        return self.endonucleases.get(endonuclease, (None, None, None, None, None))
