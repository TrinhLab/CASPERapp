import os
import yaml
from dotenv import load_dotenv
from collections import OrderedDict
from PyQt6.QtCore import QObject, pyqtSignal

def represent_ordereddict(dumper, data):
    return dumper.represent_mapping('tag:yaml.org,2002:map', data.items())

yaml.add_representer(OrderedDict, represent_ordereddict)

class NewlineDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()

class ConfigManager(QObject):
    config_updated = pyqtSignal()
    env_file_created = pyqtSignal()

    def __init__(self, app_dir_path, logger):
        super().__init__()
        self.app_dir_path = app_dir_path
        self.endonucleases = OrderedDict()
        self.logger = logger
        self.config_dir_path = os.path.join(self.app_dir_path, 'config')
        self.casper_info_path = os.path.join(self.config_dir_path, 'CASPERinfo')
        self.env_path = os.path.join(self.app_dir_path, '.env')
        self.config = self._initialize_config()
        self.load_endonucleases_data()

    def _initialize_config(self):
        config_path = os.path.join(self.config_dir_path, 'config.yml')
        if not os.path.exists(config_path):
            self.logger.info(f"Config file not found. Creating default config at {config_path}")
            self._create_default_config(config_path)
        return self._load_config(config_path)

    def _create_default_config(self, config_path):
        config = OrderedDict([
            ('paths', OrderedDict([
                ('ui', 'ui'),
                ('controllers', 'controllers'),
                ('assets', 'assets'),
                ('casper_info', 'CASPERinfo')
            ])),
            ('algorithms', ['Azimuth 2.0'])
        ])
        
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w') as f:
            yaml.dump(config, f, Dumper=NewlineDumper, default_flow_style=False)
        self.logger.info(f"Default config created at {config_path}")

    def _load_config(self, config_path):
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            self.logger.error(f"Error loading config from {config_path}: {e}")
            return None

    def save_config(self):
        config_path = os.path.join(self.config_dir_path, 'config.yml')
        config = OrderedDict([
            ('paths', OrderedDict([
                ('ui', self.config['paths']['ui']),
                ('controllers', self.config['paths']['controllers']),
                ('assets', self.config['paths']['assets']),
                ('casper_info', self.config['paths']['casper_info'])
            ])),
            ('algorithms', self.config.get('algorithms', ["Azimuth 2.0"]))
        ])
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f, Dumper=NewlineDumper, default_flow_style=False)
        
        self.logger.info(f"Config saved to {config_path}")
        self.config_updated.emit()

    def load_env(self):
        if not os.path.exists(self.env_path):
            self._create_default_env()
        load_dotenv(self.env_path)

    def _create_default_env(self):
        with open(self.env_path, 'w') as f:
            f.write(f'APP_DIR="{self.app_dir_path}"\n')
            f.write('CSPR_DB=""\n')
            f.write('LOG_LEVEL=DEBUG\n')
            f.write('FIRST_TIME_START=TRUE\n')  # Add this line
        self.logger.info(f"Default .env file created at {self.env_path}")
        self.env_file_created.emit()

    def write_to_env(self, key, value):
        """Write or update a key-value pair in the .env file"""
        try:
            # Read all lines
            if os.path.exists(self.env_path):
                with open(self.env_path, 'r') as f:
                    lines = f.readlines()
            else:
                lines = []

            # Find if key exists
            key_exists = False
            new_lines = []
            for line in lines:
                if line.strip().startswith(f'{key}='):
                    new_lines.append(f'{key}="{value}"\n')  # Always use double quotes
                    key_exists = True
                else:
                    new_lines.append(line)

            # If key doesn't exist, append it
            if not key_exists:
                new_lines.append(f'{key}="{value}"\n')

            # Write back to file
            with open(self.env_path, 'w') as f:
                f.writelines(new_lines)

            # Update environment variable in memory
            os.environ[key] = value
            
            self.logger.info(f"Updated {key}={value} in .env file")
            
            # Verify the write
            with open(self.env_path, 'r') as f:
                content = f.read()
                if f'{key}="{value}"' not in content:
                    self.logger.error(f"Failed to verify {key}={value} in .env file")
                    raise Exception("Failed to verify environment variable update")

        except Exception as e:
            self.logger.error(f"Error writing to .env file: {str(e)}")
            raise

    def get_env_value(self, key, default=None):
        return os.getenv(key, default)

    def set_env_value(self, key, value):
        """Set and save an environment variable"""
        try:
            # Write to .env file first
            self.write_to_env(key, value)
            
            # Update runtime environment
            os.environ[key] = value
            
            # Verify the update
            if os.getenv(key) != value:
                self.logger.error(f"Failed to verify environment variable update for {key}")
                raise Exception(f"Environment variable verification failed for {key}")
                
            self.logger.info(f"Successfully set environment variable: {key}={value}")
            
        except Exception as e:
            self.logger.error(f"Error setting environment variable {key}: {str(e)}")
            raise

    def get_config_value(self, key, default=None):
        keys = key.split('.')
        value = self.config
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        return value

    def set_config_value(self, key, value):
        keys = key.split('.')
        config = self.config
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        config[keys[-1]] = value
        self.save_config()
        self.logger.info(f"Set config value: {key}")

    def ensure_config_dir_exists(self):
        os.makedirs(self.config_dir_path, exist_ok=True)
        self.logger.info(f"Ensured config directory exists: {self.config_dir_path}")

    def load_endonucleases_data(self):
        self.endonucleases.clear()  # Clear existing data before reloading
        try:
            with open(self.casper_info_path, 'r') as casper_info_file:
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
        print(fields)
        if len(fields) == 10:
            endonuclease_abbreviation = fields[0]
            pam_sequence = fields[1]
            five_prime_length, seed_length, three_prime_length, direction = fields[2:6]
            endonuclease_name = fields[6]
            endonuclease_CRISPR_type = fields[7]
            endonuclease_on_target_matrix = fields[8]
            endonuclease_off_target_matrix = fields[9]

            endonuclease_key = f"{endonuclease_abbreviation} - PAM: {pam_sequence}"
            endonuclease_value = {
                'endonuclease_abbreviation': endonuclease_abbreviation,
                'endonuclease_pam_sequence': pam_sequence,
                'endonuclease_five_prime_length': five_prime_length,
                'endonuclease_seed_length': seed_length,
                'endonuclease_three_prime_length': three_prime_length,
                'endonuclease_direction': direction,
                'endonuclease_organism': endonuclease_name,
                'endonuclease_CRISPR_type': endonuclease_CRISPR_type,
                'endonuclease_on_target_scoring': endonuclease_on_target_matrix,
                'endonuclease_off_target_scoring': endonuclease_off_target_matrix
            }
            self.endonucleases[endonuclease_key] = endonuclease_value

    def _get_primary_pam(self, pam_field):
        return pam_field.split(',')[0] if ',' in pam_field else pam_field
    
    def get_endonucleases(self):
        return self.endonucleases
