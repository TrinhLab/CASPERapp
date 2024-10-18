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
        self.logger = logger
        self.config_dir_path = os.path.join(self.app_dir_path, 'config')
        self.env_path = os.path.join(self.app_dir_path, '.env')
        self.config = self._initialize_config()

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
        with open(self.env_path, 'r') as f:
            lines = f.readlines()

        with open(self.env_path, 'w') as f:
            for line in lines:
                if line.startswith(f'{key}='):
                    f.write(f'{key}="{value}"\n')  # Always use double quotes
                else:
                    f.write(line)
        self.logger.info(f"Updated {key} in .env file")

    def get_env_value(self, key, default=None):
        return os.getenv(key, default)

    def set_env_value(self, key, value):
        self.write_to_env(key, value)
        os.environ[key] = value
        self.logger.info(f"Set environment variable: {key}")

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
        os.makedirs(self.config_dir, exist_ok=True)
        self.logger.info(f"Ensured config directory exists: {self.config_dir}")
