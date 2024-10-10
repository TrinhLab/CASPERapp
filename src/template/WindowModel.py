import os
import glob
from typing import List, Dict, Any
from utils.ui import show_error

class WindowModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self.data = {}
        self.settings = {}

    def load_data(self) -> None:
        """Load initial data required for the window."""
        try:
            # Example: Load data from a file or database
            self._load_settings()
            self._load_items()
        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")
            show_error(self.global_settings, "Error Loading Data", str(e))

    def _load_settings(self) -> None:
        """Load settings for the window."""
        try:
            # Example: Load settings from a configuration file
            settings_file = os.path.join(self.global_settings.get_config_dir(), "window_settings.json")
            # Implement actual loading logic here
            self.settings = {"setting1": "value1", "setting2": "value2"}
        except Exception as e:
            self.logger.error(f"Error loading settings: {str(e)}")

    def _load_items(self) -> None:
        """Load items for combo boxes or lists."""
        try:
            # Example: Load items from a directory
            items_dir = os.path.join(self.global_settings.get_data_dir(), "items")
            self.data["items"] = [f for f in os.listdir(items_dir) if f.endswith(".item")]
        except Exception as e:
            self.logger.error(f"Error loading items: {str(e)}")

    def get_items(self) -> List[str]:
        """Get the list of items."""
        return self.data.get("items", [])

    def process_input(self, input_text: str) -> Dict[str, Any]:
        """Process input from the view."""
        try:
            # Example: Process the input text
            words = input_text.split()
            word_count = len(words)
            char_count = len(input_text)
            return {"word_count": word_count, "char_count": char_count}
        except Exception as e:
            self.logger.error(f"Error processing input: {str(e)}")
            return {}

    def save_data(self, data: Dict[str, Any]) -> bool:
        """Save data to a file or database."""
        try:
            # Example: Save data to a file
            output_file = os.path.join(self.global_settings.get_data_dir(), "output.json")
            # Implement actual saving logic here
            self.logger.info(f"Data saved to {output_file}")
            return True
        except Exception as e:
            self.logger.error(f"Error saving data: {str(e)}")
            return False

    def update_setting(self, key: str, value: Any) -> None:
        """Update a setting."""
        self.settings[key] = value
        self._save_settings()

    def _save_settings(self) -> None:
        """Save current settings."""
        try:
            # Example: Save settings to a configuration file
            settings_file = os.path.join(self.global_settings.get_config_dir(), "window_settings.json")
            # Implement actual saving logic here
            self.logger.info("Settings saved successfully")
        except Exception as e:
            self.logger.error(f"Error saving settings: {str(e)}")

    def validate_data(self, data: Dict[str, Any]) -> bool:
        """Validate data before processing or saving."""
        # Implement validation logic
        return True

    # Add any additional methods for data processing, validation, or business logic
