from models.AnnotationParser import AnnotationParser
import os

class BaseModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        
        # Initialize the annotation parser
        self._initialize_annotation_parser()
        
        # Connect to annotation file changes
        self.global_settings.annotation_file_changed.connect(self._on_annotation_file_changed)

    def _initialize_annotation_parser(self) -> None:
        """Initialize or get existing annotation parser from global settings"""
        try:
            # Try to get existing parser from global settings
            if not hasattr(self.global_settings, 'annotation_parser'):
                # Create new parser if none exists
                annotation_file = self.global_settings.get_current_annotation_file()
                annotation_path = os.path.join(
                    self.global_settings.get_db_path(),
                    'GBFF',
                    annotation_file
                )
                self.global_settings.annotation_parser = AnnotationParser(self.global_settings)
                self.global_settings.annotation_parser.set_annotation_file(annotation_path)
                self.logger.debug("Created new annotation parser")
            
            self.annotation_parser = self.global_settings.annotation_parser
            self.annotation_path = self.annotation_parser.annotation_file_name
            
        except Exception as e:
            self.logger.error(f"Error initializing annotation parser: {str(e)}")
            raise

    def _on_annotation_file_changed(self, new_annotation_file):
        """Handle annotation file changes"""
        try:
            self.logger.debug(f"Clearing caches for new annotation file: {new_annotation_file}")
            self._clear_caches()
            self._initialize_annotation_parser()
        except Exception as e:
            self.logger.error(f"Error handling annotation file change: {str(e)}")

    def _clear_caches(self):
        """Clear model-specific caches. Override in subclasses."""
        pass

    def cleanup(self):
        """Cleanup resources. Override in subclasses if needed."""
        try:
            self.global_settings.annotation_file_changed.disconnect(self._on_annotation_file_changed)
            self.logger.debug("Disconnected from annotation file changes")
            self._clear_caches()
        except Exception as e:
            self.logger.error(f"Error in cleanup: {str(e)}") 