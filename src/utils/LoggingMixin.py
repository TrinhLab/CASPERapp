import logging
from typing import Optional

class LoggingMixin:
    """Mixin class to provide logging functionality to classes"""
    
    def __init__(self):
        self._logger: Optional[logging.Logger] = None
        self._init_logger()

    def _init_logger(self) -> None:
        """Initialize logger for the class"""
        self._logger = logging.getLogger(self.__class__.__name__)

    @property
    def logger(self) -> logging.Logger:
        """Get the logger instance"""
        if self._logger is None:
            self._init_logger()
        return self._logger

    def log_method_call(self, method_name: str, *args, **kwargs) -> None:
        """Log method calls with their arguments"""
        self.logger.debug(f"Calling {method_name} with args: {args}, kwargs: {kwargs}")

    def log_error(self, method_name: str, error: Exception) -> None:
        """Log errors with method context"""
        self.logger.error(f"Error in {method_name}: {str(error)}", exc_info=True)

    def log_info(self, message: str) -> None:
        """Log info messages"""
        self.logger.info(message)

    def log_debug(self, message: str) -> None:
        """Log debug messages"""
        self.logger.debug(message)

    def log_warning(self, message: str) -> None:
        """Log warning messages"""
        self.logger.warning(message) 