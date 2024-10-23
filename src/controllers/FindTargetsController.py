from models.FindTargetsModel import FindTargetsModel
from views.FindTargetsView import FindTargetsView
from PyQt6.QtWidgets import QMessageBox

class FindTargetsController:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.model = FindTargetsModel(global_settings)
        self.view = FindTargetsView(global_settings)
        self._connect_signals()
        self.organism = None
        self.endonuclease = None

    def _connect_signals(self):
        self.view.push_button_view_targets.clicked.connect(self.view_targets)

    def find_targets(self, input_data):
        try:
            self.global_settings.logger.debug(f"FindTargetsController received input data: {input_data}")
            self.organism = input_data['organism']
            self.endonuclease = input_data['endonuclease']
            results = self.model.find_targets(input_data)
            self.global_settings.logger.debug(f"Find targets results: {results}")
            self.view.display_results(results)
        except ValueError as e:
            error_message = str(e)
            QMessageBox.critical(self.view, "Error", error_message)
            self.global_settings.logger.error(f"ValueError in find_targets: {error_message}")
        except Exception as e:
            error_message = f"An unexpected error occurred while finding targets: {str(e)}"
            QMessageBox.critical(self.view, "Error", error_message)
            self.global_settings.logger.error(f"Unexpected error in find_targets: {str(e)}", exc_info=True)

    def view_targets(self):
        selected_targets = self.view.get_selected_targets()
        if not selected_targets:
            QMessageBox.warning(self.view, "No Selection", "Please select targets to view.")
            return
        
        view_targets_controller = self.global_settings.get_view_targets_window()
        view_targets_controller.load_targets(selected_targets, self.organism, self.endonuclease)
        self.global_settings.main_window.open_new_tab("View Targets", view_targets_controller)

    def show(self):
        self.view.show()
