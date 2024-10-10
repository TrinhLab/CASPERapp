from PyQt6 import QtWidgets
from views.NCBIRenameWindowUI import NCBIRenameWindowUI
from utils.ui import show_message
import os
import platform
from models.NCBIRenameWindowModel import NCBIRenameWindowModel

class NCBIRenameWindowController:
    def __init__(self, settings, files, parent_controller):
        self.settings = settings
        self.parent_controller = parent_controller
        self.model = NCBIRenameWindowModel(settings)
        self.model.set_files(files)
        self.view = NCBIRenameWindowUI(settings, files)
        self.logger = settings.get_logger()
        
        self.setup_connections()

    def setup_connections(self):
        self.view.submit_button.clicked.connect(self.submit_rename)
        self.view.go_back.clicked.connect(self.go_back)

    def show(self):
        self.view.show()

    def go_back(self):
        self.view.close()

    def submit_rename(self):
        try:
            new_names = self.view.get_new_names()
            results = self.model.process_rename_batch(new_names)
            
            success = all(result[0] for result in results)
            if success:
                self.view.clear_table()
                self.view.close()
                if hasattr(self.parent_controller, 'on_rename_complete'):
                    self.parent_controller.on_rename_complete()
                else:
                    self.logger.warning("Parent controller does not have on_rename_complete method")
                    # Fallback behavior
                    if hasattr(self.parent_controller, 'update_main_window'):
                        self.parent_controller.update_main_window()
                    if hasattr(self.parent_controller, 'show_navigation_prompt'):
                        self.parent_controller.show_navigation_prompt()
            else:
                error_messages = [f"{old_name} -> {new_name}: {message}" 
                                  for success, old_name, new_name, message in results if not success]
                show_message(
                    fontSize=12,
                    icon=QtWidgets.QMessageBox.Icon.Warning,
                    title="Renaming Errors",
                    message="The following errors occurred during renaming:\n\n" + "\n".join(error_messages)
                )
        except Exception as e:
            self.logger.error(f"Error in submit_rename: {str(e)}", exc_info=True)
            show_message(
                fontSize=12,
                icon=QtWidgets.QMessageBox.Icon.Critical,
                title="Renaming Error",
                message=f"An error occurred while renaming files: {str(e)}"
            )

    def exec(self):
        return self.view.exec()
