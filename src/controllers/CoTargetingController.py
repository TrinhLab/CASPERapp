from models.CoTargetingModel import CoTargetingModel
from views.CoTargetingView import CoTargetingView

class CoTargetingController:
    def __init__(self, global_settings, view_targets_controller=None):
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        self.model = CoTargetingModel(global_settings)
        self.view = CoTargetingView(global_settings)
        self.view_targets_controller = view_targets_controller
        
        # Connect signals
        self.view.push_button_cancel.clicked.connect(self.cancel)
        self.view.push_button_submit.clicked.connect(self.submit)

    def show(self):
        """Show the co-targeting window"""
        self.view.show()
        self.view.activateWindow()

    def launch(self, endo_choices, org_name):
        """Launch co-targeting analysis"""
        try:
            self.view.line_edit_organism.setText(org_name)
            self.view.populate_table(endo_choices)
            self.show()
            self.view.activateWindow()
        except Exception as e:
            self.logger.error(f"Error launching co-targeting: {str(e)}")
            self.view.show_error("Launch Error", str(e))

    def submit(self):
        try:
            selected_endos = self.view.get_selected_endonucleases()
            
            if len(selected_endos) <= 1:
                self.view.show_error(
                    "Nothing Selected",
                    "No endonucleases selected. Please select at least 2 endonucleases"
                )
                return
                
            # Validate compatibility
            if not self.model.validate_endonucleases(selected_endos):
                self.view.show_error(
                    "Invalid Endonucleases",
                    "The selected endonucleases are not compatible."
                )
                return
                
            # Update view targets controller with selected endonucleases
            if self.view_targets_controller:
                self.view_targets_controller.handle_cotargeting_result(selected_endos)
            else:
                self.logger.error("No view_targets_controller available")
                self.view.show_error("Error", "Could not update targets view")
                return
                
            self.cancel()
            
        except Exception as e:
            self.logger.error(f"Error in submit: {str(e)}")
            self.view.show_error("Submit Error", str(e))

    def cancel(self):
        """Handle cancel button click"""
        try:
            self.view.clear()
            self.view.hide()
        except Exception as e:
            self.logger.error(f"Error in cancel: {str(e)}")
