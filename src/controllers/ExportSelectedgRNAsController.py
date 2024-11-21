import os
import platform
from PyQt6.QtWidgets import QFileDialog
from utils.ui import show_message, show_error
from models.ExportSelectedgRNAsModel import ExportSelectedgRNAsModel
from views.ExportSelectedgRNAsView import ExportSelectedgRNAsView

class ExportSelectedgRNAsController:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        try:
            self.view = ExportSelectedgRNAsView(self.settings)
            self.model = ExportSelectedgRNAsModel(self.settings)
            self._setup_connections()
        except Exception as e:
            show_error(self.settings, "Error initializing ExportSelectedgRNAsController", str(e))

    def _setup_connections(self) -> None:
        try:
            self.view.push_button_export.clicked.connect(self._handle_export)
            self.view.push_button_cancel.clicked.connect(self._handle_cancel)
            self.view.push_button_browse.clicked.connect(self._handle_browse)
        except Exception as e:
            self.logger.error(f"Error setting up connections: {str(e)}")
            raise

    def show_dialog(self, selected_items: list, window_type: str) -> None:
        try:
            # Reset the form fields
            self.view.line_edit_file_name.clear()
            self.view.line_edit_leading_sequence.clear()
            self.view.line_edit_trailing_sequence.clear()
            
            # Set default path
            default_path = self.settings.get_db_path()
            self.view.set_file_path(default_path)
            
            # Set new data
            self.model.set_export_data(selected_items, window_type)
            
            # Show and bring to front
            self.view.show()
            self.view.raise_()
            self.view.activateWindow()
            
        except Exception as e:
            show_error(self.settings, "Error showing export dialog", str(e))

    def _handle_browse(self) -> None:
        try:
            directory = QFileDialog.getExistingDirectory(
                self.view,
                "Select Export Directory",
                self.settings.CSPR_DB,
                QFileDialog.Option.ShowDirsOnly
            )
            
            if directory:
                directory += "\\" if platform.system() == "Windows" else "/"
                self.view.set_file_path(directory)
        except Exception as e:
            show_error(self.settings, "Error browsing for directory", str(e))

    def _handle_export(self) -> None:
        try:
            settings = self.view.get_export_settings()
            
            if not settings['file_name']:
                settings['file_name'] = "exported_gRNAs"
                
            full_path = self.model.get_full_path(
                settings['file_path'],
                settings['file_name'],
                settings['delimiter']
            )

            self._write_export_file(full_path, settings)

            show_message(
                "Export Complete",
                f"Export to {full_path} was successful."
            )
            
            self._handle_cancel()
            
        except PermissionError:
            show_error(self.settings, "Permission Error", 
                      "Cannot access the file. Please ensure it is not open elsewhere.")
        except Exception as e:
            show_error(self.settings, "Export Error", str(e))

    def _write_export_file(self, full_path: str, settings: dict) -> None:
        with open(full_path, 'w') as output_file:
            delimiter = "\t" if settings['delimiter'] == r"\t" else settings['delimiter']
            headers = self.model.get_headers()
            output_file.write(delimiter.join(headers) + "\n")
            self._write_data_rows(output_file, headers, settings)

    def _write_data_rows(self, output_file, headers: list, settings: dict) -> None:
        """Write data rows to the output file"""
        try:
            delimiter = "\t" if settings['delimiter'] == r"\t" else settings['delimiter']
            
            for item in self.model.data['selected_items']:
                row_data = []
                
                if self.model.data['window_type'] == "Multitargeting":
                    for header in headers:
                        if header == "Seed":
                            row_data.append(str(item['seed']))
                        elif header == "Total Repeats":
                            row_data.append(str(item['total_repeats']))
                        elif header == "Avg. Repeats/Scaffold":
                            row_data.append(str(item['avg_repeats_or_scaffold']))
                        elif header == "Consensus Sequence":
                            row_data.append(str(item['consensus_sequence']))
                        elif header == "Full Sequence":
                            sequence = str(item['consensus_sequence'])
                            full_sequence = (settings['leading_sequence'] + 
                                          sequence + 
                                          settings['trailing_sequence'])
                            row_data.append(full_sequence)
                        elif header == "% Consensus":
                            row_data.append(str(item['percent_consensus']))
                        elif header == "Score":
                            row_data.append(str(item['score']))
                        elif header == "PAM":
                            row_data.append(str(item['pam']))
                        elif header == "Strand":
                            row_data.append(str(item['strand']))
                        else:
                            row_data.append("")
                        
                elif self.model.data['window_type'] == "Population Analysis":
                    # Handle Population Analysis specific data format
                    for header in headers:
                        if header == "Seed":
                            row_data.append(str(item['seed']))
                        elif header == "% Coverage":
                            row_data.append(str(item['percent_coverage']))
                        elif header == "Total Repeats":
                            row_data.append(str(item['total_repeats']))
                        elif header == "Avg. Repeats/Scaffold":
                            row_data.append(str(item['avg_repeats_or_scaffold']))
                        elif header == "Consensus Sequence":
                            row_data.append(str(item['consensus_sequence']))
                        elif header == "Full Sequence":
                            sequence = str(item['consensus_sequence'])
                            full_sequence = (settings['leading_sequence'] + 
                                          sequence + 
                                          settings['trailing_sequence'])
                            row_data.append(full_sequence)
                        elif header == "% Consensus":
                            row_data.append(str(item['percent_consensus']))
                        elif header == "Score":
                            row_data.append(str(item['score']))
                        elif header == "PAM":
                            row_data.append(str(item['pam']))
                        elif header == "Strand":
                            row_data.append(str(item['strand']))
                        else:
                            row_data.append("-")
                else:
                    # Handle View Targets window
                    for header in headers:
                        if header == "Location":
                            row_data.append(str(item['location']))
                        elif header == "Endonuclease":
                            row_data.append(str(item['endonuclease']))
                        elif header == "Sequence":
                            row_data.append(str(item['sequence']))
                        elif header == "Full Sequence":
                            sequence = str(item['sequence'])
                            full_sequence = (settings['leading_sequence'] + 
                                          sequence + 
                                          settings['trailing_sequence'])
                            row_data.append(full_sequence)
                        elif header == "Strand":
                            row_data.append(str(item['strand']))
                        elif header == "PAM":
                            row_data.append(str(item['pam']))
                        elif header == "Score":
                            row_data.append(str(item['score']))
                        elif header == "Off-Target":
                            row_data.append(str(item.get('off_target', '--.--')))
                        elif header == "Locus_Tag":
                            row_data.append(str(item.get('locus_tag', '')))
                        elif header == "Gene_Name":
                            row_data.append(str(item.get('gene_name', '')))
                        else:
                            row_data.append("") 
                
                output_file.write(delimiter.join(row_data) + "\n")
                
        except Exception as e:
            self.logger.error(f"Error writing data rows: {str(e)}")
            self.logger.error(f"Data item causing error: {item}")
            self.logger.error(f"Headers: {headers}")
            raise

    def _handle_end_of_row(self, tmp_list: list, output_file, delimiter: str) -> None:
        """Handle end of row processing"""
        tmp_list.append(self.model.data['selected_items'][-1].text())
        if self.model.data['has_locus_tag']:
            tmp = self.settings.mainWindow.Results.comboBoxGene.currentText().split(":")
            tmp_list.extend([tmp[0].strip(), tmp[-1].strip()])
        elif self.model.data['has_gene_name']:
            tmp_list.append(self.settings.mainWindow.Results.comboBoxGene.currentText().strip())
        
        delimiter = "\t" if delimiter == r"\t" else delimiter
        output_file.write(delimiter.join(tmp_list) + "\n")

    def _handle_sequence_column(self, tmp_list: list, item, settings: dict) -> None:
        sequence = item.get('sequence', '') if isinstance(item, dict) else item.text()
        tmp_list.append(sequence)
        full_sequence = (settings['leading_sequence'] + 
                        sequence + 
                        settings['trailing_sequence'])
        tmp_list.append(full_sequence)

    def _handle_cancel(self) -> None:
        self.view.hide()
