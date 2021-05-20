import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore


########################################################################################################################
# Class: export_csv_window
# This class opens a window for the user to select where they want the CSV file exported to, and the name of the file
# It takes the highlighted data from the Results page, and creates a CSV file from that
########################################################################################################################
class export_csv_window(QtWidgets.QDialog):
    # init function. Sets all of the buttons
    def __init__(self):
        # qt stuff
        super(export_csv_window, self).__init__()
        uic.loadUi(GlobalSettings.appdir + 'export_to_csv_window.ui', self)
        self.setWindowTitle("Export Guides")
        self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + "cas9image.png"))

        # button connections
        self.browse_button.clicked.connect(self.browseForFolder)
        self.cancel_button.clicked.connect(self.cancel_function)
        self.export_button.clicked.connect(self.export_function)

        # variables
        self.location = self.fileLocation_line_edit.text()
        self.selected_table_items = []
        self.num_columns = []

    # launch function. Called in Results.
    # parameter expect: a list of the items selected from the window.
    def launch(self, select_items, num_columns):
        self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "/")
        self.selected_table_items = select_items
        self.num_columns = num_columns
        self.show()

    # export function
    # Takes the path and file name and combines them
    # Writes the header line, as well as ever line selected to that file
    # calls the cancel function when it's done
    def export_function(self):
        # get the full path ( path and file name)
        file_name = self.filename_line_edit.text()
        if file_name == "":
            file_name = "exported_gRNAs"
        self.location = self.fileLocation_line_edit.text() + '/'
        full_path = ""
        if '.csv' in file_name:
            full_path = self.location + file_name
        else:
            full_path = self.location + file_name + '.csv'

        # try to do it
        try:
            #write the table headers
            if self.num_columns == 8: ###Change headers for multitargeting table export
                output_data = open(full_path, 'w')
                output_data.write('Seed,Total Repeats,Avg. Repeats/Scaffold,Consensus Sequence,% Consensus,Score,PAM,Strand\n')
            elif self.num_columns == 9:
                output_data = open(full_path, 'w')
                output_data.write('Seed,% Coverage,Total Repeats,Avg. Repeats/Scaffold,Consensus Sequence,% Consensus,Score,PAM,Strand\n')
            else: ###Change headers for view results export
                output_data = open(full_path, 'w')
                output_data.write('Location,Endonuclease,Sequence,Strand,PAM,Score,Off_Target\n')

            # loop through and write the other data
            i = 0
            for item in self.selected_table_items:
                if i == self.num_columns-1:
                    output_data.write(item.text())
                    output_data.write('\n')
                    i = 0
                else:
                    output_data.write(item.text())
                    output_data.write(',')
                    i += 1
        # catch the permission exception
        except PermissionError:
            QtWidgets.QMessageBox.question(self, "File Cannot Open",
                                           "This file cannot be opened. Please make sure that the file is not opened elsewhere and try again.",
                                           QtWidgets.QMessageBox.Ok)
            return
        # catch any other exception
        except Exception as e:
            print(e)
            return

        # close the window
        self.cancel_function()

    # Resets everything to the init funciton
    # then closes the window
    def cancel_function(self):
        self.fileLocation_line_edit.setText(GlobalSettings.CSPR_DB + "/")
        self.filename_line_edit.setText("")
        self.location = ""
        self.hide()

    # browse for folder function
    # allows user to browse for a folder where to store the CSV file
    def browseForFolder(self):
        # get the folder
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                       GlobalSettings.CSPR_DB, QtWidgets.QFileDialog.ShowDirsOnly)
        if(os.path.isdir(mydir) == False):
            return

        # make sure to append the '/' to the folder path
        self.fileLocation_line_edit.setText(mydir + "/")
        self.location = mydir + "/"
