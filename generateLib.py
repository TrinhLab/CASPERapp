import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore
from CSPRparser import CSPRparser

########################################################################################################################
# Class Name: genLibrary
# this class is a window that allows the user to select the settings for Generate Library
# When the user clicks Generate Library, it goes ahead and gets the Annotation Data needed
#   Then the user can select the settings they want, and then hit submit.
#   It creates a txt file with the data
########################################################################################################################


class genLibrary(QtWidgets.QDialog):
    def __init__(self):
        # qt stuff
        super(genLibrary, self).__init__()
        uic.loadUi('library_prompt.ui', self)
        self.setWindowTitle('Generate Library')
        self.setWindowIcon(Qt.QIcon('cas9image.png'))

        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.BrowseButton.clicked.connect(self.browse_function)

        # variables
        self.anno_data = dict()
        self.cspr_file = ''
        self.parser = CSPRparser('')
        self.kegg_nonKegg = ''

        # set the numbers for the num genes combo box item
        for i in range(10):
            self.numGenescomboBox.addItem(str(i + 1))

        # set the numbers for the minOn combo box
        for i in range(19, 50):
            self.minON_comboBox.addItem(str(i + 1))

    # this function launches the window
    # Parameters:
    #       annotation_data: a dictionary that has the data for the annotations searched for
    #           currently MainWindow's searches dict is passed into this
    #       org_file: the cspr_file that pertains to the organism that user is using at the time
    #       anno_type: whether the user is using KEGG or another type of annotation file
    def launch(self, annotation_data, org_file, anno_type):

        self.cspr_file = org_file
        self.anno_data = annotation_data
        self.kegg_nonKegg = anno_type
        self.parser.fileName = self.cspr_file

        # setting the path and file name fields
        index = self.cspr_file.find('.')
        tempFileName = self.cspr_file[:index] + '_lib.txt'
        self.filename_input.setText(self.cspr_file[:index] + '_lib.txt')
        self.output_path.setText(GlobalSettings.CSPR_DB + "/")

        for data in self.anno_data:
            print(data)
            for item in self.anno_data[data]:
                print('\t', item)
                for piece in self.anno_data[data][item]:
                    print('\t\t', piece)


        print(self.kegg_nonKegg)
        self.show()

    # cancel function
    # clears everything and hides the window
    def cancel_function(self):
        self.cspr_file = ''
        self.anno_data = dict()
        self.kegg_nonKegg = ''

        self.filename_input.setText('')
        self.output_path.setText('')

        self.hide()


    # browse function
    # allows the user to browse for a folder
    # stores their selection in the output_path line edit
    def browse_function(self):
        # get the folder
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                       GlobalSettings.CSPR_DB, QtWidgets.QFileDialog.ShowDirsOnly)
        if(os.path.isdir(mydir) == False):
            return

        # make sure to append the '/' to the folder path
        self.output_path.setText(mydir + "/")
