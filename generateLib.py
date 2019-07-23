import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
import re

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
        self.submit_button.clicked.connect(self.submit_data)

        # variables
        self.anno_data = dict()
        self.cspr_file = ''
        self.parser = CSPRparser('')
        self.kegg_nonKegg = ''
        self.gen_lib_dict = dict()
        self.S = SeqTranslate()
        self.cspr_data = list()
        self.Output = dict()

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
        self.filename_input.setText(self.cspr_file[:index] + '_lib.txt')
        self.output_path.setText(GlobalSettings.CSPR_DB + "/")

        # testing:
        #for data in self.anno_data:
         #   print(data)
          #  for item in self.anno_data[data]:
           #     print('\t', item)
            #    for piece in self.anno_data[data][item]:
             #       print('\t\t', piece)
        # print(self.kegg_nonKegg)

        # depending on the type of file, build the dictionary accordingly
        if self.kegg_nonKegg == 'kegg':
            self.build_dict_kegg_version()
        else:
            self.build_dict_non_kegg()

        # get the data from the cspr file
        self.cspr_data = self.parser.gen_lib_parser()
        #self.generate(5, 200000000000, 15, "mybsulibrary2.txt")

        self.show()


    # this is here in case the user clicks 'x' instead of cancel. Just calls the cancel function
    def closeEvent(self, event):
        self.cancel_function()
        event.accept()


    # submit function
    # this function takes all of the input from the window, and calls the generate function
    # Still need to add the checks for 5' seq, and the percentage thing
    def submit_data(self):
        output_file = self.output_path.text() + self.filename_input.text()
        minScore = int(self.minON_comboBox.currentText())
        num_targets = int(self.numGenescomboBox.currentText())
        fiveseq = ''

        # if they enter nothing, default to 15 and also make sure it's actually a digit
        if self.space_line_edit.text() == '':
            spaceValue = 15
        elif self.space_line_edit.text().isdigit():
            spaceValue = int(self.space_line_edit.text())
        elif not self.space_line_edit.text().isdigit():
            QtWidgets.QMessageBox.question(self, "Error", "Please enter integers only for space between guides.",
                                           QtWidgets.QMessageBox.Ok)
            return



        ###
        # need to check that Max off target is numbers only and that space value is numbers only
        # also need to check that targeting range are numbers only, and that the values match the rules
        ###

        # if space value is more than 200, default to 200
        if spaceValue > 200:
            spaceValue = 200
        elif spaceValue < 0:
            QtWidgets.QMessageBox.question(self, "Error", "Please enter a space-value that is 0 or greater.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # get the fiveprimseq data
        if self.fiveprimeseq.text() != '' and self.fiveprimeseq.text().isalpha():
            fiveseq = self.fiveprimeseq.text()

        self.generate(num_targets, minScore, spaceValue, output_file, fiveseq)

        self.cancel_function()

    # cancel function
    # clears everything and hides the window
    def cancel_function(self):
        self.cspr_file = ''
        self.anno_data = dict()
        self.kegg_nonKegg = ''

        self.filename_input.setText('')
        self.output_path.setText('')

        self.gen_lib_dict.clear()
        self.cspr_data.clear()
        self.Output.clear()

        self.start_target_range.setText('0')
        self.end_target_range.setText('100')
        self.space_line_edit.setText('15')
        self.find_off_Checkbox.setChecked(False)
        self.modifyParamscheckBox.setChecked(False)
        self.maxOFF_comboBox.setText('')
        self.fiveprimeseq.setText('')

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


    # this function builds the dictionary that is used in the generate function
    # this is the version that builds it from the KEGG data
    # builds it exactly as Brian built it in the files given
    def build_dict_kegg_version(self):
        for search in self.anno_data:
            for gene in self.anno_data[search]:
                for i in range(len(self.anno_data[search][gene])):
                    self.gen_lib_dict[gene] = [self.anno_data[search][gene][i][0], self.anno_data[search][gene][i][2], self.anno_data[search][gene][i][3], self.anno_data[search][gene][i][1]]

    # this function builds the dictionary that is used in the generate function
    # this is the version that builds it from data from feature_table, gbff, or gff
    # builds it exactly as Brian built it in the files given
    def build_dict_non_kegg(self):
        for search in self.anno_data:
            for gene in self.anno_data[search]:

                descript = gene.split(';')
                temp_descript = descript[0]
                if temp_descript == 'hypothetical protein':
                    temp_descript = temp_descript + " " + str(self.anno_data[search][gene][0][3])

                self.gen_lib_dict[temp_descript] = [self.anno_data[search][gene][0][1], self.anno_data[search][gene][0][3], self.anno_data[search][gene][0][4], self.anno_data[search][gene][0][5]]


    # generate function taken from Brian's code
    def generate(self,num_targets_per_gene, score_limit, space, output_file, fiveseq):
        for gene in self.gen_lib_dict:
            chrom_list = self.cspr_data[self.gen_lib_dict[gene][0] - 1]  # Gets the chromosome the gene is on
            j = 0
            k = 0 #  This keeps track of the index of the chrom list to start at
            l = 0 # This keeps track of the index of the chrom list to end at
            # this loop sets j and k to the be indeces of the start and stop targets
            while j < self.gen_lib_dict[gene][1]:
                j = chrom_list[k][0]  # k is the index of the item, 0 is the location
                k += 1
                l = k
            while j < self.gen_lib_dict[gene][2]:
                j = chrom_list[l][0]
                l += 1
            target_list = chrom_list[k:l+1]
            # Reverse the target list if the gene is on negative strand:
            if self.gen_lib_dict[gene][3] == "-":
                target_list.reverse()
            # Filter out the guides with low scores and long strings of T's

            for i in range(len(target_list) - 1, -1, -1):

                # check for the fiveseq one here
                if fiveseq != '':
                    if not target_list[i][1].startswith(fiveseq.upper()):
                        target_list.pop(i)
                elif target_list[i][3] < score_limit:
                    target_list.pop(i)
                    # del target_list[target_list[i]]
                elif re.search("T{5,10}", target_list[i][1]) is not None:
                    target_list.pop(i)
                    # del target_list[target_list[i]]
                # Now generating the targets

            self.Output[gene] = list()
            i = 0
            vec_index = 0
            prev_target = (0, "xyz", 'abc', 1, "-")
            while i < num_targets_per_gene:
                # select the first five targets with the score and space filter that is set in the beginning
                if len(target_list) == 0 or vec_index >= len(target_list):
                    break
                while abs(target_list[vec_index][0] - prev_target[0]) < space:
                    if target_list[vec_index][3] > prev_target[3]:
                        self.Output[gene].remove(prev_target)
                        self.Output[gene].append(target_list[vec_index])
                        prev_target = target_list[vec_index]
                    vec_index += 1
                # Add the new target to the output and add another to i
                self.Output[gene].append(target_list[vec_index])
                prev_target = target_list[vec_index]
                i += 1
                vec_index += 1
        # Now output to the file
        f = open(output_file, 'w')
        for essential in self.Output:
            i = 0
            for target in self.Output[essential]:
                tag_id = essential + "-" + str(i + 1)
                i += 1
                f.write(tag_id + "," + target[1] + "\n")
        f.close()







