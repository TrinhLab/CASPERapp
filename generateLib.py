import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore
from Algorithms import SeqTranslate
from functools import partial
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
        self.cspr_data = dict()
        self.Output = dict()
        self.off_tol = .05
        self.off_max_misMatch = 4
        self.process = QtCore.QProcess()


        # set the numbers for the num genes combo box item
        for i in range(10):
            self.numGenescomboBox.addItem(str(i + 1))

        # set the numbers for the minOn combo box
        for i in range(19, 70):
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
        index1 = self.cspr_file.find('.')
        index2 = self.cspr_file.rfind('/')
        self.filename_input.setText(self.cspr_file[index2 + 1:index1] + '_lib.txt')
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
        self.cspr_data = self.parser.gen_lib_parser(self.gen_lib_dict, GlobalSettings.mainWindow.endoChoice.currentText())
        #self.generate(5, 200000000000, 15, "mybsulibrary2.txt")

        #for i in range(len(self.cspr_data)):
         #   for j in range(len(self.cspr_data[i])):
          #      print(self.cspr_data[i][j])
           # print('\n\n')

        self.show()


    # this is here in case the user clicks 'x' instead of cancel. Just calls the cancel function
    def closeEvent(self, event):
        self.cancel_function()
        event.accept()


    # code that is block commented out will be used for the Off_target
    """
    def compress_file_off(self):
        f = open(GlobalSettings.appdir + "/off_compressed.txt", 'w')

        for i in range(len(self.cspr_data)):
            for j in range(len(self.cspr_data[i])):
                loc = self.S.compress(self.cspr_data[i][j][0], 64)
                seq = self.S.compress(self.cspr_data[i][j][1], 64)
                pam = self.S.compress(self.cspr_data[i][j][2], 64)
                score = self.S.compress(self.cspr_data[i][j][3], 64)
                strand = self.S.compress(self.cspr_data[i][j][4], 64)

                output = str(loc) + ',' + str(seq) + str(strand) + str(pam) + ',' + score
                f.write(output + '\n')
        f.close()
    """
    """
    def get_offTarget_data(self):
        def parse_off_data():
        app_path = GlobalSettings.appdir
        exe_path = app_path + '/OffTargetFolder/CasperOffTargetWindows'
        exe_path = '"' + exe_path + '" '
        data_path = '"' + GlobalSettings.appdir + "/off_compressed.txt" + '" '
        compressed = r' True '  ##
        cspr_path = '"' + self.cspr_file + '" '
        output_path = '"' + GlobalSettings.appdir + '/temp_off.txt" '
        filename = output_path
        filename = filename[:len(filename) - 1]
        filename = filename[1:]
        filename = filename.replace('"', '')
        CASPER_info_path = r' "' + app_path + '\\CASPERinfo' + '" '
        num_of_mismathes = self.off_max_misMatch
        tolerance = self.off_tol  # create command string

        detailed_output = " False "
        avg_output = "True"

        cmd = exe_path + data_path + compressed + cspr_path + output_path + CASPER_info_path + str(
            num_of_mismathes) + ' ' + str(tolerance) + detailed_output + avg_output

        print(cmd)
        QtCore.QTimer.singleShot(100, partial(self.process.start, cmd))
        self.process.finished.connect(parse_off_data)
    """

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

        """
        if self.find_off_Checkbox.isChecked():
            if self.maxOFF_comboBox.text() == '' or not self.maxOFF_comboBox.text().isdigit():
                QtWidgets.QMessageBox.question(self, "Error", "Please enter only numbers for Maximum Off-Target Score. It cannot be left blank",
                                               QtWidgets.QMessageBox.Ok)
                return
            else:
                self.compress_file_off()
                self.get_offTarget_data()
        """

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
        elif self.fiveprimeseq.text() != '' and not self.fiveprimeseq.text().isalpha():
            QtWidgets.QMessageBox.question(self, "Error", "Please make sure only the letters A, T, G, or C are added into 5' End specificity box.",
                                           QtWidgets.QMessageBox.Ok)
            return

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
        deletedDict = dict()

        index = 0
        for gene in self.gen_lib_dict:
            target_list = self.cspr_data[gene]  # Gets the chromosome the gene is on

            index += 1

            """
            j = 0 # j is the location
            k = 0 #  This keeps track of the index of the chrom list to start at
            l = 0 # This keeps track of the index of the chrom list to end at
            # this loop sets j and k to the be indeces of the start and stop targets
            print(self.gen_lib_dict[gene])
            while j < self.gen_lib_dict[gene][1]:
                j = chrom_list[k][0]  # k is the index of the item, 0 is the location
                k += 1
                l = k
            while j < self.gen_lib_dict[gene][2] and l != len(chrom_list) - 1:
                print('in second while loop')
                j = chrom_list[l][0]
                l += 1
            """
            #target_list = chrom_list[k:l+1]
            # Reverse the target list if the gene is on negative strand:
            if self.gen_lib_dict[gene][3] == "-":
                target_list.reverse()


            # Filter out the guides with low scores and long strings of T's
            # Also check for the fiveseq if selected
            # also store the ones deleted if the user selects 'modify search parameters'
            if self.modifyParamscheckBox.isChecked():
                deletedDict[gene] = list()
            for i in range(len(target_list) - 1, -1, -1):

                # check for the fiveseq one here
                if fiveseq != '':
                    if not target_list[i][1].startswith(fiveseq.upper()):
                        if self.modifyParamscheckBox.isChecked():
                            deletedDict[gene].append(target_list[i])
                        target_list.pop(i)
                elif target_list[i][3] < score_limit:
                    if self.modifyParamscheckBox.isChecked():
                        deletedDict[gene].append(target_list[i])
                    target_list.pop(i)
                    # del target_list[target_list[i]]
                elif re.search("T{5,10}", target_list[i][1]) is not None:
                    if self.modifyParamscheckBox.isChecked():
                        deletedDict[gene].append(target_list[i])
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

        # if the user selects modify search parameters, go through and check to see if each one has the number of targets that the user wanted
        # if not, append from the deletedDict until they do
        if self.modifyParamscheckBox.isChecked():
            for gene in self.Output:
                if len(self.Output[gene]) < num_targets_per_gene:
                    for i in range(len(deletedDict[gene])):
                        if len(self.Output[gene]) == num_targets_per_gene:
                            break
                        else:
                            loc = deletedDict[gene][i][0]
                            seq = deletedDict[gene][i][1]
                            pam = deletedDict[gene][i][2]
                            score = deletedDict[gene][i][3]
                            strand = deletedDict[gene][i][4]
                            endo = deletedDict[gene][i][5] + '*'
                            self.Output[gene].append((loc, seq, pam, score, strand, endo))

        # Now output to the file
        f = open(output_file, 'w')
        for essential in self.Output:
            i = 0
            for target in self.Output[essential]:
                if '*' in target[5]:
                    tag_id = "**" + essential + "-" + str(i + 1)
                else:
                    tag_id = essential + "-" + str(i + 1)
                i += 1
                f.write(tag_id + "," + target[1] + "\n")
        f.close()








