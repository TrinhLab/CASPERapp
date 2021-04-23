import GlobalSettings
import os
from PyQt5 import QtWidgets, Qt, uic, QtCore
from functools import partial
from CSPRparser import CSPRparser
import re
import platform

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
        uic.loadUi(GlobalSettings.appdir + 'library_prompt.ui', self)
        self.setWindowTitle('Generate Library')
        self.setWindowIcon(Qt.QIcon(GlobalSettings.appdir + 'cas9image.png'))

        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(197,224,173);
                        border-radius: 9px;
                        font: 14pt "Arial";
                        font: bold;
                        margin-top: 10px;}"""
        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2").replace("rgb(197,224,173)", "rgb(111,181,110)"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1", "Step3").replace("rgb(197,224,173)", "rgb(77,158,89)"))
        self.Step4.setStyleSheet(groupbox_style.replace("Step1", "Step4").replace("rgb(197,224,173)", "rgb(53,121,93)"))


        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.BrowseButton.clicked.connect(self.browse_function)
        self.submit_button.clicked.connect(self.submit_data)
        self.progressBar.setValue(0)

        # variables
        self.anno_data = dict()
        self.kegg_nonKegg = ''
        self.gen_lib_dict = dict()
        self.cspr_data = dict()
        self.Output = dict()
        self.off_tol = .05
        self.off_max_misMatch = 4
        self.off_target_running = False
        self.parser = CSPRparser("")

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
        self.db_file = org_file[:org_file.find('.')] + '_repeats.db'
        self.anno_data = annotation_data
        self.kegg_nonKegg = anno_type
        self.process = QtCore.QProcess()
        self.parser.fileName = org_file


        # setting the path and file name fields
        index1 = self.cspr_file.find('.')
        if platform.system() == "Windows":
            index2 = self.cspr_file.rfind('\\')
        else:
            index2 = self.cspr_file.rfind('/')

        self.filename_input.setText(self.cspr_file[index2 + 1:index1] + '_lib')


        if platform.system() == "Windows":
            self.output_path.setText(GlobalSettings.CSPR_DB + "\\")
        else:
            self.output_path.setText(GlobalSettings.CSPR_DB + "/")

        # depending on the type of file, build the dictionary accordingly
        self.build_dict_non_kegg()

        # get the gRNA data from the cspr file
        self.cspr_data = self.parser.gen_lib_parser(self.gen_lib_dict, GlobalSettings.mainWindow.endoChoice.currentText())
        self.show()


    # this is here in case the user clicks 'x' instead of cancel. Just calls the cancel function
    def closeEvent(self, event):
        closeWindow = self.cancel_function()

        # if the user is doing OT and does not decide to cancel it ignore the event
        if closeWindow == -2:
            event.ignore()
        else:
            event.accept()


    # this function takes all of the cspr data and compresses it again for off-target usage
    def compress_file_off(self):
        if platform.system() == "Windows":
            file = GlobalSettings.CSPR_DB + "\\off_input.txt"
        else:
            file = GlobalSettings.CSPR_DB + "/off_input.txt"
        f = open(file, 'w')
        for gene in self.cspr_data:
            for j in range(len(self.cspr_data[gene])):
                loc = self.cspr_data[gene][j][0]
                seq = self.cspr_data[gene][j][1]
                pam = self.cspr_data[gene][j][2]
                score = self.cspr_data[gene][j][3]
                strand = self.cspr_data[gene][j][4]
                output = str(loc) + ';' + str(seq) + ';' + str(pam) + ';' + str(score) + ';' + str(strand)
                f.write(output + '\n')
        f.close()


    # this function parses the temp_off file, which holds the off-target analysis results
    # it also updates each target in the cspr_data dictionary to replace the endo with the target's results in off-target
    def parse_off_file(self):
        if platform.system() == "Windows":
            file = GlobalSettings.CSPR_DB + "\\temp_off.txt"
        else:
            file = GlobalSettings.CSPR_DB + "/temp_off.txt"
        f = open(file, "r")
        file_data = f.read().split('\n')
        f.close()
        scoreDict = dict()

        # get the data from the file
        for i in range(len(file_data)):
            if file_data[i] == 'AVG OUTPUT':
                continue
            elif file_data[i] != '':
                buffer = file_data[i].split(':')
                scoreDict[buffer[0]] = buffer[1]

        # update cspr_Data
        for gene in self.cspr_data:
            for i in range(len(self.cspr_data[gene])):
                tempTuple = (self.cspr_data[gene][i][0], self.cspr_data[gene][i][1], self.cspr_data[gene][i][2], self.cspr_data[gene][i][3], self.cspr_data[gene][i][4], scoreDict[self.cspr_data[gene][i][1]])
                self.cspr_data[gene][i] = tempTuple


    # this function runs the off_target command
    # NOTE: some changes may be needed to get it to work with other OS besides windows
    def get_offTarget_data(self, num_targets, minScore, spaceValue, output_file, fiveseq):
        self.perc = False
        self.bool_temp = False
        self.running = False

        # when finished, parse the off file, and then generate the lib
        def finished():
            if self.off_target_running:
                self.progressBar.setValue(100)
                self.parse_off_file()
                did_work = self.generate(num_targets, minScore, spaceValue, output_file, fiveseq)
                self.off_target_running = False
                #self.process.kill()
                if did_work != -1:
                    self.cancel_function()
                    #os.remove(GlobalSettings.CSPR_DB + '/off_input.txt')
                    #os.remove(GlobalSettings.CSPR_DB + '/temp_off.txt')

        # as off-targeting outputs things, update the off-target progress bar
        def progUpdate(p):
            line = str(self.process.readAllStandardOutput())
            line = line[2:]
            line = line[:len(line) - 1]
            if platform.system() == 'Windows':
                for lines in filter(None, line.split(r'\r\n')):
                    if (lines.find("Running Off Target Algorithm for") != -1 and self.perc == False):
                        self.perc = True
                    if (self.perc == True and self.bool_temp == False and lines.find(
                            "Running Off Target Algorithm for") == -1):
                        lines = lines[32:]
                        lines = lines.replace("%", "")
                        if (float(lines) <= 99.5):
                            num = float(lines)
                            self.progressBar.setValue(num)
                        else:
                            self.bool_temp = True
            else:
                for lines in filter(None, line.split(r'\n')):
                    if (lines.find("Running Off Target Algorithm for") != -1 and self.perc == False):
                        self.perc = True
                    if (self.perc == True and self.bool_temp == False and lines.find(
                            "Running Off Target Algorithm for") == -1):
                        lines = lines[32:]
                        lines = lines.replace("%", "")
                        if (float(lines) <= 99.5):
                            num = float(lines)
                            self.progressBar.setValue(num)
                        else:
                            self.bool_temp = True

        app_path = GlobalSettings.appdir.replace('\\','/')
        if platform.system() == 'Windows':
            exe_path = app_path + r'OffTargetFolder/OT_Win.exe'
        elif platform.system() == 'Linux':
            exe_path = app_path + r'OffTargetFolder/OT_Lin' 
        else:
            exe_path = app_path + r'OffTargetFolder/OT_Mac' 
        exe_path = '"' + exe_path + '" '
        data_path = '"' + GlobalSettings.CSPR_DB + "/off_input.txt" + '" '
        cspr_path = '"' + self.cspr_file + '" '
        db_path = '"' + self.db_file + '" '
        output_path = '"' + GlobalSettings.CSPR_DB + '/temp_off.txt" '
        filename = output_path
        filename = filename[:len(filename) - 1]
        filename = filename[1:]
        filename = filename.replace('"', '')
        CASPER_info_path = '"' + app_path + 'CASPERinfo' +'" '
        num_of_mismathes = self.off_max_misMatch
        tolerance = self.off_tol  # create command string

        detailed_output = " False "
        avg_output = "True"
        # set the off_target_running to true, to keep the user from closing the window while it is running
        self.off_target_running = True

        cmd = exe_path + data_path + cspr_path + db_path + output_path + CASPER_info_path + str(
            num_of_mismathes) + ' ' + str(tolerance) + detailed_output + avg_output
        
        print(cmd)

        if platform.system() == 'Windows':
            cmd = cmd.replace('/', '\\')
        self.process.readyReadStandardOutput.connect(partial(progUpdate, self.process))
        self.process.readyReadStandardError.connect(partial(progUpdate, self.process))
        self.progressBar.setValue(0)
        QtCore.QTimer.singleShot(100, partial(self.process.start, cmd))
        self.process.finished.connect(finished)


    # submit function
    # this function takes all of the input from the window, and calls the generate function
    # Still need to add the checks for 5' seq, and the percentage thing
    def submit_data(self):
        if self.off_target_running:
            return
        output_file = self.output_path.text() + self.filename_input.text()

        minScore = int(self.minON_comboBox.currentText())
        num_targets = int(self.numGenescomboBox.currentText())
        fiveseq = ''

        # error check for csv files
        if output_file.endswith('.txt'):
            output_file = output_file.replace('.txt', '.csv')
        elif not output_file.endswith('.txt') and not output_file.endswith('.csv'):
            output_file = output_file + '.csv'

        # error checking for the space value
        # if they enter nothing, default to 15 and also make sure it's actually a digit
        if self.space_line_edit.text() == '':
            spaceValue = 15
        elif self.space_line_edit.text().isdigit():
            spaceValue = int(self.space_line_edit.text())
        elif not self.space_line_edit.text().isdigit():
            QtWidgets.QMessageBox.question(self, "Error", "Please enter integers only for space between guides.",
                                           QtWidgets.QMessageBox.Ok)
            return
        # if space value is more than 200, default to 200
        if spaceValue > 200:
            spaceValue = 200
        elif spaceValue < 0:
            QtWidgets.QMessageBox.question(self, "Error", "Please enter a space-value that is 0 or greater.",
                                            QtWidgets.QMessageBox.Ok)
            return

        if self.find_off_Checkbox.isChecked():
            self.compress_file_off()



        # get the fiveprimseq data and error check it
        if self.fiveprimeseq.text() != '' and self.fiveprimeseq.text().isalpha():
            fiveseq = self.fiveprimeseq.text()
        elif self.fiveprimeseq.text() != '' and not self.fiveprimeseq.text().isalpha():
            QtWidgets.QMessageBox.question(self, "Error", "Please make sure only the letters A, T, G, or C are added into 5' End specificity box.",
                                           QtWidgets.QMessageBox.Ok)
            return


        # get the targeting range data, and error check it here
        if not self.start_target_range.text().isdigit() or not self.end_target_range.text().isdigit():
            QtWidgets.QMessageBox.question(self, "Error",
                                           "Error: Please make sure that the start and end target ranges are numbers only."
                                           " Please make sure that start is 0 or greater, and end is 100 or less. ",
                                           QtWidgets.QMessageBox.Ok)
            return
        elif int(self.start_target_range.text()) >= int(self.end_target_range.text()):
            QtWidgets.QMessageBox.question(self, "Error",
                                           "Please make sure that the start number is always less than the end number",
                                           QtWidgets.QMessageBox.Ok)
            return


        # if they check Off-Targeting
        if self.find_off_Checkbox.isChecked():
            # make sure its a digit
            if self.maxOFF_comboBox.text() == '' or not self.maxOFF_comboBox.text().isdigit() and '.' not in self.maxOFF_comboBox.text():
                QtWidgets.QMessageBox.question(self, "Error", "Please enter only numbers for Maximum Off-Target Score. It cannot be left blank",
                                               QtWidgets.QMessageBox.Ok)
                return
            else:
                # make sure it between 0 and .5
                if not 0.0 < float(self.maxOFF_comboBox.text()) <= .5:
                    QtWidgets.QMessageBox.question(self, "Error",
                                                   "Please enter a max off-target score between 0 and 0.5!",
                                                   QtWidgets.QMessageBox.Ok)
                    return
                # compress the data, and then run off-targeting
                self.compress_file_off()
                self.get_offTarget_data(num_targets, minScore, spaceValue, output_file, fiveseq)
        else:
        # actually call the generate function
            did_work = self.generate(num_targets, minScore, spaceValue, output_file, fiveseq)

            if did_work != -1:
                self.cancel_function()


    # cancel function
    # clears everything and hides the window
    def cancel_function(self):
        if self.off_target_running:
            error = QtWidgets.QMessageBox.question(self, "Off-Targeting is running",
                                            "Off-Targetting is running. Closing this window will cancel that process, and return to the main window. .\n\n"
                                            "Do you wish to continue?",
                                            QtWidgets.QMessageBox.Yes |
                                            QtWidgets.QMessageBox.No,
                                            QtWidgets.QMessageBox.No)
            if (error == QtWidgets.QMessageBox.No):
                    return -2
            else:
                self.off_target_running = False
                self.process.kill()

        QtWidgets.QMessageBox.information(self, "Library Generated!",
                                   "CASPER has finished generating your library!",
                                   QtWidgets.QMessageBox.Ok)

        self.cspr_file = ''
        self.anno_data = dict()

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
        self.off_target_running = False
        self.progressBar.setValue(0)

        self.hide()

    # browse function
    # allows the user to browse for a folder
    # stores their selection in the output_path line edit
    def browse_function(self):
        if self.off_target_running:
            return
        # get the folder
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                       GlobalSettings.CSPR_DB, QtWidgets.QFileDialog.ShowDirsOnly)
        if(os.path.isdir(mydir) == False):
            return

        # make sure to append the '/' to the folder path
        self.output_path.setText(mydir + "/")


    # this function builds the dictionary that is used in the generate function
    # this is the version that builds it from data from feature_table, gbff, or gff
    # builds it exactly as Brian built it in the files given
    def build_dict_non_kegg(self):
        for gene in self.anno_data:
            temp = gene.split(';')
            gene_id = temp[0]
            description = temp[-1]
            gene_name = temp[len(temp)-2]
            ### Order: chromosome number, gene start, gene end, dir of gene, gene description, gene name/locus tag
            self.gen_lib_dict[gene_id] = [self.anno_data[gene][0][1], self.anno_data[gene][0][3], self.anno_data[gene][0][4], self.anno_data[gene][0][5],description,gene_name]


    # generate function taken from Brian's code
    def generate(self,num_targets_per_gene, score_limit, space, output_file, fiveseq):
        deletedDict = dict()

        # check and see if we need to search based on target_range
        startNum = float(self.start_target_range.text())
        endNum = float(self.end_target_range.text())
        checkStartandEndBool = False
        if startNum != 0.0 or endNum != 100.0:
            if startNum >= 0.0 and endNum <= 100.0:
                startNum = startNum / 100
                endNum = endNum / 100
                checkStartandEndBool = True
            else:
                QtWidgets.QMessageBox.question(self, "Invalid Targeting Range:",
                                           "Please select a targeting range between 0 and 100.",
                                           QtWidgets.QMessageBox.Ok)
                return -1

        for gene in self.gen_lib_dict:
            target_list = self.cspr_data[gene]  # Gets the gRNAs for given gene 

            #target_list = chrom_list[k:l+1]
            # Reverse the target list if the gene is on negative strand:
            if self.gen_lib_dict[gene][3] == "-":
                target_list.reverse()

            # Filter out the guides with low scores and long strings of T's
            # also store the ones deleted if the user selects 'modify search parameters'
            if self.modifyParamscheckBox.isChecked():
                deletedDict[gene] = list()
            for i in range(len(target_list) - 1, -1, -1): ### Start at end and move backwards through list
                # check the target_range here
                if int(target_list[i][3]) < int(score_limit):
                    if self.modifyParamscheckBox.isChecked():
                        deletedDict[gene].append(target_list[i])
                    target_list.pop(i)
                # check for gRNAs with poly T regions here
                elif re.search("T{5,10}", target_list[i][1]) is not None:
                    if self.modifyParamscheckBox.isChecked():
                        deletedDict[gene].append(target_list[i])
                    target_list.pop(i)

            # check for the fiveseq
            if fiveseq != '':
                for i in range(len(target_list) - 1, -1, -1): ### Start at end and move backwards through list
                    if not target_list[i][1].startswith(fiveseq.upper()):
                        if self.modifyParamscheckBox.isChecked():
                            deletedDict[gene].append(target_list[i])
                        target_list.pop(i)
            # check the target range here
            if checkStartandEndBool:
                for i in range(len(target_list) - 1, -1, -1):
                    totalDistance = self.gen_lib_dict[gene][2] - self.gen_lib_dict[gene][1]
                    target_loc = abs(int(target_list[i][0])) - int(self.gen_lib_dict[gene][1])
                    myRatio = target_loc / totalDistance

                    if not (startNum <= myRatio <= endNum):
                        if self.modifyParamscheckBox.isChecked():
                            deletedDict[gene].append(target_list[i])
                        target_list.pop(i)
            # if the user selected off-targeting, check to see that the targets do not exceed the selected max score
            if self.find_off_Checkbox.isChecked():
                maxScore = float(self.maxOFF_comboBox.text())
                for i in range(len(target_list) - 1, -1, -1):
                    if maxScore < float(target_list[i][5]):
                        if self.modifyParamscheckBox.isChecked():
                            deletedDict[gene].append(target_list[i])
                        target_list.pop(i)
            # Now generating the targets
            self.Output[gene] = list()
            i = 0
            vec_index = 0
            prev_target = (0, "xyz", 'abc', 1, "-")
            while i < num_targets_per_gene:
                # select the first five targets with the score and space filter that is set in the beginning
                if len(target_list) == 0 or vec_index >= len(target_list):
                    break
                while abs(int(target_list[vec_index][0]) - int(prev_target[0])) < int(space):
                    if target_list[vec_index][3] > prev_target[3] and prev_target != (0,"xyz", "abc", 1, "-"):
                        self.Output[gene].remove(prev_target)
                        self.Output[gene].append(target_list[vec_index])
                        prev_target = target_list[vec_index]
                    vec_index += 1
                    # check and see if there will be a indexing error
                    if vec_index >= len(target_list) - 1:
                        vec_index = vec_index - 1
                        break
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
                            strand = deletedDict[gene][i][4] + '*'
                            endo = deletedDict[gene][i][5]
                            self.Output[gene].append((loc, seq, pam, score, strand, endo))

        # Now output to the file
        try:
            f = open(output_file, 'w')
            # if OT checked
            if self.find_off_Checkbox.isChecked():
                f.write('Gene Name,Sequence,On-Target Score,Off-Target Score,Location,PAM,Strand\n')
            elif not self.find_off_Checkbox.isChecked():
                f.write('Gene Name,Sequence,On-Target Score,Location,PAM,Strand\n')
            
            for gene in self.Output:
                i = 0
                gene_name = self.gen_lib_dict[gene][-1]
                for target in self.Output[gene]:
                    # check to see if the target did not match the user's parameters and they selected 'modify'
                    # if the target has an error, put 2 asterisks in front of the target sequence
                    if '*' in target[4]:
                        tag_id = "**" + gene_name + "-" + str(i + 1)
                    else:
                        tag_id = gene_name + "-" + str(i + 1)
                    i += 1

                    tag_id = tag_id.replace(',', '')

                    # if OT checked
                    if self.find_off_Checkbox.isChecked():
                        f.write(tag_id + ',' + target[1] + ',' + str(target[3]) + ',' + str(target[5]) + ',' + str(target[0]) + ',' + target[2] + ',' + target[4][0] + '\n')
                    # if OT not checked
                    elif not self.find_off_Checkbox.isChecked():
                        f.write(tag_id + ',' + target[1] + ',' + str(target[3]) + ',' + str(target[0]) + ',' + target[2] + ',' + target[4][0] + '\n')

            f.close()
        except PermissionError:
            QtWidgets.QMessageBox.question(self, "File Cannot Open",
                                           "This file cannot be opened. Please make sure that the file is not opened elsewhere and try again.",
                                           QtWidgets.QMessageBox.Ok)
            return -1
        except Exception as e:
            print(e)
            return