import sys, os
import GlobalSettings
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QDir, QUrl
from Bio import Entrez, SeqIO
from bioservices.kegg import KEGG
import platform

class genomebrowser():
	def __init__(self, parent=None):
		self.loadingWindow = loading_window()

	def splitStringNCBI(self, longString):
		try:
			return (longString).split(':')[2]
		except:
			print(longString.split(':'))
			#print("Please search for a genome")


	def splitStringLocal(self, longString):
		try:
			#print(longString)
			#print(longString.split('/'))
			return (longString.split('/').pop()).split('.')[0]
		except:
			#print(longString.split('/').pop()).split('.')[0]
			print("Please search for a genome")


	def ncbiAPI(self, filename):

		gb_file = filename
		gb_file = GlobalSettings.filedir + "/" + gb_file
		if platform.system() == 'Windows':
			gb_file = str(gb_file).replace("/","\\")
		print(gb_file)
		genomeList = []
		for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
			genomeList.append(gb_record.id)

		return genomeList


	def createHtml(self, genomeList):
		print(genomeList)
		htmlString1 = """
		<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
		<ihtml>
		<head>
			<title>NCBI Sequence Viewer - programmatic initialization</title>
			<script type="text/javascript" src="https://www.ncbi.nlm.nih.gov/projects/sviewer/js/sviewer.js"></script>
			<script>
				function loadSV(id) {
					var svapp = SeqView.App.findAppByDivId("mySeqViewer1");
					if (!svapp)
						svapp = new SeqView.App("mySeqViewer1");
					params = 'appname=UniversityofTennesseeCasper&amp;id=' + id;
					svapp.reload(params);
				}
			</script>
		</head>
		<body>
			Select a chromosome from the list:<br/>
	
		<select onchange="loadSV(event.target.value);">
			<option value="">-</option>
		"""

		htmlString2 = """
			<option value="{}">Chromosome {}</option>
		"""

		htmlString3 = """
			</select>
			<br/>
		
			<div id="mySeqViewer1" class="SeqViewerApp" data-autoload> 
				<a href="?embedded=true&appname=testapp1&id={}"></a> 
			</div>	
			
		</body>
		</html>
		""".format(genomeList[0])

		# Find file path for template
		# seek to beginning and truncate
		genomeBrowserTemplateFilePath = GlobalSettings.appdir + "genomeBrowserTemplate.html"
		raw = open(genomeBrowserTemplateFilePath, "r+")
		raw.seek(0)
		raw.truncate()

		#write the 3 part string format
		raw.write(htmlString1)

		for index,genome in enumerate(genomeList):
			raw.write(htmlString2.format(genome,index+1))

		raw.write(htmlString3)


	def createGraph(self, p):
		#launch loading window
		self.loadingWindow.show()
		QtCore.QCoreApplication.processEvents()

		annotationWindow = p.findChild(QtWidgets.QComboBox,'Annotations_Organism')
		NCBIFileBoolean = p.findChild(QtWidgets.QRadioButton, 'NCBI_Select')
		localFileBoolean = p.findChild(QtWidgets.QRadioButton, 'Annotation_Ownfile')

		selectedGenome = p.annotation_files.currentText()

		#returns gci String value that can be used to retrieve the file
		#print("selected Genome is ", selectedGenome)
		gciVariable = ""

		# if(localFileBoolean.isChecked() == True):
		# 	print("local")
		gciVariable = self.splitStringLocal(selectedGenome)
		# elif(NCBIFileBoolean.isChecked() == True):
		# 	print("NCBI")
		# 	gciVariable = splitStringNCBI(selectedGenome)

		#print("gci variable is ", gciVariable)

		if(gciVariable == None):
			print("Please select NCBI Genome")
			return

		directory = GlobalSettings.CSPR_DB
		#print('*****')
		#print(directory)
		# onlyfiles = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
		# for file in onlyfiles:
		# 	print("One file is ",file)
		# 	if gciVariable in file and ".fna" not in file:
		# 		fileToSearch = file
		# 		print("searching file:", fileToSearch)
		# 		break

		fileToSearch = GlobalSettings.mainWindow.annotation_files.currentText()
		#print(fileToSearch)
		if(str(fileToSearch).find(".gbff") == -1):
			#print("filetype not valid")
			QtWidgets.QMessageBox.information(p, "Genomebrowser Error", "Filetype must be GBFF.", QtWidgets.QMessageBox.Ok)
			return


		try:
			#print("file to search is ", fileToSearch)
			genomeList = self.ncbiAPI(fileToSearch)

		except:
			#print("No gbff file found")
			QtWidgets.QMessageBox.question(p, "GBFF_FileNotFound", "GBFF file is not in selected directory", QtWidgets.QMessageBox.Ok)
			return

		self.createHtml(genomeList)
		self.browser = QWebEngineView()
		file_path = GlobalSettings.appdir + "genomeBrowserTemplate.html"
		local_url = QUrl.fromLocalFile(file_path)
		self.browser.load(local_url)
		self.browser.show()
		self.loadingWindow.hide()



#progress bar gui
class loading_window(QtWidgets.QWidget):
	def __init__(self):
		super(loading_window, self).__init__()
		uic.loadUi(GlobalSettings.appdir + "loading_data_form.ui", self)
		self.loading_bar.hide()
		self.setWindowTitle("Loading Genome Browser")
		self.info_label.setText("Loading Genome Browser")
		self.hide()