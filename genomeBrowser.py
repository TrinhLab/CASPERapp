import sys, os
import GlobalSettings
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
import PyQt5.QtWebEngineWidgets
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QDir, QUrl
from Bio import Entrez, SeqIO
import PyQt5.QtNetwork as QtNetwork
import platform
import glob

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

class WebEnginePage(PyQt5.QtWebEngineWidgets.QWebEnginePage):
	def certificateError(self, certificateError):
		print("ssl error")

class genomebrowser(QtWidgets.QWidget):
	def __init__(self, parent=None):
		self.loadingWindow = loading_window()
		default_config = QtNetwork.QSslConfiguration.defaultConfiguration()
		default_config.setProtocol(QtNetwork.QSsl.TlsV1_2)
		QtNetwork.QSslConfiguration.setDefaultConfiguration(default_config)

	def splitStringNCBI(self, longString):
		try:
			return (longString).split(':')[2]
		except:
			pass


	def splitStringLocal(self, longString):
		try:
			return (longString.split('/').pop()).split('.')[0]
		except:
			pass


	def ncbiAPI(self, filename):

		if platform.system() == 'Windows':
			filename = str(filename).replace("/","\\")

		genomeList = []
		for gb_record in SeqIO.parse(open(filename, "r"), "genbank"):
			genomeList.append(gb_record.id)

		return genomeList


	def createHtml(self, genomeList):
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
		selectedGenome = p.annotation_files.currentText()
		gciVariable = self.splitStringLocal(selectedGenome)

		if(gciVariable == None):
			return

		fileToSearch = GlobalSettings.mainWindow.annotation_files.currentText()
		for file in glob.glob(GlobalSettings.CSPR_DB + "/**/*.gbff", recursive=True):
			if file.find(fileToSearch) != -1:
				fileToSearch = file
				break

		if(str(fileToSearch).find(".gbff") == -1):
			QtWidgets.QMessageBox.information(p, "Genomebrowser Error", "Filetype must be GBFF.", QtWidgets.QMessageBox.Ok)
			return


		try:
			genomeList = self.ncbiAPI(fileToSearch)
		except:
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