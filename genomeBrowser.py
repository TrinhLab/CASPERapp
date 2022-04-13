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
import traceback
import webbrowser

#global logger
logger = GlobalSettings.logger

ssl._create_default_https_context = ssl._create_unverified_context

class WebEnginePage(PyQt5.QtWebEngineWidgets.QWebEnginePage):
	def certificateError(self, certificateError):
		print("ssl error")


class genomebrowser(QtWidgets.QWidget):
	def __init__(self, parent=None):
		try:
			default_config = QtNetwork.QSslConfiguration.defaultConfiguration()
			default_config.setProtocol(QtNetwork.QSsl.TlsV1_2)
			QtNetwork.QSslConfiguration.setDefaultConfiguration(default_config)
		except Exception as e:
			logger.critical("Error initializing genomebrowser class.")
			logger.critical(e)
			logger.critical(traceback.format_exc())
			msgBox = QtWidgets.QMessageBox()
			msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
			msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
			msgBox.setWindowTitle("Fatal Error")
			msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
			msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
			msgBox.exec()
			exit(-1)

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
		try:
			if platform.system() == 'Windows':
				filename = str(filename).replace("/","\\")

			genomeList = []
			for gb_record in SeqIO.parse(open(filename, "r"), "genbank"):
				genomeList.append(gb_record.id)

			return genomeList

		except Exception as e:
			logger.critical("Error in ncbiAPI() in genomebrowser.")
			logger.critical(e)
			logger.critical(traceback.format_exc())
			msgBox = QtWidgets.QMessageBox()
			msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
			msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
			msgBox.setWindowTitle("Fatal Error")
			msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
			msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
			msgBox.exec()

			exit(-1)

	def createHtml(self, genomeList):
		try:
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
		except Exception as e:
			logger.critical("Error in createHtml() in genomebrowser.")
			logger.critical(e)
			logger.critical(traceback.format_exc())
			msgBox = QtWidgets.QMessageBox()
			msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
			msgBox.setWindowTitle("Fatal Error")
			msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
			msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
			msgBox.exec()

			exit(-1)

	def createGraph(self, p):
		try:
			selectedGenome = p.annotation_files.currentText()
			print(selectedGenome)
			gciVariable = self.splitStringLocal(selectedGenome)
			print(gciVariable)

			if(gciVariable == None):
				return

			fileToSearch = GlobalSettings.mainWindow.annotation_files.currentText()
			for file in glob.glob(GlobalSettings.CSPR_DB + "/**/*.gb*", recursive=True):
				if file.find(fileToSearch) != -1:
					fileToSearch = file
					break

			if(str(fileToSearch).find(".gb") == -1):
				QtWidgets.QMessageBox.information(p, "Genomebrowser Error", "Filetype must be GenBank format.", QtWidgets.QMessageBox.Ok)
				return


			try:
				genomeList = self.ncbiAPI(fileToSearch)
			except:
				QtWidgets.QMessageBox.question(p, "GenBank_FileNotFound", "GenBank file is not in selected directory", QtWidgets.QMessageBox.Ok)
				return

			self.createHtml(genomeList)
			self.browser = QWebEngineView()

			file_path = GlobalSettings.appdir + "genomeBrowserTemplate.html"
			if platform.system() == "Darwin":
				file_path = "file:///"+file_path
			### Add logic for Linux ?
			webbrowser.open(file_path, new=2)
		except Exception as e:
			logger.critical("Error in createGraph() in genomebrowser.")

			logger.critical(e)
			logger.critical(traceback.format_exc())
			msgBox = QtWidgets.QMessageBox()
			msgBox.setStyleSheet("font: " + str(self.fontSize) + "pt 'Arial'")
			msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
			msgBox.setWindowTitle("Fatal Error")
			msgBox.setText("Fatal Error:\n"+str(e)+ "\n\nFor more information on this error, look at CASPER.log in the application folder.")
			msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Close)
			msgBox.exec()

			exit(-1)