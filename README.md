# CASPERapp

Thank you for your interest in CASPER.  Please visit the Release page to download the CASPERapp executable for your OS; otherwise, you can download the source code, and run the CASPERapp using Python3.

### To run CASPERapp using python3:
1) Clone the repository
2) Ensure that you have Python3 installed on your computer (Python 3.6 recommended).
3) We have coded in a requirements file such that dependencies CASPER operates on should be automatically downloaded to your computer.
If this does not work and you have a dependency issue, please use pip or anaconda to install the appropriate packages given in
the warnings.
4) Download the executables for OffTarget and SeqFinder for your operating system from the Google Drive links below (note: SeqFinder is not available currently):
  
    Windows: https://drive.google.com/drive/folders/1n45GwEN4glZ5M3apBDNXlPLvw3q5-7YT?usp=sharing <br />
    Mac: https://drive.google.com/drive/folders/1EuX9KnQEW2AHQSty9fzZYfYwYbkNmdiE?usp=sharing <br />
    Linux: https://drive.google.com/drive/folders/1y3z-ty_C7YnpY73iw4zYdYFvvZvJNYDO?usp=sharing

5) Place the OffTarget executable in the OffTargetFolder folder within the CASPERapp directory.
6) Place the SeqFinder executable in the CASPERapp directory (does not need to go into a specific folder).
6) Run CASPER!  Simply cd into the CASPERapp-master directory (wherever you put it) and run the command: "python3 main.py"
CASPER will launch and you will be good to go!

Any problems please email Brian Mendoza at bmendoz1@vols.utk.edu

### Python Dependencies Needed
- PyQt5
- Bioservices
- biopython
- pyqtgraph
- PyQtChart
- matplotlib
- matplotlib-venn
- gffutils
- webbrowser
- gzip
- sqlite3

### Sample CSPR Files to use with the CASPERapp
- At the Google Drive link below, simply download the entire folder, extract all of the files to your CASPERapp database folder on your computer. Then when you launch
CASPERapp, provide the appropriate database folder path.
https://drive.google.com/drive/folders/16HPXOM0k1wfNn1gnbFcMSOFQZHBghGty?usp=sharing
