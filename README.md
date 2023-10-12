# CASPERapp

![](CASPER-logo.jpg)

Software for Enhanced guide-RNA Design and Targeting Analysis for Precise CRISPR Genome Editing of Single and Consortia of Industrially Relevant and Non-model Organisms. This software is a GUI interface for the following paper:

Brian Mendoza and Cong T. Trinh. 2017. Enhanced Guide-RNA Design And Targeting Analysis For Precise CRISPR Genome Editing Of Single And Consortia Of Industrially Relevant And Non-Model Organisms. DOI: https://doi.org/10.1093/bioinformatics/btx564

Thank you for your interest in CASPER. Our packaged releases for Windows 10 and MacOS (Ventura 13.3,1a or greater) are located [here](https://github.com/TrinhLab/CASPERapp/releases). Backwards compatability is not guaranteed. If a release does not work for you, CASPER may also be ran on any OS using a terminal and python3. Please refer to the instructions below.

### How to run CASPERapp using python3:
1) Clone the repository using the command: ```git clone https://github.com/TrinhLab/CASPERapp```
2) Ensure that you have Python3.8 installed on your computer (Python 3.8.13 recommended).
3) cd into the CASPERapp directory: ```cd CASPERapp```
4) _Recommended: Create and enter virtual environment:_ ```virtualenv <env_name> && source <env_name>/bin/activate```
6) Install setuptools with pip (important to do this before installing remaining dependencies): ```pip install setuptools==58.0.0```
7) Install dependencies: ```pip install -r requirements.txt```
8) Run CASPER! ```python3 main.py```

CASPER will launch and you will be good to go! If you have any problems, please email David Dooley at ddooley2@vols.utk.edu

NOTE: CASPER may take a long time to launch for the first time due to initialization in the background. After the first launch, it should load much faster.
