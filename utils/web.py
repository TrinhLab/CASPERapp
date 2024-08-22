import webbrowser
from utils.ui import show_error

def ncbi_blast_page():
    open_webpage("https://blast.ncbi.nlm.nih.gov/Blast.cgi")

def ncbi_page():
    open_webpage("https://www.ncbi.nlm.nih.gov/")

def repo_page():
    open_webpage("https://github.com/TrinhLab/CASPERapp")

def open_webpage(url):
    try:
        webbrowser.open(url, new=2)
    except Exception as e:
        show_error("Error in open_webpage()", e)