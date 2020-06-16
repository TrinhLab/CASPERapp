"""
This is a setup.py script for cx_freeze running on Windows

Usage:
    python setupwindows.py build
    python setupwindows.py bdist_msi
"""

import sys
from cx_freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {'include_files':
                    ['ylicodontable.txt', 'cthcodontable.txt', 'ecocodontable.txt',
                     'ctxcodontable.txt', 'codon_opt_gui.ui', 'myicon.ico']}

# GUI applications require a different base on Windows (the default is for a console application)
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(
    name = "CodonOpt",
    version = "1.0",
    description = "simple codon optimizer",
    author = "Brian Mendoza",
    options = {"build_exe": build_exe_options},
    executables = [Executable("CodonOpt.py", base=base, icon="myicon.ico")]
)
