# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 18:42 2020

@author: Dion Engels
MBx Python Data Analysis

Setup code

----------------------------

v1: Full setup: 26/07/2020

"""
_ = r'-b C:\Users\s150127\Downloads\__build_test'

import os
import shutil
import sys
from cx_Freeze import setup, Executable

#  os.environ['TCL_LIBRARY'] = r'C:\bin\Python37-32\tcl\tcl8.6'
#  os.environ['TK_LIBRARY'] = r'C:\bin\Python37-32\tcl\tk8.6'

__version__ = '1.0.0'
base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

include_files = []
includes = ['tkinter']
excludes = ['matplotlib.tests', 'numpy.random._examples']
packages = ['numpy', 'scipy', 'matplotlib']

setup(
    name='MbxPython',
    description='Fast localization from .nd2 data',
    version=__version__,
    executables=[Executable('main_gui.py', base=base)],
    options={'build_exe': {
        'packages': packages,
        'includes': includes,
        'include_files': include_files,
        'include_msvcr': True,
        'excludes': excludes,
    }},
)

path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
build_path = os.path.join(path, 'build', 'exe.win32-3.7')
#  shutil.copy(r'C:\bin\Python37-32\DLLs\tcl86t.dll', build_path)
#  shutil.copy(r'C:\bin\Python37-32\DLLs\tk86t.dll', build_path)
