# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 18:42 2020

@author: Dion Engels
MBx Python Data Analysis

Setup code

This piece of code allows you to compile Mbx Python to an .exe

----------------------------

v0.1: Full setup: 26/07/2020
v0.2: minor changes
v0.2.1: trying to get lower file size, no progress
v0.3: adding in dependencies to ensure working
v0.4: also copy spectral corrections
v1.0: initial release done
v1.2: icon

"""

import os
import sys
from cx_Freeze import setup, Executable


r'-b C:\Users\s150127\Downloads\___MBx\build'  # take this control location of build

__version__ = '1.4'

include_files = ['spectral_corrections/', 'ico.ico']
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

base = None
if sys.platform == 'win32':
    base = 'Win32GUI'
    DLLS_FOLDER = os.path.join(PYTHON_INSTALL_DIR, 'Library', 'bin')

    dependencies = ['libiomp5md.dll', 'mkl_core.dll', 'mkl_def.dll', 'mkl_intel_thread.dll']

    for dependency in dependencies:
        include_files.append(os.path.join(DLLS_FOLDER, dependency))


includes = ['tkinter']
excludes = ['matplotlib.tests', 'numpy.random._examples']
packages = ['numpy', 'matplotlib', 'multiprocessing', 'scipy']

setup(
    name='MbxPython',
    description='MbxPython',
    version=__version__,
    executables=[Executable('main_gui.py', targetName="MbxPython.exe", base=base, icon='ico.ico')],
    options={'build_exe': {
        'packages': packages,
        'includes': includes,
        'include_files': include_files,
        'include_msvcr': True,
        'excludes': excludes,
    }},
)

path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
build_path = os.path.join(path, 'build')
