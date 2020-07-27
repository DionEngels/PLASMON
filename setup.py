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

"""

import os
import sys
from cx_Freeze import setup, Executable

r'-b C:\Users\s150127\Downloads\__build_test'  # take this control location of build

__version__ = '1.0'
base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

include_files = []
includes = ['tkinter']
excludes = ['matplotlib.tests', 'numpy.random._examples']
packages = ['numpy', 'matplotlib', 'multiprocessing', 'scipy']

setup(
    name='MbxPython',
    description='MbxPython',
    version=__version__,
    executables=[Executable('main_gui.py', targetName="MbxPython.exe", base=base)],
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