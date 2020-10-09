# -*- coding: utf-8 -*-
"""
Created on Fri Oct 09

----------------------------

@author: Dion Engels
MBx Python Data Analysis

warnings
----------------------

Holds some warnings I can raise during runtime
"""


class InputWarning(UserWarning):  # class for input warnings, when the user inputs something stupid
    pass


class DataWarning(UserWarning):  # class for dubious resources, such as missing metadata
    pass
