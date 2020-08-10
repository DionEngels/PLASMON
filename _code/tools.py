# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
MBx Python Data Analysis

Tools

Some additional tools used by MBx Python

----------------------------

v0.1: Save to CSV & Mat: 31/05/2020
v0.2: also switch array: 04/06/2020
v0.3: cleaned up: 24/07/2020
v0.4: settings and results text output: 25/07/2020
v0.5: own ND2Reader class to prevent warnings: 29/07/2020
v0.6: save drift and save figures: 03/08/2020
v0.6.1: better MATLAB ROI and Drift output
v0.7: metadata v2
v0.8: metadata v3
v1.0: bug fixes
v1.0.1: metadata goldmine found, to be implemented
v1.1: emptied out: 09/08/2020

"""
from numpy import zeros


def switch_to_matlab_coordinates(array, height):

    array = switch_axis(array)
    array[:, 1] = height - array[:, 1]

    return array


def switch(array):
    """
    Switches a single arrays values

    Parameters
    ----------
    array : array to switch

    Returns
    -------
    new : switched array

    """
    new = zeros(array.shape)
    new[:, 1] = array[:, 0]
    new[:, 0] = array[:, 1]
    return new
