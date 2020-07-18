# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:09:00 2020

@author: s150127
"""

import numpy as np
import time

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

roi = makeGaussian(9, center=[4, 4])


num_loop = 1000000

loops = list(range(0, num_loop))

start = time.time()

for loop in loops:
    minimum = np.min(roi)

print('Time taken NP ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

import MBx_FORTRAN_TOOLS_v2 as fortran

start = time.time()
    
for loop in loops:
    minimum2 = fortran.min9(roi)

print('Time taken FORTRAN: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))


