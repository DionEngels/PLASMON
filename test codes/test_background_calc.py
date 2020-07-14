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
    roi_bg1 = np.mean(np.append(np.append(np.append(roi[:, 0], 
    roi[:, -1]), np.transpose(roi[0, 1:-1])), np.transpose(roi[-1, 1:-1])))

print('Time taken build-in: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

def determine_background(empty_background, my_roi):
        roi_background = empty_background
        roi_background[0:9] = my_roi[:, 0]
        roi_background[9:9*2] = my_roi[:, -1]
        roi_background[9*2:9*2+9-2] = my_roi[0, 1:-1]
        roi_background[9*2+9-2:] = my_roi[-1, 1:-1]
        
        return np.mean(roi_background)

start = time.time()

empty_background = np.zeros(9*2+(9-2)*2)
    
for loop in loops:
    roi_bg2 = determine_background(empty_background, roi)

print('Time taken build-in v2: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

import MBx_FORTRAN_TEST_v1 as fortran

start = time.time()
    
for loop in loops:
    roi_bg3 = fortran.calc_max(roi)

print('Time taken FORTRAN: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))


