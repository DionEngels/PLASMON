# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 14:32:55 2020

@author: s150127
"""

import numpy as np
import MBx_FORTRAN_TEST_v6 as fortran

size = 7

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

roi = makeGaussian(size, center=[4, 4])

if size == 9:
    maximum = lambda roi: fortran.max9(roi)
else:
    maximum = lambda roi: fortran.max7(roi)
    
print(maximum(roi))