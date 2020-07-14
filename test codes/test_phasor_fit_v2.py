# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 12:21:19 2020

@author: s150127
"""
from scipy.fftpack import ifftn
import numpy as np
from math import pi
import math
import cmath
import matplotlib.pyplot as plt
import time # for timekeeping

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
 
import pyfftw
    
def build_in_v4(roi, num_loop):
    
    roi_bb = pyfftw.empty_aligned(roi.shape, dtype='float64')
    for i in range(0,size*10):
        for j in range(0, size*10):
            try:
                roi_bf = pyfftw.empty_aligned((num_loop, i, j), dtype='complex128')
                fft_values_list = pyfftw.FFTW(roi_bb, roi_bf,axes=(1,2),flags=('FFTW_MEASURE',), direction='FFTW_FORWARD')
                print(str(i) + str(j))
                roi_bb=roi
                fft_values_list = fft_values_list(roi_bb)
            except:
                continue
            
    return 0, 0, 0
       
    # for fft_values in fft_values_list:
        
    
    #     roi_size = roi.shape[1];
    #     ang_x = cmath.phase(fft_values[0, 1])
    #     if ang_x>0:
    #         ang_x=ang_x-2*pi
    
    #     pos_x = abs(ang_x)/(2*pi/roi_size)
    
    #     ang_y = cmath.phase(fft_values[1,0])
    
    #     if ang_y >0:
    #         ang_y = ang_y - 2*pi
    
    #     pos_y = abs(ang_y)/(2*pi/roi_size)
    
    #     if pos_x > 8.5:
    #         pos_x -= roi_size
    #     if pos_y > 8.5:
    #         pos_y -= roi_size

    # return [pos_x, pos_y, fft_values]


roi = makeGaussian(size)

num_loop = 1000

loops = list(range(0, num_loop))

start = time.time()

roi_total = np.ones((num_loop,size,size))

for loop in loops:
    roi_total[loop,:,:] = roi
    
start = time.time()
    
pos_x5, pos_y5, fft_values5 = build_in_v4(roi_total, num_loop)

print('Time taken build-in v4: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))


plt.imshow(roi, extent=[0,roi.shape[0],roi.shape[0],0], aspect='auto')
plt.scatter(pos_x5+.5, pos_y5+.5, s=500, c='red', marker='x', alpha=0.5)
plt.title("Test")
plt.show()