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


def build_in(roi):
    fft_values = ifftn(roi)

    roi_size = roi.shape[0];
    ang_x = cmath.phase(fft_values[0, 1])
    if ang_x>0:
        ang_x=ang_x-2*pi

    pos_x = roi_size - abs(ang_x)/(2*pi/roi_size)

    ang_y = cmath.phase(fft_values[1,0])

    if ang_y >0:
        ang_y = ang_y - 2*pi

    pos_y = roi_size - abs(ang_y)/(2*pi/roi_size)

    if pos_x > 8.5:
        pos_x -= roi_size
    if pos_y > 8.5:
        pos_y -= roi_size

    return [pos_x, pos_y]


def self_made(roi):
    roi_int = int(roi.shape[0])
    roi_double = float(roi.shape[0])

    fitomega = (np.array(range(0,roi_int))+1)*2*pi/roi_double
    fitcos = np.cos(fitomega)
    fitsin = np.sin(fitomega)

    tot_x_loop_re = np.zeros(fitomega.shape)
    tot_x_loop_im = np.zeros(fitomega.shape)
    tot_y_loop_re = np.zeros(fitomega.shape)
    tot_y_loop_im = np.zeros(fitomega.shape)

    first_fourier_coeff_x_re = 0
    first_fourier_coeff_x_im = 0
    first_fourier_coeff_y_re = 0
    first_fourier_coeff_y_im = 0

    for i, value in enumerate(tot_x_loop_re):
        for j, value in enumerate(tot_y_loop_re):
            tot_x_loop_re[j] +=  fitcos[j]*roi[i, j];
            tot_x_loop_im[j] -= fitsin[j]*roi[i, j];
            tot_y_loop_re[i] += fitcos[i]*roi[i, j];
            tot_y_loop_im[i] -= fitsin[i]*roi[i, j];

    first_fourier_coeff_x_re = np.sum(tot_x_loop_re)
    first_fourier_coeff_x_im = np.sum(tot_x_loop_im)
    first_fourier_coeff_y_re = np.sum(tot_y_loop_re)
    first_fourier_coeff_y_im = np.sum(tot_y_loop_im)

    ang_x = math.atan2(first_fourier_coeff_x_im,first_fourier_coeff_x_re)
    if ang_x>0:
        ang_x=ang_x-2*pi

    pos_x = abs(ang_x)/(2*pi/roi_double) - 1

    ang_y = math.atan2(first_fourier_coeff_y_im,first_fourier_coeff_y_re)

    if ang_y >0:
        ang_y = ang_y - 2*pi

    pos_y = abs(ang_y)/(2*pi/roi_double) - 1

    return [pos_x, pos_y]



roi = makeGaussian(9, center=[3, 3])

loops = list(range(0, 1))

for loop in loops:
    [pos_x1, pos_y1] = build_in(roi)


for loop in loops:
    [pos_x2, pos_y2] = self_made(roi)



plt.imshow(roi, extent=[0,roi.shape[0],roi.shape[0],0], aspect='auto')
plt.scatter(pos_x1+.5, pos_y1+.5, s=500, c='red', marker='x', alpha=0.5)
plt.title("Test")
plt.show()