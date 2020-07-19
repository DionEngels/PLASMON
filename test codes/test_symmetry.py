# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:09:00 2020

@author: s150127
"""

import numpy as np

def determine_background(size, my_roi):
       
        roi_background = np.zeros(size*2+(size-2)*2, dtype=np.uint16)
        roi_background[0:size] = my_roi[:, 0]
        roi_background[size:size*2] = my_roi[:, -1]
        roi_background[size*2:size*2+size-2] = my_roi[0, 0:-2]
        roi_background[size*2+size-2:] = my_roi[-1, 0:-2]
        
        return np.mean(roi_background)

size = 9

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

my_roi = 600*makeGaussian(9, center=[4, 4])

my_roi +=200

my_roi += 200*makeGaussian(9, center=[1, 4])

x_range = np.asarray(range(0, size))+0.5
y_range = np.asarray(range(0, size))+0.5

my_roi_bg = determine_background(size, my_roi)
my_roi = my_roi - my_roi_bg

x_sum = np.sum(my_roi, axis=0)
y_sum = np.sum(my_roi, axis=1)
total = np.sum(x_sum)

mu_x = np.sum(x_range*x_sum)/total
mu_y = np.sum(y_range*y_sum)/total

x_range = x_range - mu_x
y_range = y_range - mu_y
x_mesh, y_mesh = np.meshgrid(x_range, y_range)

roi_symmetry_xy = np.sum(my_roi*x_mesh*y_mesh)/(total-1)
roi_symmetry_xx = np.sum(my_roi*x_mesh*x_mesh)/(total-1)
roi_symmetry_yx = np.sum(my_roi*y_mesh*x_mesh)/(total-1)
roi_symmetry_yy = np.sum(my_roi*y_mesh*y_mesh)/(total-1)

C = np.array([[roi_symmetry_xx, roi_symmetry_xy], 
     [roi_symmetry_yx, roi_symmetry_yy]])

eigenvalues = np.linalg.eigvals(C)

roi_symmetry = np.sqrt(np.min(eigenvalues))/np.sqrt(np.max(eigenvalues))

