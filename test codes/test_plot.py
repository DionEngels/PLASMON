# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 13:13:07 2020

@author: s150127
"""

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    p2 = optimize.least_squares(errorfunction, params)

    return [p2.x, p2.nfev]







# Create the gaussian data
Xin, Yin = np.mgrid[0:9, 0:9]
data = gaussian(3, 2, 4, 2, 2)(Xin, Yin) #+ np.random.random(Xin.shape)

loops = range(0,1)

for loop in loops:
    params, nfev = fitgaussian(data)


fit = gaussian(*params)

#plt.matshow(data)
plt.imshow(data, extent=[0,data.shape[0],data.shape[0],0], aspect='auto')
plt.scatter(params[2]+.5, params[1]+.5, s=500, c='red', marker='x', alpha=0.5)
plt.title("Test")
plt.show()