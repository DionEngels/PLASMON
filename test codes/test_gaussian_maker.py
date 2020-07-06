# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 15:55:56 2020

@author: s150127
"""

import numpy as np
import time
from functools import lru_cache


#%% Settings 

p = [1000, 4, 4, 1, 1]

data = np.zeros((9,9));

num_loop = 100000

loops = list(range(0, num_loop))

#%% v1

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


start = time.time()

indices = np.indices(data.shape)
errorfunction1 = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)))
    
for loop in loops:
    roi1 = errorfunction1(p)

roi1 = np.reshape(roi1,(9,9))

print('Time taken v1: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v2

start = time.time()

indices = np.indices(data.shape)
errorfunction2 = lambda p: np.ravel(gaussian(*p)(*indices))  
    
for loop in loops:
    roi2 = errorfunction2(p)
    
roi2 = np.reshape(roi2,(9,9))

print('Time taken v2: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v3

@lru_cache(maxsize=1000)
def gaussian_v2(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

start = time.time()

indices = np.indices(data.shape)
errorfunction3 = lambda p: np.reshape(gaussian_v2(*p)(*indices),-1)  
    
for loop in loops:
    roi3 = errorfunction3(p)
    
roi3 = np.reshape(roi2,(9,9))

print('Time taken v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))


print(errorfunction1 == errorfunction2)

#%% Not functions but actual gaussians
# def makeGaussian(size, fwhm = 3, center=None):

#     x = np.arange(0, size, 1, float)
#     y = x[:,np.newaxis]
    
#     x0 = center[0]
#     y0 = center[1]

#     return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

# start = time.time()
    
# for loop in loops:
#     roi1 = makeGaussian(9, center=[4, 4])

# print('Time taken v1: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

# #%% v2

# def makeGaussian_v2(size, fwhm = 3, center=None):

#     x=np.arange(size)[None].astype(np.float)
#     y=x.T
    
#     x0 = center[0]
#     y0 = center[1]

#     return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

# for loop in loops:
#     roi2 = makeGaussian_v2(9, center=[4, 4])

# print('Time taken v2: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

# #%% v3

# def makeGaussian_v3(size, fwhm = 3, center=None):

#     x=np.arange(size)[None].astype(np.float)
#     y=x.T
    
#     xx,yy=np.meshgrid(x,y)
    
#     x0 = center[0]
#     y0 = center[1]

#     return np.exp(-4*np.log(2) * ((xx-x0)**2 +(yy-y0)**2)/fwhm**2) #x and y vectors

# for loop in loops:
#     roi3 = makeGaussian_v3(9, center=[4, 4])

# print('Time taken v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))