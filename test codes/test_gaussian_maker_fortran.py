# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 15:55:56 2020

@author: s150127
"""

import numpy as np
import time


#%% Settings 

p = [1, 4, 4, 3/(2*np.sqrt(2*np.log(2))), 3/(2*np.sqrt(2*np.log(2)))]

data = np.zeros((9,9));

num_loop = 100000

loops = list(range(0, num_loop))

#%% Not functions but actual gaussians
def makeGaussian(size, fwhm = 3, center=None):

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    x0 = center[0]
    y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

start = time.time()
    
for loop in loops:
    roi1 = makeGaussian(9, center=[4, 4])

print('Time taken v1: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v2

def makeGaussian_v2(size, fwhm = 3, center=None):

    x=np.arange(size)[None].astype(np.float)
    y=x.T
    
    x0 = center[0]
    y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

start = time.time()

for loop in loops:
    roi2 = makeGaussian_v2(9, center=[4, 4])

print('Time taken v2: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v3

# def makeGaussian_v3(size, fwhm = 3, center=None):

#     x=np.arange(size)[None].astype(np.float)
#     y=x.T
    
#     xx,yy=np.meshgrid(x,y)
    
#     x0 = center[0]
#     y0 = center[1]

#     return np.exp(-4*np.log(2) * ((xx-x0)**2 +(yy-y0)**2)/fwhm**2) #x and y vectors

# start = time.time()

# for loop in loops:
#     roi3 = makeGaussian_v3(9, center=[4, 4])

# print('Time taken v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v4
import gauss4 as gauss

start = time.time()

for loop in loops:
    
    roi4 = gauss.gaussian(*p)
    
print('Time taken v4: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v5
import gauss_full3 as gauss2

start = time.time()

for loop in loops:
    
    roi5 = gauss2.gaussian(*p, 9)
    
roi5 = np.reshape(roi5, (9,9))

print('Time taken v5: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

#%% v6
import gauss_full20 as gauss3

size = 9

x=np.arange(size)[None].astype(np.int)
x_res = np.zeros(81, dtype=int)
y_res = np.zeros(81, dtype=int)


for i in range(0,size):
    y = np.ones(size)*i
    x_res[i*size:(i+1)*size] = x
    y_res[i*size:(i+1)*size] = y
    
start = time.time()

for loop in loops:
    
    roi6 = gauss3.gaussian(*p, 9)
    
roi6 = np.reshape(roi6, (9,9))
print('Time taken v6: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))