# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 21:44:46 2020

@author: s150127
"""

import numpy as np
from scipy import signal

ROI_size= 4

# Inputs
arr = np.random.rand(10,10)

# Define kernel for convolution                                         
kernel = np.array([[1,1,1],
                   [1,1,1],
                   [1,1,1]]) 

# Perform 2D convolution with input data and kernel 
out = signal.convolve2d(arr, kernel, boundary='fill', mode='same')/kernel.sum()

row_min= ROI_size
row_max = arr.shape[0]-ROI_size
column_min = ROI_size
column_max = arr.shape[1]-ROI_size

# rows = range(ROI_size,arr.shape[0]-ROI_size)
# columns = range(ROI_size,arr.shape[1]-ROI_size)

maxima = np.zeros((arr.shape[0], arr.shape[1],8), dtype=bool)

#test3 = signal.argrelextrema(out,np.greater,axis=1)

maxima[row_min:row_max,column_min:column_max,0] = out[row_min:row_max, column_min:column_max] > out[row_min:row_max,column_min+1:column_max+1]
maxima[row_min:row_max,column_min:column_max,1] = out[row_min:row_max, column_min:column_max] >= out[row_min:row_max,column_min-1:column_max-1]
maxima[row_min:row_max,column_min:column_max,2] = out[row_min:row_max, column_min:column_max] > out[row_min+1:row_max+1,column_min:column_max]
maxima[row_min:row_max,column_min:column_max,3] = out[row_min:row_max, column_min:column_max] >= out[row_min-1:row_max-1,column_min:column_max]
maxima[row_min:row_max,column_min:column_max,4] = out[row_min:row_max, column_min:column_max] >= out[row_min-1:row_max-1,column_min-1:column_max-1]
maxima[row_min:row_max,column_min:column_max,5] = out[row_min:row_max, column_min:column_max] >= out[row_min-1:row_max-1,column_min+1:column_max+1]
maxima[row_min:row_max,column_min:column_max,6] = out[row_min:row_max, column_min:column_max] > out[row_min+1:row_max+1,column_min-1:column_max-1]
maxima[row_min:row_max,column_min:column_max,7] = out[row_min:row_max, column_min:column_max] >= out[row_min+1:row_max+1,column_min+1:column_max+1]

mask = maxima.all(axis=2)
#%%
indices = np.where(mask == True)
indices = np.asarray([x for x in zip(indices[0],indices[1])])
values = [[value] for value in out[mask]]

result = np.append(indices,values, axis=1)