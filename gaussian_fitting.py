# -*- coding: utf-8 -*-
"""
Created on Sun May 31 22:58:08 2020

@author: Dion Engels
MBx Python Data Analysis

Gaussian fitting

----------------------------

v0.1: Setup: 31/05/2020

"""
#%% Imports
import math
import numpy as np
from scipy import signal

#%% Code

def doastorm_3D():
    pass

#https://github.com/ZhuangLab/storm-analysis/tree/master/storm_analysis/daostorm_3d
    
class rainSTORM_Dion():
    
    def __init__(self, metadata, ROI_size, wavelength):
        self.result = []
        self.ROI_size = ROI_size
        self.ROI_size_1D = int((self.ROI_size-1)/2)
        self.allowSig = [0.5, self.ROI_size_1D+1]
        self.maxIts = 60
        self.maxtol = 0.2
        self.mintol = 0.005
        self.allowX = 1
        self.initSig = wavelength/(2*metadata['NA']*math.sqrt(8*math.log(2)))/(metadata['calibration_um']*1000)
        
        
    def determine_threshold(self, frame, ROI_locations, threshold):
        boolean = np.ones(frame.shape,dtype=int)
        
        for x,y in ROI_locations:
            if x < self.ROI_size_1D:
                x_min = 0
            else:
                x_min = int(x-self.ROI_size_1D)
            if x > frame.shape[0]-self.ROI_size_1D:
                x_max = frame.shape[0]
            else:
                x_max = int(x+self.ROI_size_1D)
            if y < self.ROI_size_1D:
                y_min = 0
            else:
                y_min = int(y-self.ROI_size_1D)
            if y > frame.shape[1]-self.ROI_size_1D:
                y_max = frame.shape[0]
            else:
                y_max = int(y+self.ROI_size_1D)
            boolean[x_min:x_max,y_min:y_max]  =np.zeros([x_max-x_min,y_max-y_min])
        self.ROI_zones = np.array(np.ones(boolean.shape)-boolean,dtype=bool)
        self.bg_mean = np.mean(frame[boolean==1])
        self.threshold = self.bg_mean + threshold*math.sqrt(self.bg_mean)
            
    
    def main(self, ROI_locations, frames):
        for frame_index, frame in enumerate(frames):
            frame_result = self.main_loop(frame, frame_index)
            self.result.append(frame_result)
            
        
    def main_loop(self, frame, frame_index):
        
        peaks = self.find_peaks(frame)
        peaks = peaks[peaks[:,2]> self.threshold]
        
        return self.fitter(frame_index, frame, peaks)
        
        
        
    def find_peaks(self, frame):
        kernel = np.array([[1,1,1],
                   [1,1,1],
                   [1,1,1]]) 
        
        out = signal.convolve2d(frame, kernel, boundary='fill', mode='same')/kernel.sum()
        
        row_min = self.ROI_size_1D
        row_max = frame.shape[0]-self.ROI_size_1D
        column_min = self.ROI_size_1D
        column_max = frame.shape[1]-self.ROI_size_1D
        
        maxima = np.zeros((frame.shape[0], frame.shape[1],8), dtype=bool)
        
        maxima[row_min:row_max,column_min:column_max,0] = out[row_min:row_max, column_min:column_max] > out[row_min:row_max,column_min+1:column_max+1]
        maxima[row_min:row_max,column_min:column_max,1] = out[row_min:row_max, column_min:column_max] >= out[row_min:row_max,column_min-1:column_max-1]
        maxima[row_min:row_max,column_min:column_max,2] = out[row_min:row_max, column_min:column_max] > out[row_min+1:row_max+1,column_min:column_max]
        maxima[row_min:row_max,column_min:column_max,3] = out[row_min:row_max, column_min:column_max] >= out[row_min-1:row_max-1,column_min:column_max]
        maxima[row_min:row_max,column_min:column_max,4] = out[row_min:row_max, column_min:column_max] >= out[row_min-1:row_max-1,column_min-1:column_max-1]
        maxima[row_min:row_max,column_min:column_max,5] = out[row_min:row_max, column_min:column_max] >= out[row_min-1:row_max-1,column_min+1:column_max+1]
        maxima[row_min:row_max,column_min:column_max,6] = out[row_min:row_max, column_min:column_max] > out[row_min+1:row_max+1,column_min-1:column_max-1]
        maxima[row_min:row_max,column_min:column_max,7] = out[row_min:row_max, column_min:column_max] >= out[row_min+1:row_max+1,column_min+1:column_max+1]
        
        mask = maxima.all(axis=2)
        mask = mask*self.ROI_zones
        indices = np.where(mask == True)
        indices = np.asarray([x for x in zip(indices[0],indices[1])])
        values = [[value] for value in out[mask]]
        
        return np.append(indices,values, axis=1)
        
    
    def fitter(self, frame_index, frame, peaks):
        pass
        


def MaxBergkamp():
    pass

