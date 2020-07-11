# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:20:29 2020

@author: Dion Engels
MBx Python Data Analysis

Frame Analysis

----------------------------

v1.0, ROI detection:  31/05/2020
v1.1, conventional naming: 04/06/2020

"""
import time # to sleep to ensure that matlab licence is registered
import numpy as np # for linear algebra
import matlab.engine # to run matlab
import scipy.io as sio # to load in .mat


 #%% Adapted from FindParticles.m from SPectraA4

def ROI_finder(image, ROI_size):
    """
    Finds ROIs

    Parameters
    ----------
    Image : Image to find ROIs for
    ROI_size : Size of ROI

    Returns
    -------
    Population : ROI info

    """

    threshold = [0.05, 0.75]

    eng = matlab.engine.start_matlab()

    image_list = image.tolist()
    image_matlab = matlab.double(image_list)

    time.sleep(3)

    [corrected, background] = background_correction_wavelet_set(image_matlab, 4, eng)
    image = image - background

    image_list = image.tolist()
    image_matlab = matlab.double(image_list)

    bead_filter_dict = sio.loadmat('Typical_BeadFilter_for_Finding_AuNPs.mat')
    bead_filter = bead_filter_dict['BeadFilter']
    bead_filter_list = bead_filter.tolist()
    bead_filter_matlab = matlab.double(bead_filter_list)

    population_dict_matlab = eng.LocateBeadsProfile(image_matlab, bead_filter_matlab, False, threshold[0], threshold[1], 20)

    population_matlab = population_dict_matlab['Location']
    population = np.array(population_matlab)

    eng.quit()

    return population




 #%% Runs BackgroundCorrection_WaveletSet Version 1.2 by Emiel Visser

def background_correction_wavelet_set(frame_data, level, eng):
    """
    Runs BackgroundCorrection_WaveletSet Version 1.2 by Emiel Visser

    Parameters
    ----------
    FrameData : Info about frame
    Level : Some parameters
    eng : Matlab engine

    Returns
    -------
    list : hipass en lowpass info

    """

    hi_pass_matlab =  eng.BackgroundCorrection_WaveletSet(frame_data, level)

    hi_pass = np.array(hi_pass_matlab)
    low_pass = frame_data - hi_pass

    return [hi_pass, low_pass]

#%% Python ROI finder

from scipy.stats import norm
import scipy.ndimage.filters as filters

class roi_finder():
    
    def __init__(self, intensity_min, intensity_max, sigma_min, sigma_max,
                 shape, symmetry, roi_size):
        self.intensity_min = intensity_min
        self.intensity_max = intensity_max
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        self.symmetry = symmetry
        self.shape = shape
        
        self.roi_size = int(roi_size)
        self.roi_size_1d = int((roi_size-1)/2)
        self.roi_locations = []
        
    def determine_threshold_min(self, frame):
        
        frame_ravel = np.ravel(frame)
        for _ in range(4):
            mean, std = norm.fit(frame_ravel)
            frame_ravel = frame_ravel[frame_ravel < mean+std*5]
        
        return mean+self.threshold_sigma*std
        
    def determine_standard_values(self, frame):
        
        self.threshold_sigma = 5
        self.intensity_min = self.determine_threshold_min(frame)
        self.intensity_max = np.inf
        self.sigma_min = 0
        self.sigma_max = np.inf
        self.symmetry = 0.8
        self.shape = 0.3
        
    def find_within_intensity_range(self, frame):
        
        boolean_int_min = frame > self.intensity_min
        boolean_int_max = frame < self.intensity_max
        
        return boolean_int_max & boolean_int_min
        
    def find_local_maximum(self, frame):
        
        data_max = filters.maximum_filter(frame, self.roi_size)
        maxima = (frame == data_max)
        
        return maxima
        
    def adjacent_or_boundary_rois(self, roi_boolean):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            
            y = int(roi[0])
            x = int(roi[1])

            try:            
                my_roi = roi_boolean[y-self.roi_size_1d:y+self.roi_size_1d+1, x-self.roi_size_1d:x+self.roi_size_1d+1]
            except:
                keep_boolean[roi_index] = False # if this fails, the roi is on the boundary
                continue
            
            trues_in_roi = np.transpose(np.where(my_roi == True))
            
            if trues_in_roi.shape[0] > 1:
                keep_boolean[roi_index] = False
                
        self.roi_locations = self.roi_locations[keep_boolean, :]
        
    def remove_boundary_rois(self, frame):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        x_low = self.roi_size
        x_high = frame.shape[1] - self.roi_size
        
        y_low = self.roi_size
        y_high = frame.shape[0] - self.roi_size
        
        for roi_index, roi in enumerate(self.roi_locations):
            
            y = int(roi[0])
            x = int(roi[1]) 
            if x < x_low or x > x_high:
                keep_boolean[roi_index] = False
            if y < y_low or y > y_high:
                keep_boolean[roi_index] = False
        
        self.roi_locations = self.roi_locations[keep_boolean, :] 
    
    def main(self, frame):
        
        roi_boolean = self.find_within_intensity_range(frame)
        
        local_maxima = self.find_local_maximum(frame)
        
        roi_boolean = roi_boolean & local_maxima
        
        self.roi_locations = np.transpose(np.where(roi_boolean == True))
        
        self.adjacent_or_boundary_rois(roi_boolean) # check if any rois too close
        
        self.remove_boundary_rois(frame) # remove rois too close to edge
        
        #self.roi_symmetry(frame) # check symmetry of roi
        
        #self.roi_shape(frame) # check shape of roi
        
        
        
        
        
        
        
        return self.roi_locations
        