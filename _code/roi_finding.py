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
    
    def determine_threshold_min(self, frame):
        
        frame_ravel = np.ravel(frame)
        for _ in range(4):
            mean, std = norm.fit(frame_ravel)
            frame_ravel = frame_ravel[frame_ravel < mean+std*5]
        
        return mean+self.threshold_sigma*std
    
    def __init__(self, roi_size, frame, intensity_min = None, intensity_max = None, 
                 sigma_min = 0, sigma_max = 10):
        
        self.threshold_sigma = 5
        
        if intensity_min == None:
            self.intensity_min = self.determine_threshold_min(frame)
        else:
            self.intensity_min = intensity_min
            
        if intensity_max == None:
            self.intensity_max = 65535 # max of unsigned integer 16 bits
        else:
            self.intensity_max = intensity_max
            
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        
        self.roi_size = int(roi_size)
        self.roi_size_1d = int((roi_size-1)/2)
        self.roi_locations = []
        
        self.empty_background = np.zeros(self.roi_size*2+(self.roi_size-2)*2, dtype=np.uint16)
             
    def change_settings(self, intensity_min = None, 
                        intensity_max = None,
                        sigma_min = 0, sigma_max = 10):
        if intensity_min != None:
            self.intensity_min = intensity_min
        if intensity_max != None:
            self.intensity_max = intensity_max
        
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        
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
                my_roi = roi_boolean[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                                     x-self.roi_size_1d:x+self.roi_size_1d+1]
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
        
    def determine_background(self, my_roi):
       
        roi_background = self.empty_background
        roi_background[0:self.roi_size] = my_roi[:, 0]
        roi_background[self.roi_size:self.roi_size*2] = my_roi[:, -1]
        roi_background[self.roi_size*2:self.roi_size*2+self.roi_size-2] = my_roi[0, 0:-2]
        roi_background[self.roi_size*2+self.roi_size-2:] = my_roi[-1, 0:-2]
        
        return np.mean(roi_background)
        
    def roi_symmetry(self, frame):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            x_range = np.asarray(range(0, self.roi_size))+0.5
            y_range = np.asarray(range(0, self.roi_size))+0.5
            
            y = int(roi[0])
            x = int(roi[1]) 
            
            my_roi = frame[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                           x-self.roi_size_1d:x+self.roi_size_1d+1]
            
            #my_roi_bg = self.determine_background(my_roi)
            #my_roi = my_roi - my_roi_bg
            
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
            
            if roi_symmetry < self.symmetry:
                keep_boolean[roi_index] = False
            
        self.roi_locations = self.roi_locations[keep_boolean, :] 
            
    def roi_shape(self, frame):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        shape = (self.roi_locations.shape[0], self.roi_size, self.roi_size)
        stack = np.zeros(shape)
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            my_roi = frame[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                           x-self.roi_size_1d:x+self.roi_size_1d+1]
            
            my_roi_bg = self.determine_background(my_roi)
            stack[roi_index, :, :] = my_roi - my_roi_bg
        
        mean = np.mean(stack, axis=0)
        mean_sum = np.sum(mean)
        mean = mean/np.max(mean) #normalize
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            my_roi = frame[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                           x-self.roi_size_1d:x+self.roi_size_1d+1]
            
            my_roi_bg = self.determine_background(my_roi)
            my_roi = my_roi - my_roi_bg
            
            remain = my_roi - mean*np.max(my_roi)
            
            remain_sum = np.sum(remain)
            shape_factor = remain_sum / mean_sum
            
            if shape_factor > self.shape:
                keep_boolean[roi_index] = False
            
        
        self.roi_locations = self.roi_locations[keep_boolean, :] 
        
    def sigma_limit(self, frame, fitter):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            my_roi = frame[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                           x-self.roi_size_1d:x+self.roi_size_1d+1]
            my_roi_bg = self.determine_background(my_roi)
            my_roi = my_roi - my_roi_bg
        
            result, its, success = fitter.fitgaussian(my_roi, roi_index)
            
            if result[4] < self.sigma_min or result[3] < self.sigma_min:
                keep_boolean[roi_index] = False
            if result[4] > self.sigma_max or result[3] > self.sigma_max:
                keep_boolean[roi_index] = False
                
        self.roi_locations = self.roi_locations[keep_boolean, :] 
    
    def main(self, frame, fitter):
        
        roi_boolean = self.find_within_intensity_range(frame)
        
        local_maxima = self.find_local_maximum(frame)
        
        roi_boolean = roi_boolean & local_maxima
        
        self.roi_locations = np.transpose(np.where(roi_boolean == True))
        
        self.adjacent_or_boundary_rois(roi_boolean) # check if any rois too close
        
        self.remove_boundary_rois(frame) # remove rois too close to edge
        
        #self.roi_symmetry(frame) # check symmetry of roi
        
        #self.roi_shape(frame) # check shape of roi
        
        self.sigma_limit(frame, fitter)        
        
        return self.roi_locations
        