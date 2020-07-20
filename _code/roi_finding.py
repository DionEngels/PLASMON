# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:20:29 2020

@author: Dion Engels
MBx Python Data Analysis

Frame Analysis

----------------------------

v1.0, ROI detection:  31/05/2020
v1.1, conventional naming: 04/06/2020
v2.0: self-made ROI finding: 10/07/2020
v2.1: tweaks with standard min intensity: 14/07/2020
v2.2: tweaked max sigma and max intensity: 17/07/2020
v2.3: return sigmas option for histograms
v2.4: tweaking removing ROIs with nearby ROIs: 18/07/2020
v2.5: worked on main_v2, which uses correlation of 2D gaussian
v3.0: clean up and based on SPectrA; correlation, pixel_int, sigma and int

"""
import numpy as np # for linear algebra
from scipy.signal import medfilt, convolve2d
from scipy.stats import norm
from scipy.ndimage.filters import maximum_filter

#%% Python ROI finder

class roi_finder():
    
    def determine_threshold_min(self, frame):
        
        background = medfilt(frame, kernel_size = self.roi_size)
        background[background ==  0] = np.min(background[background > 0])       
        frame_filter = frame.astype('float') - background
        
        frame_ravel = np.ravel(frame)
        frame_filter_ravel = np.ravel(frame_filter)
        mean, std = norm.fit(frame_filter_ravel)
        
        for _ in range(10):
            mean, std = norm.fit(frame_ravel)
            frame_ravel = frame_ravel[frame_ravel < mean+std*5]
        
        return mean+self.threshold_sigma*std
    
    def __init__(self, roi_size, frame, pixel_min = None, corr_min = 0.05, 
                 sigma_min = 0, sigma_max = 5, int_min = 0, int_max = np.inf):
        
        self.roi_size = int(roi_size)
        self.roi_size_1d = int((roi_size-1)/2)
        self.side_distance = self.roi_size + 2
        self.roi_distance = self.roi_size_1d + 2
        
        
        self.threshold_sigma = 3
        if pixel_min == None:
            self.pixel_min = self.determine_threshold_min(frame)
        else:
            self.pixel_min = pixel_min
            
        self.int_min = int_min
        self.int_max = int_max
            
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        
        self.corr_min = corr_min
            
        self.roi_locations = []
        self.sigma_list = []
        self.int_list = []
             
    def change_settings(self, intensity_min = None, 
                        intensity_max = None,
                        sigma_min = 0, sigma_max = 10):
        if intensity_min != None:
            self.intensity_min = intensity_min
        if intensity_max != None:
            self.intensity_max = intensity_max
        
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        
    def adjacent_or_boundary_rois_base(self, roi_boolean):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            
            y = int(roi[0])
            x = int(roi[1])

            my_roi = roi_boolean[y-self.side_distance:y+self.side_distance+1, 
                                     x-self.side_distance:x+self.side_distance+1]
            if my_roi.shape != (self.side_distance*2+1, self.side_distance*2+1):
                keep_boolean[roi_index] = False # if this fails, the roi is on the boundary
                continue
                                
            my_roi = roi_boolean[y-self.roi_distance:y+self.roi_distance+1, 
                                     x-self.roi_distance:x+self.roi_distance+1]
            
            trues_in_roi = np.transpose(np.where(my_roi == True))
            
            if trues_in_roi.shape[0] > 1:
                keep_boolean[roi_index] = False
                
        self.roi_locations = self.roi_locations[keep_boolean, :]
        
    def int_sigma_limit(self, frame, fitter, return_int, return_sigmas):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            my_roi = frame[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                           x-self.roi_size_1d:x+self.roi_size_1d+1]
        
            result, its, success = fitter.fitgaussian(my_roi, roi_index)
            
            if return_sigmas:
                self.sigma_list.append(result[4])
                self.sigma_list.append(result[3])
            if return_int:
                self.int_list.append(result[0])
            
            if result[4] < self.sigma_min or result[3] < self.sigma_min:
                keep_boolean[roi_index] = False
            if result[4] > self.sigma_max or result[3] > self.sigma_max:
                keep_boolean[roi_index] = False
            if result[0] < self.int_min or result[0] > self.int_max:
                keep_boolean[roi_index] = False
                
        self.roi_locations = self.roi_locations[keep_boolean, :] 
        
    def make_gaussian(self, size, fwhm = 3, center=None):
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
        
    def find_particles(self, frame):
        
        background = medfilt(frame, kernel_size = self.roi_size)
        background[background ==  0] = np.min(background[background > 0])
        
        frame = frame.astype('float') - background
        
        compare = self.make_gaussian(self.roi_size)
        
        frame_convolution = convolve2d(frame, compare, mode='same')
        
        fp = np.ones((3,3), dtype=bool)
        local_peaks = maximum_filter(frame_convolution, footprint=fp)
        local_peaks_bool = (frame_convolution == local_peaks)
        
        max_convolution = np.max(frame_convolution)
        
        beads = (frame_convolution * local_peaks_bool) > self.corr_min*max_convolution
                
        return frame, beads
    
    def peak_int_min(self, frame):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            peak_int = frame[y, x]
            
            if peak_int < self.pixel_min:
                keep_boolean[roi_index] = False
                
        self.roi_locations = self.roi_locations[keep_boolean, :] 
    
    def main(self, frame, fitter, return_int = False, return_sigmas = False):
        
        frame, roi_boolean = self.find_particles(frame)
        
        self.roi_locations = np.transpose(np.where(roi_boolean == True))
        
        self.adjacent_or_boundary_rois_base(roi_boolean)
        
        self.peak_int_min(frame)
                       
        self.int_sigma_limit(frame, fitter, return_int, return_sigmas)
        
        if return_sigmas:
            return self.sigma_list
        elif return_int:
            return self.int_list
        else:
            return self.roi_locations 
        
        
        
        