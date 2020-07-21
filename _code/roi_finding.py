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
v3.1: adaptations to GUI
v3.2: further adaptations for GUI

"""
import numpy as np # for linear algebra
from scipy.signal import medfilt, convolve2d
from scipy.stats import norm
from scipy.ndimage.filters import maximum_filter

#%% Python ROI finder

class roi_finder():
    
    def determine_threshold_min(self):
        
        frame_ravel = np.ravel(self.frame)
        
        for _ in range(10):
            mean, std = norm.fit(frame_ravel)
            frame_ravel = frame_ravel[frame_ravel < mean+std*5]
        
        return mean+self.threshold_sigma*std
    
    def __init__(self, filter_size, frame, fitter, pixel_min = None, corr_min = 0.05, 
                 sigma_min = 0, sigma_max = None, int_min = None, int_max = None,
                 roi_size = 7):
        
        self.filter_size = int(filter_size)
        self.roi_size = roi_size
        self.roi_size_1d = int((self.roi_size-1)/2)
        self.side_distance = 11
        self.roi_distance = 6
        
        self.base_frame = frame
        
        background = medfilt(frame, kernel_size = self.filter_size)
        background[background ==  0] = np.min(background[background > 0])       
        self.frame = frame.astype('float') - background
        
        self.threshold_sigma = 5
        if pixel_min == None:
            self.pixel_min = self.determine_threshold_min()
        else:
            self.pixel_min = pixel_min
            
        self.sigma_list = []
        self.int_list = []
        self.corr_list = []
        
        self.corr_min = corr_min
            
        self.roi_locations = []
    
        self.int_min = int_min
        self.sigma_min = sigma_min
        
        if sigma_max == None or int_max == None or int_min == None:
            self.sigma_max = 5
            self.int_max = np.inf
            self.int_min = 0
            self.roi_locations = self.main(fitter)
            if int_max == None or int_min == None:
                self.int_sigma_limit(fitter, True, False)
                if int_max == None:
                    self.int_max = np.max(self.int_list)*1.05 # 5% margin
                if int_min == None:
                    self.int_min = np.min(self.int_list)/1.05 # 5% margin
            if sigma_max == None:
                self.int_sigma_limit(fitter, False, True)
                self.sigma_max = np.max(self.sigma_list)*1.05 # 5% margin                
        else:
            self.sigma_max = sigma_max
            self.int_max = int_max
        
    def change_settings(self, int_min, int_max, corr_min, pixel_min,
                        sigma_min, sigma_max, roi_size,
                        filter_size, roi_side, inter_roi):

        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        
        self.corr_min = corr_min
        self.pixel_min = pixel_min
        
        self.int_max = int_max
        self.int_min = int_min
        
        self.roi_size = roi_size
        self.roi_size_1d = int((self.roi_size-1)/2)
        
        self.side_distance = roi_side
        self.roi_distance = inter_roi
        
        if filter_size != self.filter_size:
            self.filter_size = filter_size
            background = medfilt(self.base_frame, kernel_size = self.filter_size)
            background[background ==  0] = np.min(background[background > 0])       
            self.frame = self.base_frame.astype('float') - background
        
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
        
    def int_sigma_limit(self, fitter, return_int, return_sigmas):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            my_roi = self.frame[y-self.roi_size_1d:y+self.roi_size_1d+1, 
                           x-self.roi_size_1d:x+self.roi_size_1d+1]
        
            result, its, success = fitter.fitgaussian(my_roi, roi_index)
            
            if return_sigmas:
                self.sigma_list.append(result[4])
                self.sigma_list.append(result[3])
            elif return_int:
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
        
    def find_particles(self, return_corr):
        
        if return_corr:
            self.corr_min = self.corr_min/2
        
        compare = self.make_gaussian(self.roi_size)
        
        frame_convolution = convolve2d(self.frame, compare, mode='same')
        
        fp = np.ones((3,3), dtype=bool)
        local_peaks = maximum_filter(frame_convolution, footprint=fp)
        local_peaks_bool = (frame_convolution == local_peaks)
        
        max_convolution = np.max(frame_convolution)
        
        beads = (frame_convolution * local_peaks_bool) > self.corr_min*max_convolution
        
        locations = np.transpose(np.where(beads == True))
        
        if return_corr:
            corr = frame_convolution * local_peaks_bool / max_convolution
            
            for roi in locations:
                y = int(roi[0])
                x = int(roi[1]) 
                
                value = corr[y, x]
                self.corr_list.append(value)
                
        return beads, locations
    
    def peak_int_min(self):
        
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)
        
        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1]) 
            
            peak_int = self.frame[y, x]
            
            if peak_int < self.pixel_min:
                keep_boolean[roi_index] = False
                
        self.roi_locations = self.roi_locations[keep_boolean, :] 
        
    def make_corr_list(self):
        pass
        
    
    def main(self, fitter, return_int = False, return_sigmas = False, 
             return_corr = False):
        
        roi_boolean, self.roi_locations = self.find_particles(return_corr)
        
        if return_corr:
            return self.corr_list
        
        self.adjacent_or_boundary_rois_base(roi_boolean)
        
        self.peak_int_min()
                       
        self.int_sigma_limit(fitter, return_int, return_sigmas)
        
        if return_sigmas:
            return self.sigma_list
        elif return_int:
            return self.int_list
        else:
            return self.roi_locations 
        
        
        
        