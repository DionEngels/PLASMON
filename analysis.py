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

class roi_finder():
    
    def __init__(self, intensity_min, intensity_max, sigma_min, sigma_max, ):
        self.intensity_min = intensity_min
        self.intensity_max = intensity_max
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        self.shape = 0.8
        
        
    def main(self, frame):
        
        pass
        
        