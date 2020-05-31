# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:20:29 2020

@author: Dion Engels
MBx Python Data Analysis

Frame Analysis

----------------------------

v1.0, ROI detection:  31/05/2020

"""

import numpy as np # for linear algebra
import matlab.engine # to run matlab
import scipy.io as sio # to load in .mat
import time # to sleep to ensure that matlab licence is registered

 #%% Adapted from FindParticles.m from SPectraA4

def ROI_finder(Image, ROI_size):
    
    threshold = [ 0.005 ,0.75 ]
    
    eng = matlab.engine.start_matlab()
    time.sleep(1)
    
    ImageList = Image.tolist()
    ImageMatlab = matlab.double(ImageList)
    
    
    [Corrected, Background] = BackgroundCorrection_WaveletSet(ImageMatlab,4,eng)
    Image = Image - Background
    
    ImageList = Image.tolist()
    ImageMatlab = matlab.double(ImageList)
    
    BeadFilterDict = sio.loadmat('Typical_BeadFilter_for_Finding_AuNPs.mat')
    BeadFilter = BeadFilterDict['BeadFilter']
    BeadFilterList = BeadFilter.tolist()
    BeadFilterMatlab = matlab.double(BeadFilterList)
    
    PopulationDictMatlab = eng.LocateBeadsProfile(ImageMatlab, BeadFilterMatlab, False, threshold[0], threshold[1],20)
    
    PopulationMatlab = PopulationDictMatlab['Location']
    Population = np.array(PopulationMatlab)
    
    eng.quit()
    
    return Population




 #%% Runs BackgroundCorrection_WaveletSet Version 1.2 by Emiel Visser

def BackgroundCorrection_WaveletSet(FrameData,Level,eng):
    
    HiPassMatlab =  eng.BackgroundCorrection_WaveletSet(FrameData,Level)
    
    HiPass = np.array(HiPassMatlab)
    LowPass = FrameData - HiPass
    
    return [HiPass, LowPass]
   