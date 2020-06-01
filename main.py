# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:44:02 2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

----------------------------

v0.1, Loading & ROIs & Saving: 31/05/2020
 
"""

## GENERAL IMPORTS
import os # to get standard usage
import math # for generic math

## Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib.pyplot as plt

## GUI
from tkinter import Tk # for GUI
from tkinter.filedialog import askopenfilenames # for popup that asks to select .nd2's

## ND2 related
from pims import ND2_Reader # reader of ND2 files



## Own code
import analysis
import gaussian_fitting
import tools

#%% Initialization
ROI_size = 9
wavelength = 637 #nm
threshold = 5 # X*sigma

filetypes=[("ND2", ".nd2")]

# Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
# filenames = askopenfilenames(filetypes = filetypes,title = "Select file", initialdir =os.getcwd()) # show an "Open" dialog of only .nd2's starting at code location

filenames = ("C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2",)

#%% Main loop cell
            
for name in filenames:
    with ND2_Reader(name) as ND2:
        ## create folder for output
        basedir = os.getcwd()
        directory = name.split(".")[0].split("/")[-1]
        path = os.path.join(basedir, directory)
        try:
            os.mkdir(path)
        except:
            pass
        
        ## parse ND2 info
        frames = ND2
        metadata = ND2.metadata
        
        ## Find ROIs
        
        #ROI_locations = analysis.ROI_finder(frames[0],ROI_size)
        ROI_locations = np.array([[1, 2], [3, 4]])

        
        ## Fit Gaussians
        
        rainSTORM = gaussian_fitting.rainSTORM_Dion(metadata, ROI_size, wavelength)
        rainSTORM.determine_threshold(frames[0], ROI_locations, threshold)
        rainSTORM.main(ROI_locations,frames)
        
        
        
        
        ## Plot frames
        # for index, frame in enumerate(frames):
        #     plt.matshow(frame,fignum=index)
        #     plt.title("Frame number: " + str(index))
        #     plt.show()

        
        ## filter metadata
        metadata_filtered = {k:v for k,v in metadata.items() if v is not None}
        del metadata_filtered['time_start']
        del metadata_filtered['time_start_utc']
        
        ## ROI_locations dict
        ROI_locations_dict = dict(zip(['x','y'], ROI_locations.T))
        
#%% save everything
        tools.SaveToCsvMat('metadata', metadata_filtered, directory)
        tools.SaveToCsvMat('ROI_locations', ROI_locations_dict, directory)
        
    

