# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:44:02 2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

----------------------------

v0.1, Loading & ROIs & Saving: 31/05/2020
v0.2, rainSTORM inspired v1: 03/06/2020
 
"""

## GENERAL IMPORTS
import os # to get standard usage
import math # for generic math
import time # for timekeeping

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
        frames = frames[0:10]
        ## Find ROIs
        
        print('Starting to find ROIs')
        
        #ROI_locations = analysis.ROI_finder(frames[0],ROI_size)
        ROI_locations = np.load('ROI_locations.npy')
        #ROI_locations = np.array([[1, 2], [3, 4]])
        
        plt.matshow(frames[0],fignum=0)
        plt.scatter(ROI_locations[:,0],ROI_locations[:,1], s = 2,c = 'red', marker='x', alpha=0.5)
        plt.title("ROI locations")
        plt.show()

        
        ## Fit Gaussians
        print('Starting to prepare fitting')
        rainSTORM = gaussian_fitting.rainSTORM_Dion(metadata, ROI_size, wavelength)
        rainSTORM.determine_threshold(frames[0], ROI_locations, threshold)
        print('Starting fitting')
        start = time.time()
        results = rainSTORM.main(ROI_locations,frames, metadata)
        print('Time taken: ' + str(round(time.time() - start,3)) + ' s. Fits done: ' + str(results.shape[0]))
        
        print('Starting saving')
        
        #%% Plot frames
        for index, frame in enumerate(frames):
            toPlot = results[results[:,0] == index]
            plt.matshow(frame,fignum=index)
            plt.scatter(toPlot[:,2],toPlot[:,3], s = 2,c = 'red', marker='x', alpha=0.5)
            plt.title("Frame number: " + str(index))
            plt.show()

        
        #%% filter metadata
        metadata_filtered = {k:v for k,v in metadata.items() if v is not None}
        del metadata_filtered['time_start']
        del metadata_filtered['time_start_utc']
        
        ## ROI_locations dict
        ROI_locations_dict = dict(zip(['x','y'], ROI_locations.T))
        
        ## Localization dict
        results_dict = {'Localizations': results}
        
#%% save everything
        tools.SaveToCsvMat('metadata', metadata_filtered, directory)
        tools.SaveToCsvMat('ROI_locations', ROI_locations_dict, directory)
        tools.SaveToCsvMat('Localications', results_dict,directory)
        
    

