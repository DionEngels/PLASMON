# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:44:02 2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

----------------------------

v0.1, Loading & ROIs & Saving: 31/05/2020
v0.2, rainSTORM inspired v1: 03/06/2020
v0.3, rainSTORM inspired working v1: 04/06/2020
v1.0, main working for .nd2 loading, no custom ROI fitting: 05/06/2020
v1.1. MATLAB loading: 15/06/2020

"""

## GENERAL IMPORTS
import os # to get standard usage
#import math # for generic math
import time # for timekeeping

## Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib.pyplot as plt

## GUI
from tkinter import Tk # for GUI
from tkinter.filedialog import askopenfilenames # for popup that asks to select .nd2's

## ND2 related
from pims import ND2_Reader # reader of ND2 files
#from pims import Frame # for converting ND2 generator to one numpy array

## Multiprocessing
import multiprocessing as mp

## comparison to other fitters
import scipy.io

## Own code
import analysis
import gaussian_fitting
import tools

#%% Inputs
ROI_SIZE = 9 
WAVELENGTH = 637 #nm
THRESHOLD = 5 # X*sigma

#%% Initializations

FILETYPES = [("ND2", ".nd2")]

# Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
# filenames = askopenfilenames(filetypes = FILETYPES,title = "Select file", initialdir =os.getcwd()) # show an "Open" dialog of only .nd2's starting at code location

filenames = ("C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2",)

METHOD = "ScipyPhasorGuess"
DATASET = "MATLAB" # "MATLAB" OR "YUYANG"
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
        
        if DATASET == "MATLAB":
            ## Load in MATLAB data
            
            frames = scipy.io.loadmat('Data_1000f_14_06_pure_matlab_bg_600')['frame']
            frames = np.swapaxes(frames,1,2)
            frames = np.swapaxes(frames,0,1)
            metadata = {'NA' : 1, 'calibration_um' : 0.2, 'sequence_count' : frames.shape[0], 'time_start' : 3, 'time_start_utc': 3}
            #frames = frames[0:100,:,:]
        elif DATASET == "YUYANG":
            ## parse ND2 info
            frames = ND2
            metadata = ND2.metadata
            frames = frames[0:10]

        #%% Find ROIs (for standard NP2 file)
        print('Starting to find ROIs')
        
        if DATASET == "MATLAB":
            for i in range(20):
                for j in range(10):
                    if i == 0 and j == 0:
                        ROI_locations = [19, 19]
                    else:
                        ROI_locations = np.vstack((ROI_locations, [19+j*20, 19+i*20]))
        
        elif DATASET == "YUYANG":
            #ROI_locations = analysis.ROI_finder(frames[0],ROI_size)
            ROI_locations = np.load('ROI_locations.npy')

            ROI_locations = ROI_locations - 1
            #ROI_locations = ROI_locations[0:2,:]

            ## switch array columns since MATLAB gives x,y. Python likes y,x
            ROI_locations = tools.switch(ROI_locations)

        plt.imshow(frames[0], extent=[0,frames[0].shape[1],frames[0].shape[0],0], aspect='auto')
        #takes x,y hence the switched order
        plt.scatter(ROI_locations[:,1], ROI_locations[:,0], s=2, c='red', marker='x', alpha=0.5)
        plt.title("ROI locations")
        plt.show()


        #%% Fit Gaussians
        print('Starting to prepare fitting')
        start = time.time()

        if METHOD == "Sum":
            summation = gaussian_fitting.summation(metadata, ROI_SIZE, WAVELENGTH, THRESHOLD, ROI_locations, METHOD)
            results = summation.main(frames,metadata)
        elif METHOD == "rainSTORM":
            rainSTORM = gaussian_fitting.rainSTORM_Dion(metadata, ROI_SIZE, WAVELENGTH, THRESHOLD, ROI_locations, METHOD)
            results = rainSTORM.main(frames, metadata)
        elif METHOD == "PhasorOnly":
            phasor_only = gaussian_fitting.phasor_only(metadata, ROI_SIZE, WAVELENGTH, THRESHOLD, ROI_locations, METHOD)
            results = phasor_only.main(frames, metadata)
        elif METHOD == "ScipyPhasorGuess":
            scipy_phasor = gaussian_fitting.scipy_phasor(metadata, ROI_SIZE, WAVELENGTH, THRESHOLD, ROI_locations, METHOD)
            results = scipy_phasor.main(frames, metadata)
        elif METHOD == "ScipyPhasorGuessROILoop":
            scipy_phasor_guess_roi = gaussian_fitting.scipy_phasor_guess_roi(metadata, ROI_SIZE, WAVELENGTH, THRESHOLD, ROI_locations, METHOD)
            results = scipy_phasor_guess_roi.main(frames, metadata)
        elif METHOD == "ScipyLastFitGuessROILoop":
            scipy_last_fit_guess_roi = gaussian_fitting.scipy_last_fit_guess_roi(metadata, ROI_SIZE, WAVELENGTH, THRESHOLD, ROI_locations, METHOD)
            results = scipy_last_fit_guess_roi.main(frames, metadata)
        
        

        print('Time taken: ' + str(round(time.time() - start, 3)) + ' s. Fits done: ' + str(results.shape[0]))




        print('Starting saving')
        #%% Plot frames
        # for index, frame in enumerate(frames):
        #     toPlot = results[results[:, 0] == index]
        #     plt.imshow(frames[0], extent=[0,frame.shape[1],frame.shape[0],0], aspect='auto')
        #     #takes x,y hence the switched order
        #     plt.scatter(toPlot[:,3], toPlot[:,2], s=2, c='red', marker='x', alpha=0.5)
        #     plt.title("Frame number: " + str(index))
        #     plt.show()
        #%% filter metadata
        metadata_filtered = {k:v for k, v in metadata.items() if v is not None}
        del metadata_filtered['time_start']
        del metadata_filtered['time_start_utc']

        ## ROI_locations dict
        ROI_locations_dict = dict(zip(['x', 'y'], ROI_locations.T))

        ## Localization dict
        results_dict = {'Localizations': results}

#%% save everything
        tools.save_to_csv_mat('metadata', metadata_filtered, directory)
        tools.save_to_csv_mat('ROI_locations', ROI_locations_dict, directory)
        tools.save_to_csv_mat('Localizations', results_dict, directory)
