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
v1.2: own ROI finder: 11/07/2020
v1.3: 7x7 and 9x9 ROIs: 13/07/2020
v1.4: removed any wavelength dependency
v1.5: removed MATLAB ROI finding
v1.6: MATLAB v3 loading
v1.7: cleanup
v1.8: rejection options

 """

# GENERAL IMPORTS
import os  # to get standard usage
# import math # for generic math
import time  # for timekeeping

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
# import matplotlib.pyplot as plt

# ND2 related
from pims import ND2_Reader  # reader of ND2 files

# comparison to other fitters
import scipy.io

# Own code
import _code.roi_finding as roi_finding
import _code.fitters as fitting
import _code.tools as tools

# %% Inputs
ROI_SIZE = 7  # 7 or 9
FILTER_SIZE = 9

# %% Initializations

filenames = ("C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2",)

METHOD = "Gaussian"
DATASET = "MATLAB_v3"  # "MATLAB_v2, "MATLAB_v3" OR "YUYANG"
THRESHOLD_METHOD = "Loose"  # "Strict", "Loose", or "None"

# %% Main loop cell

for name in filenames:
    with ND2_Reader(name) as ND2:
        # create folder for output
        basedir = os.getcwd()
        directory = name.split(".")[0].split("/")[-1]
        path = os.path.join(basedir, directory)
        try:
            os.mkdir(path)
        except:
            pass

        if DATASET == "MATLAB_v2" or "MATLAB_v3":
            # Load in MATLAB data

            if DATASET == "MATLAB_v2":
                frames = scipy.io.loadmat('Data_1000f_06_30_pure_matlab_bg_600_v2')['frame']
            elif DATASET == "MATLAB_v3":
                frames = scipy.io.loadmat('Data_1000f_06_30_pure_matlab_bg_600_v3')['frame']

            frames = np.swapaxes(frames, 1, 2)
            frames = np.swapaxes(frames, 0, 1)
            metadata = {'NA': 1, 'calibration_um': 0.120, 'sequence_count': frames.shape[0], 'time_start': 3,
                        'time_start_utc': 3}
            # frames = frames[0:100, :, :]
        elif DATASET == "YUYANG":
            # parse ND2 info
            frames = ND2
            metadata = ND2.metadata
            # frames = frames[0:200]

        # %% Find ROIs (for standard NP2 file)
        print('Starting to find ROIs')

        fitter = fitting.Gaussian(ROI_SIZE, {}, "None", "Gaussian", 5)

        roi_finder = roi_finding.RoiFinder(FILTER_SIZE, frames[0], fitter, roi_size=ROI_SIZE)

        if DATASET == "MATLAB_v3":
            roi_finder.int_min = 100
            roi_finder.pixel_min = 100
            roi_finder.corr_min = 0.001

        ROI_locations = roi_finder.main(fitter)

        tools.plot_rois(frames[0], ROI_locations, ROI_SIZE)

        # %% Fit Gaussians
        print('Starting to prepare fitting')
        start = time.time()

        thresholds = {'pixel_min': roi_finder.pixel_min, 'int_min': roi_finder.int_min, 'int_max': roi_finder.int_max,
                      'sigma_min': roi_finder.sigma_min, 'sigma_max': roi_finder.sigma_max}

        if METHOD == "Phasor":
            fitter = fitting.Phasor(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD)
        elif METHOD == "PhasorDumb":
            fitter = fitting.PhasorDumb(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD)
        elif METHOD == "PhasorSum":
            fitter = fitting.PhasorSum(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD)
        elif METHOD == "Gaussian":
            fitter = fitting.Gaussian(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD, 5)
        elif METHOD == "GaussianBackground":
            fitter = fitting.GaussianBackground(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD, 6)

        results = fitter.main(frames, metadata, ROI_locations)

        print('Time taken: ' + str(round(time.time() - start, 3)) + ' s. Fits done: ' + str(results.shape[0]))

        print('Starting saving')
        # %% Plot frames
        # for index, frame in enumerate(frames):
        #     toPlot = results[results[:, 0] == index]
        #     plt.imshow(frames[0], extent=[0,frame.shape[1],frame.shape[0],0], aspect='auto')
        #     #takes x,y hence the switched order
        #     plt.scatter(toPlot[:,3], toPlot[:,2], s=2, c='red', marker='x', alpha=0.5)
        #     plt.title("Frame number: " + str(index))
        #     plt.show()
        # %% filter metadata
        metadata_filtered = {k: v for k, v in metadata.items() if v is not None}
        del metadata_filtered['time_start']
        del metadata_filtered['time_start_utc']

        # ROI_locations dict
        ROI_locations_dict = dict(zip(['x', 'y'], ROI_locations.T))

        # Localization dict
        results_dict = {'Localizations': results}

# %% save everything
        tools.save_to_csv_mat('metadata', metadata_filtered, directory)
        tools.save_to_csv_mat('ROI_locations', ROI_locations_dict, directory)
        tools.save_to_csv_mat('Localizations', results_dict, directory)
