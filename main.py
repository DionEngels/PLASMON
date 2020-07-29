# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:44:02 2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

main

The original main, working without a GUI. Mostly used for development purposes.

----------------------------

v0.1.1, Loading & ROIs & Saving: 31/05/2020
v0.1.2, rainSTORM inspired v1: 03/06/2020
v0.1.3, rainSTORM inspired working v1: 04/06/2020
v0.2, main working for .nd2 loading, no custom ROI fitting: 05/06/2020
v0.2.1. MATLAB loading: 15/06/2020
v0.2.2: own ROI finder: 11/07/2020
v0.2.3: 7x7 and 9x9 ROIs: 13/07/2020
v0.2.4: removed any wavelength dependency
v0.2.5: removed MATLAB ROI finding
v0.2.6: MATLAB v3 loading
v0.2.7: cleanup
v0.2.8: rejection options
v0.3.0: ready for Peter review. MATLAB coordinate system output, bug fix and text output
v0.3.1: different directory for output
v0.3.2: no longer overwrites old data

 """

# GENERAL IMPORTS
from os import mkdir  # to get standard usage
import time  # for timekeeping

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np

# ND2 related
from pims import ND2_Reader, FramesSequenceND  # reader of ND2 files

# loading in .mat
from scipy.io import loadmat

# Own code
import _code.roi_finding as roi_finding
import _code.fitters as fitting
import _code.tools as tools

# %% Inputs
ROI_SIZE = 7  # 7 or 9

# %% Initializations

filenames = ("C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2",)

METHOD = "Gaussian"
DATASET = "MATLAB_v3"  # "MATLAB_v2, "MATLAB_v3" OR "YUYANG"
THRESHOLD_METHOD = "Loose"  # "Strict", "Loose", or "None"

# %% Main loop cell

for name in filenames:
    with ND2_Reader(name) as ND2:
        # create folder for output
        path = name.split(".")[0]
        if DATASET == "MATLAB_v3" or DATASET == "MATLAB_v2":
            path = r"C:\Users\s150127\Downloads\___MBx\datasets\MATLAB"

        directory_try = 0
        directory_success = False
        while not directory_success:
            try:
                mkdir(path)
                directory_success = True
            except:
                directory_try += 1
                if directory_try == 1:
                    path += "_%03d" % directory_try
                else:
                    path = path[:-4]
                    path += "_%03d" % directory_try

        if DATASET == "MATLAB_v2" or DATASET == "MATLAB_v3":
            # Load in MATLAB data

            if DATASET == "MATLAB_v2":
                frames = loadmat('Data_1000f_06_30_pure_matlab_bg_600_v2')['frame']
            elif DATASET == "MATLAB_v3":
                frames = loadmat('Data_1000f_06_30_pure_matlab_bg_600_v3')['frame']

            frames = np.swapaxes(frames, 1, 2)
            frames = np.swapaxes(frames, 0, 1)
            metadata = {'NA': 1, 'calibration_um': 0.120, 'sequence_count': frames.shape[0], 'time_start': 3,
                        'time_start_utc': 3}
            # frames = frames[0:100, :, :]
        elif DATASET == "YUYANG":
            # parse ND2 info
            frames = ND2
            metadata = ND2.metadata
            #  frames = frames[0:5]

        # %% Find ROIs (for standard NP2 file)
        print('Starting to find ROIs')

        fitter = fitting.Gaussian(ROI_SIZE, {}, "None", "Gaussian", 5)

        roi_finder = roi_finding.RoiFinder(frames[0], fitter)
        roi_finder.roi_size = ROI_SIZE
        roi_finder.roi_size_1d = int((ROI_SIZE - 1) / 2)

        if DATASET == "MATLAB_v3":
            roi_finder.int_min = 100
            roi_finder.corr_min = 0.001

        ROI_locations = roi_finder.main(fitter)

        tools.plot_rois(frames[0], ROI_locations, ROI_SIZE)

        # %% Fit Gaussians
        print('Starting to prepare fitting')
        start = time.time()

        thresholds = {'int_min': roi_finder.int_min, 'int_max': roi_finder.int_max,
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

        total_fits = results.shape[0]
        failed_fits = results[np.isnan(results[:, 3]), :].shape[0]
        successful_fits = total_fits - failed_fits

        time_taken = round(time.time() - start, 3)

        print('Time taken: ' + str(time_taken) + ' s. Fits done: ' + str(successful_fits))

        print('Starting saving')
        # %% Plot frames

        # for index, frame in enumerate(frames):
        #     toPlot = results[results[:, 0] == index]
        #     plt.imshow(frames[0], extent=[0,frame.shape[1],frame.shape[0],0], aspect='auto')
        #     #takes x,y hence the switched order
        #     plt.scatter(toPlot[:,3], toPlot[:,2], s=2, c='red', marker='x', alpha=0.5)
        #     plt.title("Frame number: " + str(index))
        #     plt.show()

        if DATASET == "YUYANG":
            # MATLAB works from bottom left as zero point, Python top left. Thus, y is switched
            results[:, 3] = frames[0].shape[0] - results[:, 3]

        # %% filter metadata
        metadata_filtered = {k: v for k, v in metadata.items() if v is not None}
        del metadata_filtered['time_start']
        del metadata_filtered['time_start_utc']

# %% save everything
        tools.save_to_csv_mat('metadata', metadata_filtered, path)
        tools.save_to_csv_mat_roi('ROI_locations', ROI_locations, frames[0].shape[0], path)
        tools.save_to_csv_mat_results('Localizations', results, METHOD, path)

        tools.text_output({}, METHOD, THRESHOLD_METHOD, "", total_fits, failed_fits, time_taken, path)
