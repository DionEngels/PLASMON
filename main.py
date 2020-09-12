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
v0.4: drift correction v1: 31/07/2020
v1.0: bugfixes and release: 07/08/2020

 """

# GENERAL IMPORTS
from os import mkdir  # to get standard usage
import time  # for timekeeping

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np

# loading in .mat
from scipy.io import loadmat

# Own code
import _code.roi_finding as roi_finding
import _code.fitters as fitting
import _code.tools as tools
import _code.drift_correction as drift_correction
import _code.figure_making as figuring
import _code.output as outputting
import _code.nd2_reading as nd2_reading
import _code.hsm as hsm

# %% Inputs
ROI_SIZE = 7  # 7 or 9

# %% Initializations

filenames = ("C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2",)
hsm_dir = ("C:/Users/s150127/Downloads/___MBx/datasets/_1nMimager_newGNRs_100mW_HSM",)

fit_options = ["Gaussian - Fit bg", "Gaussian - Estimate bg",
               "Phasor + Intensity", "Phasor + Sum", "Phasor"]

FIGURE_OPTION = "Overview" # "Overview" "All"

METHOD = "Gaussian - Estimate bg"
DATASET = "YUYANG"  # "MATLAB_v2, "MATLAB_v3" OR "YUYANG"
THRESHOLD_METHOD = "Loose"  # "Strict", "Loose", or "None"
CORRECTION = "SN_objTIRF_PFS_510-800"  # "Matej_670-890"
NM_OR_PIXELS = "nm"
SLICE = slice(0, 10)

# %% Main loop cell

for name in filenames:
    with nd2_reading.ND2ReaderSelf(name) as ND2:
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
            #  frames = frames[SLICE, :, :]
            n_frames = frames.shape[0]
        elif DATASET == "YUYANG":
            # parse ND2 info
            frames = ND2
            metadata = ND2.get_metadata()
            frames = frames[SLICE]
            n_frames = len(frames)

        # %% Find ROIs (for standard NP2 file)
        print('Starting to find ROIs')

        fitter = fitting.Gaussian(ROI_SIZE, {}, "None", "Gaussian", 5, 300)

        roi_finder = roi_finding.RoiFinder(frames[0], fitter)
        roi_finder.roi_size = ROI_SIZE
        roi_finder.roi_size_1d = int((ROI_SIZE - 1) / 2)

        if DATASET == "MATLAB_v3":
            roi_finder.int_min = 100
            roi_finder.corr_min = 0.001

        roi_finder.corr_min = 0.2
        roi_finder.int_min = 4000

        ROI_locations = roi_finder.main(fitter)
        max_its = roi_finder.find_snr(fitter)

        figuring.plot_rois(frames[0], ROI_locations, ROI_SIZE)

        # %% Fit Gaussians
        print('Starting to prepare fitting')
        start = time.time()

        thresholds = {'int_min': roi_finder.int_min, 'int_max': roi_finder.int_max,
                      'sigma_min': roi_finder.sigma_min, 'sigma_max': roi_finder.sigma_max}

        if METHOD == "Phasor + Intensity":
            fitter = fitting.Phasor(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD)
        elif METHOD == "Phasor":
            fitter = fitting.PhasorDumb(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD)
        elif METHOD == "Phasor + Sum":
            fitter = fitting.PhasorSum(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD)
        elif METHOD == "Gaussian - Estimate bg":
            fitter = fitting.Gaussian(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD, 5, max_its)
        elif METHOD == "Gaussian - Fit bg":
            fitter = fitting.GaussianBackground(ROI_SIZE, thresholds, THRESHOLD_METHOD, METHOD, 6, max_its)

        results = fitter.main(frames, metadata, ROI_locations)

        total_fits = results.shape[0]
        failed_fits = results[np.isnan(results[:, 3]), :].shape[0]
        successful_fits = total_fits - failed_fits

        time_taken = round(time.time() - start, 3)

        print('Time taken: ' + str(time_taken) + ' s. Fits done: ' + str(successful_fits))

        # %% Plot frames

        # for index, frame in enumerate(frames):
        #     toPlot = results[results[:, 0] == index]
        #     plt.imshow(frames[0], extent=[0,frame.shape[1],frame.shape[0],0], aspect='auto')
        #     #takes x,y hence the switched order
        #     plt.scatter(toPlot[:,3], toPlot[:,2], s=2, c='red', marker='x', alpha=0.5)
        #     plt.title("Frame number: " + str(index))
        #     plt.show()

        # %% drift correction

        print('Starting drift correction')
        drift_corrector = drift_correction.DriftCorrector(METHOD)
        results_drift, drift, event_or_not = drift_corrector.main(results, ROI_locations, n_frames)

        # %% HSM

        print('Starting HSM')

        hsm = hsm.HSM(hsm_dir, np.asarray(frames[0], dtype=frames[0].dtype), ROI_locations, metadata, CORRECTION)
        hsm_result, hsm_intensity = hsm.main(verbose=True)

        print('Starting saving')

        # %% Figures

        settings = {'roi_size': ROI_SIZE}

        start = time.time()

        figuring.save_graphs(frames, results, results_drift, ROI_locations, METHOD, "pixels", FIGURE_OPTION,
                             path, event_or_not, settings, metadata['timesteps'][SLICE])

        time_taken = round(time.time() - start, 3)
        print('Time taken plotting: ' + str(time_taken) + ' s. Fits done: ' + str(successful_fits))

        # %% Convert coordinate system

        if DATASET == "YUYANG":
            results = tools.switch_results_to_matlab_coordinates(results, frames[0].shape[0], METHOD,
                                                                 NM_OR_PIXELS, metadata)
            results_drift = tools.switch_results_to_matlab_coordinates(results_drift, frames[0].shape[0], METHOD,
                                                                       NM_OR_PIXELS, metadata)
            ROI_locations = tools.switch_axis_to_matlab_coordinates(ROI_locations, frames[0].shape[0])
            drift = tools.switch_axis(drift)

        # %% save everything
        outputting.save_to_csv_mat_metadata('metadata', metadata, path)
        outputting.save_to_csv_mat_roi('ROI_locations', ROI_locations, frames[0].shape[0], path)
        outputting.save_to_csv_mat_drift('Drift_correction', drift, path)
        outputting.save_to_csv_mat_results('Localizations', results, METHOD, path)
        outputting.save_to_csv_mat_results('Localizations_drift', results_drift, METHOD, path)

        outputting.text_output({}, METHOD, THRESHOLD_METHOD, "", total_fits, failed_fits, time_taken, path)
