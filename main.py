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
from scipy.ndimage import median_filter

# loading in .mat
from scipy.io import loadmat

# Own code
from src.roi_finding import RoiFinder
import src.fitters as fitting
import src.tools as tools
import src.drift_correction as drift_correction
import src.figure_making as figuring
import src.output as outputting
import src.nd2_reading as nd2_reading
import src.hsm as hsm

from src.hsm import normxcorr2, normxcorr2_large

__self_made__ = True

# %% Inputs
ROI_SIZE = 7  # 7 or 9

# %% Initializations

filenames = ("C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2",)
hsm_dir = ("C:/Users/s150127/Downloads/___MBx/datasets/_1nMimager_newGNRs_100mW_HSM",)

fit_options = ["Gaussian - Fit bg", "Gaussian - Estimate bg",
               "Phasor + Intensity", "Phasor + Sum", "Phasor"]

FIGURE_OPTION = "Overview"  # "Overview" "All"

METHOD = "Gaussian - Estimate bg"
DATASET = "YUYANG"  # "MATLAB_v2, "MATLAB_v3" OR "YUYANG"
THRESHOLD_METHOD = "Loose"  # "Strict", "Loose", or "None"
CORRECTION = "SN_objTIRF_PFS_510-800"  # "Matej_670-890"
NM_OR_PIXELS = "nm"
SLICE = slice(0, 10)

# %% Classes


class Dataset:
    def __init__(self, experiment):
        self.experiment = experiment
        self.frames = None
        self.frame_for_rois = None
        self.metadata = None
        self.fitter = None
        self.drift_corrector = None
        self.roi_offset = None
        self.active_rois = []

    @staticmethod
    def parse_start_end(start, end):
        if start == "Leave empty for start" and end == "Leave empty for end":
            return slice(None), 0
        elif start == "Leave empty for start" and end != "Leave empty for end":
            return slice(0, int(end)), 0
        elif start != "Leave empty for start" and end == "Leave empty for end":
            return slice(int(start), None), start
        else:  # start != "Leave empty for start" and end != "Leave empty for end":
            return slice(int(start), int(end)), start

    @staticmethod
    def correlate_frames(frame_old, frame_new):
            if frame_old.shape == frame_new.shape:
                corr = normxcorr2(frame_old, frame_new)
                maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
                offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
            else:
                corr = normxcorr2_large(frame_old, frame_new)
                maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
                offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
            return offset

    def correlate(self, settings):
        x_slice, x_offset = self.parse_start_end(settings['x_min'], settings['x_max'])
        y_slice, y_offset = self.parse_start_end(settings['y_min'], settings['y_max'])
        offset_crop = np.asarray([y_offset, x_offset])

        experiment_frame_shape = self.experiment.frame_for_rois.shape
        frame_shape = self.frame_for_rois.shape

        # test offset crop
        if frame_shape[0] > experiment_frame_shape[0] and frame_shape[1] > experiment_frame_shape[1]:
            small_frame = self.experiment.frame_for_rois
            cropped_frame = self.frame_for_rois(y_slice, x_slice)
            offset = self.correlate_frames(small_frame, cropped_frame) - offset_crop
        elif frame_shape[0] < experiment_frame_shape[0] and frame_shape[1] < experiment_frame_shape[1]:
            small_frame = self.frame_for_rois
            cropped_frame = self.experiment.frame_for_rois(y_slice, x_slice)
            offset = self.correlate_frames(cropped_frame, small_frame) + offset_crop
        else:
            old_frame = self.experiment.frame_for_rois
            new_frame = self.frame_for_rois
            offset = self.correlate_frames(old_frame, new_frame)
        return offset

    def find_rois(self, settings, frame_for_rois, created_by):
        self.roi_offset = self.correlate(settings)
        self.active_rois = [roi for roi in self.experiment.rois if roi.in_frame(self.frame_for_rois.shape,
                                                                                self.roi_offset)]


class TimeTrace(Dataset):
    def __init__(self, experiment, nd2):
        super().__init__(experiment)
        self.frames = nd2
        background = median_filter(np.asarray(nd2[0]), size=9)
        self.frame_for_rois = np.asarray(nd2[0]).astype('float') - background
        self.metadata = nd2.get_metadata()
        self.pixels_or_nm = None
        self.slice = None

    def prepare_run(self, settings):
        check = self.experiment.proceed_question("OK", "Cancel", "Are you sure?",
                                          "Fitting may take a while. Are you sure everything is set up correctly?")
        if not check:
            return

        if settings['#cores'] > 1 and "Phasor" in settings['method']:
            cores_check = self.experiment.proceed_question("ok_cancel", "Just a heads up",
                                                           """Phasor will be used with one core since the
            overhead only slows it down""")
            if not cores_check:
                return
            settings['#cores'] = 1

        # TO WRITE FUNCTION FOR
        max_its = 300
        self.pixels_or_nm = settings['pixels_or_nm']
        self.slice = self.parse_start_end(settings['frame_begin'], settings['frame_end'])

        if settings['method'] == "Phasor + Intensity":
            self.fitter = fitting.Phasor(settings)
        elif settings['method'] == "Phasor":
            self.fitter = fitting.PhasorDumb(settings)
        elif settings['method'] == "Gaussian - Fit bg":
            self.fitter = fitting.GaussianBackground(settings, max_its, 6)
        elif settings['method'] == "Gaussian - Estimate bg":
            self.fitter = fitting.Gaussian(settings, max_its, 5)
        else:
            self.fitter = fitting.PhasorSum(settings)

    def run(self):
        self.fitter.run(self)


class HSMDataset(Dataset):
    def __init__(self, experiment, nd2):
        super().__init__(experiment)
        self.frames = nd2
        self.metadata = nd2.get_metadata()

        self.corrected_merged, self.corrected = self.hsm_drift()
        self.frame_for_rois = self.corrected_merged


    def prepare_run(self, settings):
        self.hsm_object.change_settings(settings)

    def run(self):
        self.hsm_object.run(self)

    def hsm_drift(self):
        x =1
        y= 1
        return x, y


class Experiment:

    def __init__(self, created_by, filename, proceed_question):
        self.created_by = created_by
        self.directory = filename
        self.name = None
        self.datasets = []
        self.experiment_settings = None
        self.proceed_question = proceed_question

        nd2 = nd2_reading.ND2ReaderSelf(filename)

        if created_by == 'HSM':
            self.init_new_hsm(nd2)
            self.frame_for_rois = np.asarray(self.datasets[-1].corrected_merged)
        elif created_by == 'TT':
            self.init_new_tt(nd2)
            self.frame_for_rois = np.asarray(nd2[0])
        self.roi_finder = RoiFinder(self.frame_for_rois)
        self.rois = self.roi_finder.main()

    def init_new_hsm(self, nd2):
        hsm_object = HSMDataset(self, nd2)
        self.datasets.append(hsm_object)

    def init_new_tt(self, nd2):
        time_trace_object = TimeTrace(self, nd2)
        self.datasets.append(time_trace_object)

    def change_rois(self, settings):
        self.roi_finder.change_settings(settings)
        self.rois = self.roi_finder.main()

    def show_rois(self, experiment_or_dataset):
        if experiment_or_dataset == "Experiment":
            figuring.plot_rois(self.frame_for_rois, self.rois, self.roi_finder.roi_size)
        elif experiment_or_dataset == "Dataset":
            figuring.plot_rois(self.datasets[-1].frame_for_rois, self.datasets[-1].active_rois,
                               self.roi_finder.roi_size)

    def finalize_rois(self, name, experiment_settings):
        self.name = name
        file_dir = '/'.join(self.directory.split(".")[0].split("/")[:-1]) + '/'

        date = self.datasets[0].metadata.pop('date', None)
        if date is None:
            date = "XXXX-XX-XX"
        else:
            date_split = date.split(" ")[0].split("-")
            date = "{:04d}-{:02d}-{:02d}".format(int(date_split[2]), int(date_split[1]), int(date_split[0]))

        file_dir = file_dir + date + "_" + name

        directory_try = 0
        directory_success = False
        while not directory_success:
            try:
                mkdir(file_dir)
                directory_success = True
            except:
                directory_try += 1
                if directory_try == 1:
                    file_dir += "_%03d" % directory_try
                else:
                    file_dir = file_dir[:-4]
                    file_dir += "_%03d" % directory_try
        self.directory = file_dir
        self.experiment_settings = experiment_settings

    def find_rois_dataset(self, settings):
        self.datasets[-1].find_rois(settings, self.frame_for_rois, self.created_by)

    def add_to_queue(self, settings):
        self.datasets[-1].prepare_run(settings)

    def save(self):

        tools.convert_to_matlab(self.rois)

        results = self.rois.to_dict()

        outputting.save_results(self.dir, results)
        outputting.save_roi_pos(self.dir, self.rois)
        outputting.save_datasets(self.dir, self.datasets)

        figuring.save_overview(self.dir, self.rois)
        figuring.individual_figures(self.dir, self.rois, self.datasets)

    def run(self):

        for dataset in self.datasets:
            dataset.run()

        self.save()


# %% GUI-less specific

def proceed_question(option1, option2, title, text):
    answer = input(title + "\n" + text + "\n" + option1 + "/" + option2)
    if answer == option1:
        return True
    else:
        return False


# %% Main loop cell

for name in filenames:

    experiment = Experiment("TT", name, proceed_question)

    experiment.show_rois("Experiment")

    defaults = experiment.roi_finder.get_settings()

    settings_rois = {'int_max': np.inf, 'int_min': 0,
                     'sigma_min': 0, 'sigma_max': int((ROI_SIZE - 1) / 2),
                     'corr_min': 0.05, 'roi_size': ROI_SIZE, 'filter_size': 9,
                     'roi_side': 11, 'inter_roi': 9}

    experiment.change_rois(settings_rois)

    experiment.show_rois("Experiment")

    name = "v2_test"

    settings_experiment = {'All Figures': True}

    experiment.finalize_rois(name, settings_experiment)

    settings_correlation = {'x_min': "Leave empty for start", 'x_max': "Leave empty for end",
                            'y_min': "Leave empty for start", 'y_max': "Leave empty for end"}

    experiment.find_rois_dataset(settings_correlation)

    experiment.show_rois("Dataset")

    settings_runtime = {'method': "Gaussian", 'rejection': "Loose", '#cores': 6, "pixels_or_nm": "nm",
                        'roi_size': ROI_SIZE,
                        'frame_begin': "Leave empty for start", 'frame_end': "Leave empty for end"}
    experiment.add_to_queue(settings_runtime)

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

        # roi_finder.corr_min = 0.2
        # roi_finder.int_min = 4000

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

        hsm = hsm.HSM(hsm_dir, np.asarray(frames[0], dtype=frames[0].dtype), ROI_locations.copy(), metadata, CORRECTION)
        hsm_result, hsm_raw, hsm_intensity = hsm.main(verbose=False)

        print('Starting saving')

        # %% Figures

        settings = {'roi_size': ROI_SIZE}

        start = time.time()

        figuring.save_graphs(frames, results, results_drift, ROI_locations, METHOD, "pixels", FIGURE_OPTION,
                             path, event_or_not, settings, metadata['timesteps'][SLICE],
                             hsm_result, hsm_raw, hsm_intensity, hsm.wavelength)

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

            hsm_result, hsm_raw, hsm_intensity = tools.switch_to_matlab_hsm(hsm_result, hsm_raw, hsm_intensity)

        # %% save everything
        outputting.save_to_csv_mat_metadata('metadata', metadata, path)
        outputting.save_to_csv_mat_roi('ROI_locations', ROI_locations, path)
        outputting.save_to_csv_mat_drift('Drift_correction', drift, path)
        outputting.save_to_csv_mat_results('Localizations', results, METHOD, path)
        outputting.save_to_csv_mat_results('Localizations_drift', results_drift, METHOD, path)
        outputting.save_hsm(hsm_result, hsm_raw, hsm_intensity, path)

        outputting.text_output({}, METHOD, THRESHOLD_METHOD, "", total_fits, failed_fits, time_taken, path)
