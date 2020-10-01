# -*- coding: utf-8 -*-
"""
Created on Thu 01/10/2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

data

The data structure of v2 of the program
"""
# GENERAL IMPORTS
from os import mkdir  # to get standard usage
import numpy as np

# own code
from src.hsm import normxcorr2, normxcorr2_large  # for correlation
from src.nd2_reading import ND2ReaderSelf
from src.roi_finding import RoiFinder
import src.fitters as fitting
from src.hsm import HSMDataset
import src.tools as tools
import src.drift_correction as drift_correction
import src.figure_making as figuring
import src.output as outputting

# %% ROI


class Roi:

    def __init__(self, x, y):

        self.x = x
        self.y = y
        self.index = None

        self.results = {}

    def set_index(self, index):
        self.index = index

    def get_roi(self, frame, roi_size_1d):
        return frame[self.y - roi_size_1d:self.y + roi_size_1d + 1,
                     self.x - roi_size_1d:self.x + roi_size_1d + 1]

    def in_frame(self, shape, offset):
        if self.x + offset[1] < 0 or self.x + offset[1] > shape[1]:
            in_frame_boolean = False
        elif self.y + offset[0] < 0 or self.y + offset[0] > shape[0]:
            in_frame_boolean = False
        else:
            in_frame_boolean = True

        return in_frame_boolean
# %% Dataset


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
# %% Experiment


class Experiment:

    def __init__(self, created_by, filename, proceed_question):
        self.created_by = created_by
        self.directory = filename
        self.name = None
        self.datasets = []
        self.experiment_settings = None
        self.proceed_question = proceed_question

        nd2 = ND2ReaderSelf(filename)

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
        time_trace_object = fitting.TimeTrace(self, nd2)
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