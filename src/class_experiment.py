# -*- coding: utf-8 -*-
"""
Created on Thu 01/10/2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

class_experiment

The experiment class of v2 of the program. Holds several datasets.
"""
# GENERAL IMPORTS
from os import mkdir  # to get standard usage
import numpy as np

# OWN CODE

from src.nd2_reading import ND2ReaderSelf
from src.roi_finding import RoiFinder
import src.tt as fitting
from src.hsm import HSMDataset
import src.tools as tools
import src.drift_correction as drift_correction
import src.figure_making as figuring
import src.output as outputting

__self_made__ = True

# %% Experiment


class Experiment:

    def __init__(self, created_by, filename, proceed_question, progress_function):
        self.created_by = created_by
        self.directory = filename
        self.name = None
        self.datasets = []
        self.experiment_settings = None
        self.proceed_question = proceed_question
        self.progress_function = progress_function

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
