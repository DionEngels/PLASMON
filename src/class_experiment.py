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

from warnings import warn

# OWN CODE

from src.nd2_reading import ND2ReaderSelf
from src.roi_finding import RoiFinder
import src.tt as fitting
from src.hsm import HSMDataset
import src.tools as tools
import src.figure_making as figuring
import src.output as outputting

__self_made__ = True

# %% Experiment


class Experiment:

    def __init__(self, created_by, filename, proceed_question, error_func, progress_updater, show_rois):
        self.created_by = created_by
        self.directory = filename
        self.name = None
        self.datasets = []
        self.settings = None
        self.proceed_question = proceed_question
        self.progress_updater = progress_updater
        self.error_func = error_func
        self.show_rois_func = show_rois

        if created_by == 'HSM':
            self.init_new_hsm(filename)
        elif created_by == 'TT':
            self.init_new_tt(filename)
        self.frame_for_rois = self.datasets[-1].frame_for_rois
        self.roi_finder = RoiFinder(self.frame_for_rois)
        self.rois = self.roi_finder.main()

    def init_new_hsm(self, filename):
        nd2 = ND2ReaderSelf(filename)
        hsm_object = HSMDataset(self, nd2, filename)
        self.datasets.append(hsm_object)

    def init_new_tt(self, filename):
        nd2 = ND2ReaderSelf(filename)
        time_trace_object = fitting.TimeTrace(self, nd2, filename)
        self.datasets.append(time_trace_object)

    def change_rois(self, settings):
        self.roi_finder.change_settings(settings)
        self.rois = self.roi_finder.main()

    def show_rois(self, experiment_or_dataset, figure=None):
        if experiment_or_dataset == "Experiment":
            self.show_rois_func(self.frame_for_rois, roi_locations=self.rois,
                                roi_size=self.roi_finder.roi_size, figure=figure)
        elif experiment_or_dataset == "Dataset":
            self.show_rois_func(self.datasets[-1].frame_for_rois,
                                roi_locations=self.datasets[-1].active_rois, figure=figure,
                                roi_size=self.roi_finder.roi_size, roi_offset=self.datasets[-1].roi_offset)

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
        self.settings = experiment_settings

    def find_rois_dataset(self, settings):
        self.datasets[-1].find_rois(settings)

    def add_to_queue(self, settings):
        status = self.datasets[-1].prepare_run(settings)
        if status is False:
            return status

    def run(self):

        for dataset in self.datasets:
            self.progress_updater.new_dataset(dataset.type)
            dataset.run()

            # clear memory
            try:
                dataset.frames.close()
            except:
                pass
            dataset.frames = None

        self.save()

    def save(self):
        self.progress_updater.message("Starting saving")
        settings = self.settings_to_dict()
        outputting.save_settings(self.directory, settings)

        try:
            self.progress_updater.message("Saving overview")
            figuring.save_overview(self)
        except:
            warn("Overview figure creation failed", RuntimeWarning)

        if self.settings['All Figures'] is True:
            try:
                self.progress_updater.message("Saving individual figures")
                figuring.individual_figures(self)
            except:
                warn("Individual figure creation failed", RuntimeWarning)

        self.progress_updater.message("Converting to MATLAB coordinate system")
        tools.convert_to_matlab(self)
        results = self.rois_to_dict()
        metadata = self.metadata_to_dict()
        self.progress_updater.message("Saving to .mat")
        outputting.save_to_mat(self.directory, "Results", results)
        outputting.save_to_mat(self.directory, "Metadata", metadata)
        self.progress_updater.message("Done")


    def rois_to_dict(self):
        result_dict = {}
        for roi in self.rois:
            result_dict["ROI_{}".format(roi.index)] = roi.results

        result_dict = {'Results': result_dict}
        return result_dict

    def metadata_to_dict(self):
        metadata_dict = {}
        for dataset in self.datasets:
            metadata_dict['meta_{}'.format(dataset.name)] = dataset.metadata

        metadata_dict = {'Metadata': metadata_dict}
        return metadata_dict

    def settings_to_dict(self):
        settings_dict = {'Experiment': self.settings, 'ROIs': self.roi_finder.get_settings()}

        for dataset in self.datasets:
            settings_dict[dataset.name] = dataset.settings
            # add type, filename to the front
            settings_dict[dataset.name] = {**{'Type': dataset.type}, **{'filename': dataset.filename},
                                           **settings_dict[dataset.name]}
            settings_dict[dataset.name]['Offset'] = dataset.roi_offset

        return settings_dict
