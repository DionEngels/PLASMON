# -*- coding: utf-8 -*-
"""
Created on Thu 01/10/2020

----------------------------

@author: Dion Engels
PLASMON Data Analysis

class_experiment

The experiment class of v2 of the program. Holds several datasets.

-----------------

v2.0: part of v2.0: 15/10/2020

"""
# GENERAL IMPORTS
from os import mkdir  # to get standard usage
import logging  # for logging warnings
logger = logging.getLogger('main')
import time  # for time keeping

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
    """
    Experiment class. Holds datasets of same experiment and same ROIs.
    """
    def __init__(self, created_by, filename, proceed_question, error_func, progress_updater, show_rois):
        """
        Initialises experiment. Sets some settings and calls first dataset initialization and ROI finder
        ----------------------
        :param created_by: whether or not first dataset is TT or HSM
        :param filename: filename of first nd2
        :param proceed_question: proceed question function. Changes if GUI is used or not
        :param error_func: Error function. Also changes if GUI is used or not
        :param progress_updater: Progress updater. GUI changes this
        :param show_rois: Function to show ROIs. Also changes with GUI
        """
        self.created_by = created_by
        self.directory = filename
        self.dir_made = False
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
        self.roi_finder = RoiFinder(self.frame_for_rois, self.datasets[-1].data_type_signed)
        self.rois = self.roi_finder.main()

    def init_new_hsm(self, filename):
        """
        Add a new HSM to experiment. Loads nd2, initialises HSM class, appends to self.datasets
        -----------------------------------
        :param filename: filename of new HSM
        :return: None. Edits class.
        """
        nd2 = ND2ReaderSelf(filename)
        hsm_object = HSMDataset(self, nd2, filename)
        self.datasets.append(hsm_object)

    def init_new_tt(self, filename):
        """
        Add a new TT to experiment. Loads nd2, initialises TT class, appends to self.datasets
        :param filename: filename of new HSM
        :return: None. Edits class.
        """
        nd2 = ND2ReaderSelf(filename)
        time_trace_object = fitting.TimeTrace(self, nd2, filename)
        self.datasets.append(time_trace_object)

    def change_rois(self, settings):
        """
        Changes settings of ROI finder, also changing the ROIs of the experiment.
        ------------------------------
        :param settings: new settings of ROI finder
        :return: None. Changes class and ROIs in class
        """
        self.roi_finder.change_settings(settings)
        self.rois = self.roi_finder.main()

    def show_rois(self, experiment_or_dataset, figure=None, overwrite=False):
        """
        Show ROIs function. Depending if you want to see experiment or dataset, shows either using show_rois func
        ----------------------------------------------------
        :param experiment_or_dataset: String determines if Experiment or Dataset is shown
        :param figure: GUI inputs a figure to show to
        :param overwrite: Only called by ROI page, when the figure needs to be updated with new ROIs but same frame.
        :return: None. Calls show_rois_func which either plots to output or to figure given
        """
        if experiment_or_dataset == "Experiment":
            self.show_rois_func(self.frame_for_rois, roi_locations=self.rois,
                                roi_size=self.roi_finder.roi_size, figure=figure, overwrite=overwrite)
        elif experiment_or_dataset == "Dataset":
            self.show_rois_func(self.datasets[-1].frame_for_rois,
                                roi_locations=self.datasets[-1].active_rois, figure=figure,
                                roi_size=self.roi_finder.roi_size, roi_offset=self.datasets[-1].roi_offset)

    def finalize_rois(self, name, experiment_settings):
        """
        Finalize ROI settings. Input name and final experiment settings.
        -------------------------
        :param name: Name of experiment
        :param experiment_settings: Final experiment settings
        :return: None, changes class
        """
        self.name = name
        # get directory
        file_dir = '/'.join(self.directory.split(".")[0].split("/")[:-1]) + '/'

        # get date
        try:
            date = self.datasets[0].metadata.pop('date')
            date_split = date.split(" ")[0].split("-")
            # if not split on -, try slashes.
            if len(date_split) == 1:
                date_split = date_split[0].split("/")
            date = "{:04d}-{:02d}-{:02d}".format(int(date_split[2]), int(date_split[1]), int(date_split[0]))
        except:
            date = "XXXX-XX-XX"

        # add date to directory
        file_dir = file_dir + date + "_" + name

        # try to create directory
        directory_try = 0
        while not self.dir_made:
            try:
                mkdir(file_dir)
                self.dir_made = True
            except:
                directory_try += 1
                if directory_try == 1:
                    file_dir += "_%03d" % directory_try
                else:
                    file_dir = file_dir[:-4]
                    file_dir += "_%03d" % directory_try
        # save directory & settings
        self.directory = file_dir
        self.settings = experiment_settings

    def find_rois_dataset(self, settings):
        """
        Find ROIs in latest dataset. Uses the pre-cropping settings.
        -------------------------
        :param settings: Pre cropping settings
        :return: None. Edits dataset class by setting active ROIs in there.
        """
        self.datasets[-1].find_rois(settings)

    def add_to_queue(self, settings):
        """
        Finalize a dataset, adding it to the queue.
        ------------------------
        :param settings: Finalization settings. Given directory to latest dataset
        :return: status: whether or not addition of settings was a success. Mostly edits dataset class.
        """
        status = self.datasets[-1].prepare_run(settings)
        if status is False:
            return status

    def run(self):
        """
        Run experiment. Iterates through each dataset and runs those.
        ---------------
        :return: None. Saves to disk
        """
        start_time = time.time()
        for dataset in self.datasets:
            # iterate through datasets and update progress
            if dataset.type == "TT":
                self.progress_updater.new_dataset(dataset.type, len(dataset.active_rois),
                                                  method=dataset.settings['method'], tt_parts=len(dataset.tt_parts))
            else:
                self.progress_updater.new_dataset(dataset.type, len(dataset.active_rois), method="HSM")
            dataset.run()

            # clear memory
            try:
                dataset.frames.close()
            except:
                pass
            dataset.frames = None

        # save
        time_taken = time.time() - start_time
        self.save(time_taken)

    def save(self, time_taken):
        """
        Saves experiment results
        -------------------
        :return: None. Saves to disk
        """
        # save settings used to txt
        self.progress_updater.message("Starting saving")
        settings = self.settings_to_dict()
        outputting.save_settings(self.directory, settings, time_taken)

        # save overview
        self.progress_updater.message("Saving overview")
        figuring.save_overview(self)  # try except statement within function

        # save individual figures if selected
        if self.settings['All Figures'] is True:
            try:
                self.progress_updater.message("Saving individual figures")
                figuring.individual_figures(self)
            except Exception as e:
                logger.error("Individual figure creation failed")
                logger.info("Info about individual figure creation failed", exc_info=e)

        # convert to matlab coordinates
        self.progress_updater.message("Converting to MATLAB coordinate system")
        tools.convert_to_matlab(self)
        # convert results and metadata to dict
        results = self.rois_to_dict()
        metadata = self.metadata_to_dict()
        # save everything to .mat
        self.progress_updater.message("Saving to .mat")
        outputting.save_to_mat(self.directory, "Results", results)
        outputting.save_to_mat(self.directory, "Metadata", metadata)
        self.progress_updater.message("Done")

    def rois_to_dict(self):
        """
        Converts results to dictionary, also adds roi.x, roi.x, and roi.index
        --------------------
        :return: results_dict: dictionary of results
        """
        result_dict = {}
        for roi in self.rois:
            result_dict["ROI_{}".format(roi.index)] = roi.results
            result_dict["ROI_{}".format(roi.index)]['x'] = roi.x
            result_dict["ROI_{}".format(roi.index)]['y'] = roi.y
            result_dict["ROI_{}".format(roi.index)]['index'] = roi.index

        return result_dict

    def metadata_to_dict(self):
        """
        Converts metadata to dictionary
        ----------------------
        :return: metadata_dict: dictionary of metadata
        """
        metadata_dict = {}
        for dataset in self.datasets:
            metadata_dict['meta_{}'.format(dataset.name_result[4:])] = dataset.metadata

        return metadata_dict

    def settings_to_dict(self):
        """
        Convert settings to dictionary
        -----------------------
        :return: settings_dict: dictionary of settings
        """
        settings_dict = {'Experiment': self.settings, 'ROIs': self.roi_finder.get_settings()}

        for dataset in self.datasets:
            # get settings
            settings_dict[dataset.name] = dataset.settings
            # add type, filename to the front
            settings_dict[dataset.name] = {**{'Type': dataset.type}, **{'filename': dataset.filename},
                                           **settings_dict[dataset.name]}
            settings_dict[dataset.name]['Offset'] = dataset.roi_offset

        return settings_dict
