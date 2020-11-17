# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:44:02 2020

----------------------------

@author: Dion Engels
PLASMON Data Analysis

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
v2.0: Program v2: 15/10/2020

 """
# GENERAL IMPORTS
import sys
import logging

from os import getcwd, path, remove  # create directory, check path, and remove
from datetime import datetime  # current time

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib.pylab as plt

# Own code
from src.class_experiment import Experiment
import src.figure_making as figuring

__self_made__ = True

# %% Inputs
ROI_SIZE = 7  # 7 or 9

# %% Initializations

tt_name = "C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2"
hsm_name = "C:/Users/s150127/Downloads/___MBx/datasets/_1nMimager_newGNRs_100mW_HSM/Merged Documents.nd2"

NAME = "test_v2"

fit_options = ["Gaussian - Fit bg", "Gaussian - Estimate bg",
               "Phasor + Intensity", "Phasor + Sum", "Phasor"]

ALL_FIGURES = False
METHOD = "Gaussian - Estimate bg"
REJECTION = True  # True or False
CORRECTION = "SN_objTIRF_PFS_510-800"  # "Matej_670-890"
NM_OR_PIXELS = "nm"
FRAME_BEGIN = "Leave empty for start"  # number or "Leave empty for start"
FRAME_END = 300  # number or "Leave empty for end"
CORR_INT = 500  # "Never" or integer

# %% Proceed question


def proceed_question(title, text):
    """
    Asks user to proceed or not
    :param title: Title
    :param text: Text
    :return: True or False depending on proceed or not
    """
    option1 = "OK"
    option2 = "Cancel"
    answer = input(title + "\n" + text + "\n" + option1 + "/" + option2)
    if answer == option1:
        return True
    else:
        return False

# %% Progress updater non-GUI


class ProgressUpdater:
    """
    Updates the user on progress
    """
    def __init__(self):
        """
        Initializer of ProgressUpdater. Sets a lot of base values
        """
        self.current_type = None
        self.current_dataset = None
        self.dataset_parts = None
        self.dataset_completed = False
        self.method = None
        self.total_datasets = None
        self.current_experiment = None
        self.total_experiments = None
        self.progress = None
        self.total = None
        self.message_string = None

        self.experiment_ind_figures = None
        self.experiment_rois = None

    def start(self, experiments):
        """
        Called when run starts. Sets length of run.
        -------------------------
        :param experiments: experiments to analyze
        :return: None
        """
        self.current_type = None
        self.current_dataset = 0
        self.total_datasets = 0
        self.total_experiments = len(experiments)
        for experiment in experiments:
            self.total_datasets += len(experiment.datasets)

    def new_dataset(self, new_type, n_rois, method=None, tt_parts=None):
        """
        Called when new dataset is starting to be analyzed. Sets new type and updates dataset counter
        -------------------------
        :param new_type: new type of dataset
        :param n_rois: number of active ROIs in dataset
        :param method: method of fitters if used
        :param tt_parts: number of parts the TT has been split in. If not given, parts is set to 1
        :return: Calls update
        """
        self.current_type = new_type
        self.method = method
        self.total = n_rois
        self.progress = 0
        self.current_dataset += 1
        if tt_parts is None:
            self.dataset_parts = 1
        else:
            self.dataset_parts = tt_parts
        self.dataset_completed = False
        self.update(False, True, False)

    def new_experiment(self, exp_index, ind_figures, n_rois):
        """
        Called when new dataset is starting to be analyzed. Sets new experiment index
        --------------------------
        :param exp_index: #experiment
        :param ind_figures: Boolean whether or not individual figures are going to be printed
        :param n_rois: the number of rois to print individual figures for
        :return: Calls update
        """
        self.current_experiment = exp_index
        self.experiment_ind_figures = ind_figures
        self.experiment_rois = n_rois
        self.dataset_completed = False
        self.update(True, False, False)

    def update_progress(self):
        """
        Updates progress within dataset. Adds one to previous progress. Means that one ROI has been completed
        ---------------------------
        :return: Calls update when need be
        """
        self.progress += 1
        # if HSM or Phasor, update every ten
        if (self.method == "HSM" or "Phasor" in self.method) and \
                self.progress % round(self.total * self.dataset_parts / 10, 0) == 0 and \
                self.total * self.dataset_parts > 9:
            self.update(False, False, False)
        # if Gaussian, update every five, since it is slower than Phasor
        elif "Gaussian" in self.method and self.progress % round(self.total * self.dataset_parts / 20, 0) == 0 and \
                self.total * self.dataset_parts > 19:
            self.update(False, False, False)
        # if complete, always update. Also call when only 19 or fewer ROIs
        elif self.total == self.progress or self.total * self.dataset_parts < 20:
            self.update(False, False, False)

        # set to completed when all done
        if self.progress == self.total * self.dataset_parts:
            self.dataset_completed = True

    def message(self, message_string):
        """
        Called when you want to print a message
        :param message_string: Message to print
        :return: Calls update
        """
        self.message_string = message_string
        self.update(False, False, True)

    def update(self, new_experiment, new_dataset, message_bool):
        """
        Update function. Does the actual communication.
        :param new_experiment: Boolean. Whether or not new experiment
        :param new_dataset: Boolean. Whether or not new dataset
        :param message_bool: Boolean. Whether or not message
        :return: prints out info
        """
        if new_experiment:
            print('Starting experiment {}'.format(self.current_experiment))
        elif new_dataset:
            print('Starting dataset {} of {}. Type: {}. Part of Experiment {}'
                  .format(self.current_dataset, self.total_datasets, self.current_type, self.current_experiment))
        elif message_bool:
            print("Experiment {}: ".format(self.current_experiment) + self.message_string)
        else:
            print('{} of {} of current dataset done'.format(self.progress, self.total * self.dataset_parts))

# %% Logging


def logging_setup():
    logger = logging.getLogger('main')
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s: %(filename)s lineno: %(lineno)s:\n%(levelname)s: %(message)s\n')
    # clear log
    if path.isfile(getcwd() + '/Logging/logging.log'):
        try:
            remove(getcwd() + '/Logging/logging.log')
        except:
            pass
    file_handler = logging.FileHandler('Logging/logging.log')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger, formatter

# %% Input error


def input_error(title, text):
    """
    To show input errors to user
    :param title: Title of error
    :param text: Text of error
    :return: Prints out
    """
    print("\033[91m {}\033[00m".format(title + "\n" + text))

# %% General plot


def show_rois(frame, figure=None, roi_locations=None, roi_size=None, roi_offset=None, overwrite=False):
    """
    Shows ROIs within python
    :param frame: frame to make figure of
    :param figure: figure (only used by GUI)
    :param roi_locations: ROI locations within frame
    :param roi_size: ROI size
    :param roi_offset: Offset of ROIs within dataset
    :param overwrite: only used in GUI
    :return:
    """
    if roi_offset is None:
        roi_offset = [0, 0]
    figure, ax = plt.subplots(1)
    figuring.plot_rois(ax, frame, roi_locations, roi_size, roi_offset)
    plt.show()

# %% Run


def run(experiments, progress_updater):
    """
    The actual run command for GUI-less working
    :param experiments: experiments to analyze
    :param progress_updater: progress updater to call
    :return: Prints out and saves to disk
    """
    progress_updater.start(experiments)
    for exp_index, experiment in enumerate(experiments, start=1):
        progress_updater.new_experiment(exp_index, experiment.settings['All Figures'], len(experiment.rois))
        experiment.run()

# %% Main loop cell


if __name__ == '__main__':
    # setup
    logger, formatter = logging_setup()
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.WARNING)
    logger.addHandler(stream_handler)
    logger.info("Started new run\n-----------------------------\n")
    progress_updater = ProgressUpdater()
    experiments = []

    # create experiment
    experiment = Experiment("TT", tt_name, proceed_question, input_error, progress_updater, show_rois)
    experiment.show_rois("Experiment")

    defaults = experiment.roi_finder.get_settings()
    #settings_rois = {'int_max': np.inf, 'int_min': 0,
    #                 'sigma_min': 0, 'sigma_max': int((ROI_SIZE - 1) / 2),
    #                 'corr_min': 0.05, 'roi_size': ROI_SIZE, 'filter_size': 9,
    #                 'roi_side': 11, 'inter_roi': 9}

    # change ROI settings
    #experiment.change_rois(settings_rois)
    #experiment.show_rois("Experiment")

    # finalize experiment
    settings_experiment = {'All Figures': ALL_FIGURES}
    experiment.finalize_rois(NAME, settings_experiment)

    # correlate dataset
    settings_correlation = {'x_min': "Leave empty for start", 'x_max': "Leave empty for end",
                            'y_min': "Leave empty for start", 'y_max': "Leave empty for end"}
    experiment.find_rois_dataset(settings_correlation)
    experiment.show_rois("Dataset")

    # finalize TT dataset
    settings_runtime = {'method': METHOD, 'rejection': REJECTION, '#cores': 1, "pixels_or_nm": NM_OR_PIXELS,
                        'roi_size': ROI_SIZE, 'name': '1nMimager_newGNRs_100mW_TT', "correlation_interval": CORR_INT,
                        'frame_begin': FRAME_BEGIN, 'frame_end': FRAME_END}
    if experiment.add_to_queue(settings_runtime) is False:
        sys.exit("Did not pass check")

    # %% Add HSM
    #experiment.init_new_hsm(hsm_name)

    # correlate HSM ROIs
    #settings_correlation_hsm = {'x_min': "Leave empty for start", 'x_max': "Leave empty for end",
    #                            'y_min': "Leave empty for start", 'y_max': "Leave empty for end"}
    #experiment.find_rois_dataset(settings_correlation_hsm)
    #experiment.show_rois("Dataset")

    # finalize HSM dataset
    #settings_runtime_hsm = {'correction_file': CORRECTION, 'wavelengths': '[510:10:740]',
    #                        'name': '1nMimager_newGNRs_100mW_HSM'}

    #if experiment.add_to_queue(settings_runtime_hsm) is False:
    #    sys.exit("Did not pass check")

    # finalize experiment by adding to experiment list and run
    experiments.append(experiment)
    run(experiments, progress_updater)
