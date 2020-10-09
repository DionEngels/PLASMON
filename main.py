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
import time  # for timekeeping
import sys
import warnings  # for warning diversion

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib.pylab as plt

# v2

from src.class_experiment import Experiment
import src.figure_making as figuring
from src.warnings import InputWarning

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
METHOD = "Gaussian - Fit bg"
THRESHOLD_METHOD = "Loose"  # "Loose", or "None"
CORRECTION = "SN_objTIRF_PFS_510-800"  # "Matej_670-890"
NM_OR_PIXELS = "nm"
FRAME_BEGIN = "Leave empty for start"  # number or "Leave empty for start"
FRAME_END = 10  # number or "Leave empty for end"

# %% Proceed question


def proceed_question(option1, option2, title, text):
    answer = input(title + "\n" + text + "\n" + option1 + "/" + option2)
    if answer == option1:
        return True
    else:
        return False


# %% Progress updater non-GUI


class ProgressUpdater:
    def __init__(self):
        self.current_type = None
        self.current_dataset = None
        self.total_datasets = None
        self.progress = None
        self.total = None
        self.message_bool = False
        self.message_string = None

    def start(self, n_datasets):
        self.current_type = None
        self.current_dataset = 0
        self.total_datasets = n_datasets

    def new_dataset(self, new_type):
        self.message_bool = False
        self.current_type = new_type
        self.current_dataset += 1
        self.update(True)

    def status(self, progress, total):
        self.message_bool = False
        self.progress = progress + 1
        self.total = total
        self.update(False)

    def message(self, message_string):
        self.message_bool = True
        self.message_string = message_string
        self.update(False)

    def update(self, new_dataset):
        if new_dataset:
            print('Starting dataset {} of {}. Type: {}'.format(self.current_dataset, self.total_datasets,
                                                               self.current_type))
        elif self.message_bool:
            print(self.message_string)
        else:
            print('{} of {} of current dataset done'.format(self.progress, self.total))

# %% Divert errors


class DivertError:
    def error(self, *_):
        traceback_details = self.extract_error()
        self.show(True, traceback_details)

    def warning(self, message, category, filename, lineno, file=None, line=None):
        if category == InputWarning:
            message = "Your input is shit"
        elif "Z-levels details missing in metadata" in str(message):
            return  # Only called by metadata loader, which is annoying, thus, not show
        else:
            message = warnings.formatwarning(message, category, filename, lineno)
            message = '\n'.join(message.split('\n')[:-2])
        self.show(False, message)

    @staticmethod
    def extract_error():
        exc_type, exc_value, exc_traceback = sys.exc_info()
        while True:
            try:
                self_made = exc_traceback.tb_next.tb_frame.f_globals['__self_made__']
            except:
                self_made = None
            if self_made is None:
                break
            elif self_made:
                exc_traceback = exc_traceback.tb_next
            else:
                break
        traceback_details = {
            'filename': exc_traceback.tb_frame.f_code.co_filename,
            'lineno': exc_traceback.tb_lineno,
            'name': exc_traceback.tb_frame.f_code.co_name,
            'type': exc_type.__name__,
            'message': exc_value
        }
        return traceback_details

    @staticmethod
    def show(error, traceback_details):
        if error:
            print("\033[91m {}\033[00m".format(traceback_details))
        else:
            print("\033[93m {}\033[00m".format(traceback_details))

# %% General plot


def show_rois(frame, roi_locations=None, roi_size=None):
    fig, ax = plt.subplots(1)
    figuring.plot_rois(ax, frame, roi_locations, roi_size)
    plt.show()

# %% Main loop cell


divertor = DivertError()
warnings.showwarning = divertor.warning

experiment = Experiment("TT", tt_name, proceed_question, ProgressUpdater(), show_rois)

experiment.show_rois("Experiment")

defaults = experiment.roi_finder.get_settings()

settings_rois = {'int_max': np.inf, 'int_min': 0,
                 'sigma_min': 0, 'sigma_max': int((ROI_SIZE - 1) / 2),
                 'corr_min': 0.05, 'roi_size': ROI_SIZE, 'filter_size': 9,
                 'roi_side': 11, 'inter_roi': 9}

experiment.change_rois(settings_rois)

experiment.show_rois("Experiment")

settings_experiment = {'All Figures': ALL_FIGURES}

experiment.finalize_rois(NAME, settings_experiment)

settings_correlation = {'x_min': "Leave empty for start", 'x_max': "Leave empty for end",
                        'y_min': "Leave empty for start", 'y_max': "Leave empty for end"}

experiment.find_rois_dataset(settings_correlation)

experiment.show_rois("Dataset")

settings_runtime = {'method': METHOD, 'rejection': THRESHOLD_METHOD, '#cores': 1, "pixels_or_nm": NM_OR_PIXELS,
                    'roi_size': ROI_SIZE,
                    'frame_begin': FRAME_BEGIN, 'frame_end': FRAME_END}

status = experiment.add_to_queue(settings_runtime)
if status is False:
    sys.exit("Did not pass check")

# %% Add HSM

experiment.init_new_hsm(hsm_name)

settings_correlation_hsm = {'x_min': "Leave empty for start", 'x_max': "Leave empty for end",
                            'y_min': "Leave empty for start", 'y_max': "Leave empty for end"}

experiment.find_rois_dataset(settings_correlation_hsm)

settings_runtime_hsm = {'correction_file': CORRECTION, 'wavelengths': '[510:10:740]'}

status = experiment.add_to_queue(settings_runtime_hsm)
if status is False:
    sys.exit("Did not pass check")

experiment.run()
