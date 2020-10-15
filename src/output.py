# -*- coding: utf-8 -*-
"""
Created on Fri August 07 2020

@author: Dion Engels
MBx Python Data Analysis

Output

Everything related to saving the results

----------------------------

v1.0: split from tools: 07/08/2020
v1.1: Integrated intensity instead of peak Gaussian intensity: 27/08/2020
v2.0: Part of program v2: 15/10/2020

"""
from scipy.io import savemat  # to export for MATLAB
from datetime import datetime

__self_made__ = True

# translates the raw dictionary keys to user readable input
TRANSLATOR_DICT = {'int_max': 'Maximum Intensity', 'int_min': 'Minimum Intensity', 'sigma_min': "Minimum Sigma",
                   'sigma_max': "Maximum Sigma", 'corr_min': "Minimum Correlation",
                   'roi_size': "ROI size", 'filter_size': "Filter size",
                   'roi_side': "Side spacing", 'inter_roi': "ROI spacing",
                   'max_its': "Maximum number of iterations for fitter", 'All Figures': "All Figures",
                   'rejection': "Rejection method", '#cores': "Number of cores used",
                   'pixels_or_nm': "Pixels or nanometer output",
                   'frame_begin': "First frame fitted", 'frame_end': 'Last frame fitted', 'Type': "Type of dataset",
                   'Offset': "Offset compared to ROI finding frame", 'correction_file': "HSM spectral correction",
                   'wavelengths': "HSM wavelengths", 'filename': "Filename"}


def save_to_mat(directory, name, to_save):
    """
    Saves to matlab .mat file.
    -------------------------
    :param directory: directory to save to
    :param name: name to save with
    :param to_save: data to save
    :return: None. Saves to disk.
    """
    savemat(directory + "/" + name + '.mat', to_save, long_field_names=True)


def save_settings(directory, settings):
    """
    Save settings to .txt file
    -------------------------------
    :param directory: directory to save to
    :param settings: settings dictionary to save
    :return: None. Saves to disk.
    """
    with open(directory + "/" + "Settings" + ".txt", mode='w') as text_file:
        # get current time and save that
        now = datetime.now()
        text_file.write("Ran on: " + now.strftime('%d/%m/%Y %H:%M:%S') + "\n\n")

        # write experiment settings
        settings_experiment = settings.pop('Experiment')
        text_file.write("Experiment \n------------\n")
        for key, value in settings_experiment.items():
            text_file.write(str(TRANSLATOR_DICT[key]) + ": " + str(value) + "\n")

        # write ROI finder settings
        settings_roi_finder = settings.pop('ROIs')
        text_file.write("\n ROI Finder \n------------\n")
        for key, value in settings_roi_finder.items():
            if key == 'processed_frame':
                pass
            else:
                text_file.write(str(TRANSLATOR_DICT[key]) + ": " + str(value) + "\n")

        # write settings per dataset
        for dataset_name, dataset_settings in settings.items():
            text_file.write("\n {} \n------------\n".format(dataset_name))
            for key, value in dataset_settings.items():
                if key == "method":
                    text_file.write("Meaning of variables in Localizations output: \n")
                    if value == "Phasor + Intensity":
                        text_file.write(
                            "Frame index | ROI index | x position | y position | Pixel intensity peak | Background \n")
                    elif value == "Phasor":
                        text_file.write("Frame index | ROI index | x position | y position \n")
                    elif value == "Phasor + Sum":
                        text_file.write(
                            "Frame index | ROI index | x position | y position | Sum of ROI pixel values \n")
                    elif value == "Gaussian - Fit bg":
                        text_file.write("Frame index | ROI index | x position | y position | Integrated intensity | "
                                        "Sigma x | Sigma y | Background (fitted) | Iterations needed to converge \n")
                    else:
                        text_file.write("Frame index | ROI index | x position | y position | Integrated intensity | "
                                        "Sigma x | Sigma y | Background (estimate) | Iterations needed to converge \n")
                else:
                    text_file.write(str(TRANSLATOR_DICT[key]) + ": " + str(value) + "\n")
