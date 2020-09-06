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

"""

from csv import writer  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from numpy import savetxt
from datetime import datetime

# translates the raw dictionary keys to user readable input
TRANSLATOR_DICT = {'int_max': 'Maximum Intensity', 'int_min': 'Minimum Intensity', 'sigma_min': "Minimum Sigma",
                   'sigma_max': "Maximum Sigma", 'corr_min': "Minimum Correlation",
                   'roi_size': "ROI size", 'filter_size': "Filter size",
                   'roi_side': "Side spacing", 'inter_roi': "ROI spacing",
                   'max_its': "Maximum number of iterations for fitter"}


def save_to_csv_mat_metadata(name, values, path):
    """
    Saver to .csv and .mat for metadata

    Parameters
    ----------
    name : name to save to
    values : values to save
    path : path to save

    Returns
    -------
    None.

    """
    with open(path + "/" + name + '.csv', mode='w') as csv_file:
        writer_item = writer(csv_file)
        writer_item.writerows(values.items())

        values_dict = {name: values}

        savemat(path + "/" + name + '.mat', values_dict, long_field_names=True)


def save_to_csv_mat_results(name, results, method, path):
    """
    Save the results to .csv and .mat

    Parameters
    ----------
    name : name to save to
    results : results to save
    method: method of of fitting used
    path : path to save

    Returns
    -------
    None.

    """
    with open(path + "/" + name + '.csv', mode='w') as _:
        if method == "Phasor + Intensity":
            header = "Frame index,ROI index,x position,y position,Pixel intensity peak,Background"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)
        elif method == "Phasor":
            header = "Frame index,ROI index,x position,y position"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)
        elif method == "Phasor + Sum":
            header = "Frame index,ROI index,x position,y position,Sum of ROI pixel values"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)
        elif method == "Gaussian - Fit bg":
            header = "Frame index,ROI index,x position,y position,Integrated intensity,Sigma x,Sigma y,Background (" \
                     "fitted),Iterations needed to converge"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)
        else:
            header = "Frame index,ROI index,x position,y position,Integrated intensity,Sigma x,Sigma y,Background (" \
                     "estimate),Iterations needed to converge"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)

    results_dict = {name: results}
    savemat(path + "/" + name + '.mat', results_dict)


def save_to_csv_mat_roi(name, rois, path):
    """
    Saves the ROIs to a .mat and .csv

    Parameters
    ----------
    name : name to save to
    rois : rois to save
    path : path to save

    Returns
    -------
    None.

    """
    header = "x,y"
    savetxt(path + "/" + name + '.csv', rois, delimiter=',', header=header)

    rois_dict = {name: rois}
    savemat(path + "/" + name + '.mat', rois_dict)


def save_to_csv_mat_drift(name, drift, path):
    """
    Saves the drift correction to a .mat and .csv

    Parameters
    ----------
    name : name to save to
    drift : drift correction to save
    path : path to save

    Returns
    -------
    None.

    """
    header = "x,y"
    savetxt(path + "/" + name + '.csv', drift, delimiter=',', header=header)

    drift_dict = {name: drift}
    savemat(path + "/" + name + '.mat', drift_dict)


def text_output(settings, method, threshold_method, nm_or_pixels, total_fits, failed_fits, time_taken, path):
    """
    Outputs all settings and results to a txt for reference

    Parameters
    ----------
    settings: settings to find ROIs
    method: method of fitting used
    threshold_method: rejection method used
    nm_or_pixels: nm or pixels output
    total_fits: total fits done
    failed_fits: failed fits
    time_taken: total time taken
    path : path to save

    Returns
    -------
    None.

    """
    with open(path + "/" + "Localizations_info" + ".txt", mode='w') as text_file:
        now = datetime.now()
        text_file.write("Ran on: " + now.strftime('%d/%m/%Y %H:%M:%S') + "\n\n")
        text_file.write("Settings \n------------\n")

        hsm_directory = settings.pop('hsm_directory', "None")
        hsm_correction = settings.pop('hsm_correction', "None")

        for key, value in settings.items():
            text_file.write(str(TRANSLATOR_DICT[key]) + ": " + str(value) + "\n")

        text_file.write("\nMethod \n------------\n")
        text_file.write("Method: " + method + "\n")
        text_file.write("Rejection method: " + threshold_method + "\n")
        text_file.write("Nm or pixels: " + nm_or_pixels + "\n")

        text_file.write("\nResult \n------------\n")
        text_file.write("Total fits: " + str(total_fits) + "\n")
        text_file.write("Failed fits: " + str(failed_fits) + "\n")
        text_file.write("Successful fits: " + str(total_fits - failed_fits) + "\n")
        text_file.write("Time taken: " + str(time_taken) + "\n\n")

        text_file.write("\nHSM \n------------\n")
        text_file.write("Directory: " + str(hsm_directory) + "\n")
        text_file.write("Correction file: " + str(hsm_correction) + "\n")

        text_file.write("Meaning of variables in Localizations output: \n")
        if method == "Phasor + Intensity":
            text_file.write("Frame index | ROI index | x position | y position | Pixel intensity peak | Background \n")
        elif method == "Phasor":
            text_file.write("Frame index | ROI index | x position | y position \n")
        elif method == "Phasor + Sum":
            text_file.write("Frame index | ROI index | x position | y position | Sum of ROI pixel values \n")
        elif method == "Gaussian - Fit bg":
            text_file.write("Frame index | ROI index | x position | y position | Integrated intensity | "
                            "Sigma x | Sigma y | Background (fitted) | Iterations needed to converge \n")
        else:
            text_file.write("Frame index | ROI index | x position | y position | Integrated intensity | "
                            "Sigma x | Sigma y | Background (estimate) | Iterations needed to converge \n")
