# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
MBx Python Data Analysis

Tools

Some additional tools used by MBx Python

----------------------------

v0.1: Save to CSV & Mat: 31/05/2020
v0.2: also switch array: 04/06/2020
v0.3: cleaned up: 24/07/2020
v0.4: settings and results text output: 25/07/2020

"""

from csv import DictWriter  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from numpy import zeros, savetxt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from datetime import datetime


# translates the raw dictionary keys to user readable input
TRANSLATOR_DICT = {'int_max': 'Maximum Intensity', 'int_min': 'Minimum Intensity', 'sigma_min': "Minimum Sigma",
                   'sigma_max': "Maximum Sigma", 'corr_min': "Minimum Correlation",
                   'pixel_min': "Minimum pixel intensity", 'roi_size': "ROI size", 'filter_size': "Filter size",
                   'roi_side': "Side spacing", 'inter_roi': "ROI spacing"}


def save_to_csv_mat(name, values, path):
    """
    Basic saver to .csv and .mat, only used by metadata

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
        fieldnames = [k[0] for k in values.items()]
        writer = DictWriter(csv_file, fieldnames=fieldnames)

        #  writer.writeheader()
        writer.writerow(values)

        savemat(path + "/" + name + '.mat', values)


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
    with open(path + "/" + name + '.csv', mode='w') as csv_file:
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
            header = "Frame index,ROI index,x position,y position,Intensity Gaussian,Sigma x,Sigma y,Background (" \
                     "fitted),Iterations needed to converge"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)
        else:
            header = "Frame index,ROI index,x position,y position,Intensity Gaussian,Sigma x,Sigma y,Background (" \
                     "estimate),Iterations needed to converge"
            savetxt(path + "/" + name + '.csv', results, delimiter=',', header=header)

    results_dict = {'Localizations': results}
    savemat(path + "/" + name + '.mat', results_dict)


def save_to_csv_mat_roi(name, rois, height, path):
    """
    Saves the ROIs to a .mat and .csv

    Parameters
    ----------
    name : name to save to
    rois : rois to save
    height : height of video to switch y-axis
    path : path to save

    Returns
    -------
    None.

    """
    rois = switch(rois)  # from y,x to x,y
    rois[:, 1] = height - rois[:,1]  # MATLAB has origin in bottom left, Python top left. Switch y-axis to compensate
    header = "x,y"
    savetxt(path + "/" + name + '.csv', rois, delimiter=',', header=header)

    rois_dict = dict(zip(['x', 'y'], rois.T))
    savemat(path + "/" + name + '.mat', rois_dict)


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

        text_file.write("Meaning of variables in Localizations output: \n")
        if method == "Phasor + Intensity":
            text_file.write("Frame index | ROI index | x position | y position | Pixel intensity peak | Background \n")
        elif method == "Phasor":
            text_file.write("Frame index | ROI index | x position | y position \n")
        elif method == "Phasor + Sum":
            text_file.write("Frame index | ROI index | x position | y position | Sum of ROI pixel values \n")
        elif method == "Gaussian - Fit bg":
            text_file.write("Frame index | ROI index | x position | y position | Intensity Gaussian | "
                            "Sigma x | Sigma y | Background (fitted) | Iterations needed to converge \n")
        else:
            text_file.write("Frame index | ROI index | x position | y position | Intensity Gaussian | "
                            "Sigma x | Sigma y | Background (estimate) | Iterations needed to converge \n")


def switch(array):
    """
    Switches a single arrays values

    Parameters
    ----------
    array : array to switch

    Returns
    -------
    new : switched array

    """
    new = zeros(array.shape)
    new[:, 1] = array[:, 0]
    new[:, 0] = array[:, 1]
    return new


def plot_rois(frame, roi_locations, roi_size):
    """

    Parameters
    ----------
    frame : frame to plot
    roi_locations : locations to draw box around
    roi_size : Size of boxes to draw

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(1)
    ax.imshow(frame, extent=[0, frame.shape[1], frame.shape[0], 0], aspect='auto')
    roi_size_1d = int((roi_size - 1) / 2)

    roi_locations = roi_locations - roi_size_1d

    for roi in roi_locations:
        rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                 linewidth=0.5, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    plt.title("ROI locations")
    plt.show()
