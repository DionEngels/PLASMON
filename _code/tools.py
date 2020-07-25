# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
MBx Python Data Analysis

Tools

----------------------------

v1: Save to CSV & Mat: 31/05/2020
v2: also switch array: 04/06/2020
v3: cleaned up: 24/07/2020
v4: settings and results text output: 25/07/2020

"""

from csv import DictWriter  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from numpy import zeros
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from datetime import datetime


def save_to_csv_mat(name, values, directory):
    """

    Parameters
    ----------
    name : name to save to
    values : values to save
    directory : directory to save

    Returns
    -------
    None.

    """
    with open(directory + "/" + name + '.csv', mode='w') as csv_file:
        fieldnames = [k[0] for k in values.items()]
        writer = DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow(values)

        savemat(directory + "/" + name + '.mat', values)


def text_output(settings, method, threshold_method, nm_or_pixels, total_fits, failed_fits, time_taken, directory):
    with open(directory + "/" + "Localizations_info" + ".txt", mode='w') as text_file:
        now = datetime.now()
        text_file.write("Ran on: " + now.strftime('%d/%m/%Y %H:%M:%S') + "\n\n")
        text_file.write("Settings \n------------\n")
        for key, value in settings.items():
            text_file.write(str(key) + ": " + str(value) + "\n")

        text_file.write("\nMethod \n------------\n")
        text_file.write("Method: " + method + "\n")
        text_file.write("Threshold method: " + threshold_method + "\n")
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
