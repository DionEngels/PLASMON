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

"""

from csv import DictWriter  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from numpy import zeros
import matplotlib.pyplot as plt
import matplotlib.patches as patches


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
