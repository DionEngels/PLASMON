# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
MBx Python Data Analysis

Tools

----------------------------

v1: Save to CSV & Mat: 31/05/2020
v2: also switch array: 04/06/2020

"""

import csv # to save to csv
import scipy.io as sio #to export for MATLAB
import numpy as np
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
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow(values)

        sio.savemat(directory + "/" + name + '.mat', values)



def switch(array):
    """

    Parameters
    ----------
    array : array to switch

    Returns
    -------
    new : switched array

    """
    new = np.zeros(array.shape)
    new[:, 1] = array[:, 0]
    new[:, 0] = array[:, 1]
    return new

def plot_rois(frame, roi_locations):
    
    plt.imshow(frame, extent=[0,frame.shape[1],frame.shape[0],0], aspect='auto')
    # #takes x,y hence the switched order, and +0.5 for pixel offset
    plt.scatter(roi_locations[:,1] + 0.5, roi_locations[:,0] + 0.5, 
                 s=2, c='red', marker='x', alpha=0.5)
    plt.title("ROI locations")
    plt.show()


def plot_rois_v2(frame, roi_locations, roi_size):
    
    fig, ax = plt.subplots(1)
    ax.imshow(frame, extent=[0,frame.shape[1],frame.shape[0],0], aspect='auto')
    roi_size_1d = int((roi_size-1)/2)
    
    roi_locations = roi_locations - roi_size_1d
    
    for roi in roi_locations:
        rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                 linewidth=0.5, edgecolor='r', facecolor='none')
        ax.add_patch(rect)
        
    plt.title("ROI locations")
    plt.show()