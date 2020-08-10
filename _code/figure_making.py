# -*- coding: utf-8 -*-
"""
Created on Fri August 07 2020

@author: Dion Engels
MBx Python Data Analysis

figure_making

Everything related to making figures

----------------------------

v1.0: split from tools: 07/08/2020
v1.1: individual ROI figures

"""

from os import mkdir
from numpy import asarray, invert, concatenate

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from concurrent.futures import ThreadPoolExecutor as Executor


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


def save_graphs(frames, results, results_drift, roi_locations, method, nm_or_pixels, figures_option, path,
                event_or_not, settings):
    """
    Input results and drift-corrected results, puts out some example graphs for quick checking

    Parameters
    ----------
    frames: frame to be shown
    results: results of fitting
    results_drift: drift-corrected results of fitting
    roi_locations : locations of ROIs
    method: method of fitting
    nm_or_pixels: nm or pixels, for labels
    figures_option: all figures or only overview
    path: path in which to place graphs
    event_or_not: event or not boolean for each roi
    settings: settings used to fit with

    Returns
    -------
    None really. Outputs graphs to disk
    """
    path += "/Graphs"
    mkdir(path)

    if "Gaussian" in method:
        fig = plt.figure(constrained_layout=True, figsize=(16, 40))
        widths = [1] * 4
        heights = [1] * 10
        gs = GridSpec(10, 4, figure=fig, width_ratios=widths, height_ratios=heights)
        ax_overview = fig.add_subplot(gs[0:4, :])
        ax_sigma = fig.add_subplot(gs[4:6, :2])
        ax_int = fig.add_subplot(gs[4:6, 2:])

        overall_intensity = results[:, 4]
        ax_int.hist(overall_intensity, bins=100)
        ax_int.set_xlabel('Intensity (counts)')
        ax_int.set_ylabel('Occurrence')
        ax_int.set_title('Intensity occurrence')

        overall_sigma_x = results[:, 5]
        overall_sigma_y = results[:, 6]
        overall_sigma = concatenate((overall_sigma_x, overall_sigma_y))

        ax_sigma.hist(overall_sigma, bins=50)
        if nm_or_pixels == 'nm':
            ax_sigma.set_xlabel('Sigmas (nm)')
        else:
            ax_sigma.set_xlabel('Sigma (pixels)')
        ax_sigma.set_ylabel('Occurrence')
        ax_sigma.set_title('Sigmas occurrence')

        start_row = 6
    else:
        fig = plt.figure(constrained_layout=True, figsize=(16, 32))
        widths = [1] * 4
        heights = [1] * 8
        gs = GridSpec(8, 4, figure=fig, width_ratios=widths, height_ratios=heights)
        ax_overview = fig.add_subplot(gs[0:4, :])
        start_row = 4

    ax_overview.imshow(frames[0], extent=[0, frames[0].shape[1], frames[0].shape[0], 0], aspect='auto')
    roi_size = settings['roi_size']
    roi_size_1d = int((roi_size - 1) / 2)
    roi_locations_temp = roi_locations - roi_size_1d

    for roi in roi_locations_temp:
        rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                 linewidth=0.5, edgecolor='r', facecolor='none')
        ax_overview.add_patch(rect)

    ax_overview.set_xlabel('x (pixels)')
    ax_overview.set_ylabel('y (pixels)')
    ax_overview.set_title('Full first frame')

    to_save_list = []
    roi_list = []
    if roi_locations.shape[0] > 3:
        for roi_index in range(roi_locations.shape[0]):
            if roi_index % int(roi_locations.shape[0] / 3) == 0:
                roi_list.append(roi_index)
    else:
        for i in range(4):
            value = i % roi_locations.shape[0]
            roi_list.append(value)

    for roi_list_index, roi_index in enumerate(roi_list):
        row = start_row + int(int(roi_list_index / 2)*2)
        column = (roi_index % 2) * 2

        y = int(roi_locations[roi_index, 0])
        x = int(roi_locations[roi_index, 1])

        event_or_not_roi = event_or_not[:, roi_index]
        frame_index = [frame_index for frame_index, true_false in enumerate(event_or_not_roi) if not true_false][0]
        frame = asarray(frames[frame_index])
        my_roi = frame[y - roi_size_1d:y + roi_size_1d + 1, x - roi_size_1d:x + roi_size_1d + 1]

        ax_frame = fig.add_subplot(gs[row, column])
        ax_frame.imshow(my_roi, extent=[0, my_roi.shape[1], my_roi.shape[0], 0], aspect='auto')
        ax_frame.set_xlabel('x (pixels)')
        ax_frame.set_ylabel('y (pixels)')
        ax_frame.set_title('Zoom-in ROI ' + str(roi_index + 1))

        if "Gaussian" in method:
            ax_tt = fig.add_subplot(gs[row, column + 1])
            intensities = results[results[:, 1] == roi_index + 1, 4]
            ax_tt.plot(intensities)
            ax_tt.set_xlabel('Frames')
            ax_tt.set_ylabel('Intensity (counts)')
            ax_tt.set_title('Time trace ROI ' + str(roi_index + 1))

        x_positions = results[results[:, 1] == roi_index + 1, 2]
        y_positions = results[results[:, 1] == roi_index + 1, 3]

        ax_loc = fig.add_subplot(gs[row + 1, column])
        ax_loc.scatter(x_positions[invert(event_or_not_roi)], y_positions[invert(event_or_not_roi)], label='non-events')
        ax_loc.scatter(x_positions[event_or_not_roi], y_positions[event_or_not_roi], label='events')
        if nm_or_pixels == 'nm':
            ax_loc.set_xlabel('x-position (nm)')
            ax_loc.set_ylabel('y-position (nm)')
        else:
            ax_loc.set_xlabel('x-position (pixels)')
            ax_loc.set_ylabel('y-position (pixels)')
        ax_loc.set_title('Scatter ROI ' + str(roi_index + 1))
        ax_loc.axis('equal')
        ax_loc.legend()

        x_positions = results_drift[results_drift[:, 1] == roi_index + 1, 2]
        y_positions = results_drift[results_drift[:, 1] == roi_index + 1, 3]

        ax_loc_drift = fig.add_subplot(gs[row + 1, column + 1])
        ax_loc_drift.scatter(x_positions[invert(event_or_not_roi)], y_positions[invert(event_or_not_roi)], label='non-events')
        ax_loc_drift.scatter(x_positions[event_or_not_roi], y_positions[event_or_not_roi], label='events')
        if nm_or_pixels == 'nm':
            ax_loc_drift.set_xlabel('x-position (nm)')
            ax_loc_drift.set_ylabel('y-position (nm)')
        else:
            ax_loc_drift.set_xlabel('x-position (pixels)')
            ax_loc_drift.set_ylabel('y-position (pixels)')
        ax_loc_drift.set_title('Scatter ROI ' + str(roi_index + 1) + " post drift")
        ax_loc_drift.axis('equal')

    name = path + "/" + "_Overview.png"
    fig.savefig(name, bbox_inches='tight')

    if figures_option == "All":
        for roi_index in range(roi_locations.shape[0]):
            fig = plt.figure(figsize=(8, 8))

            y = int(roi_locations[roi_index, 0])
            x = int(roi_locations[roi_index, 1])
            event_or_not_roi = event_or_not[:, roi_index]
            frame_index = [frame_index for frame_index, true_false in enumerate(event_or_not_roi) if not true_false][0]
            frame = asarray(frames[frame_index])
            my_roi = frame[y - roi_size_1d:y + roi_size_1d + 1, x - roi_size_1d:x + roi_size_1d + 1]

            ax_frame = fig.add_subplot(2, 2, 1)
            ax_frame.imshow(my_roi, extent=[0, my_roi.shape[1], my_roi.shape[0], 0], aspect='auto')
            ax_frame.set_xlabel('x (pixels)')
            ax_frame.set_ylabel('y (pixels)')
            ax_frame.set_title('Zoom-in ROI ' + str(roi_index + 1))

            if "Gaussian" in method:
                ax_tt = fig.add_subplot(2, 2, 2)
                intensities = results[results[:, 1] == roi_index + 1, 4]
                ax_tt.plot(intensities)
                ax_tt.set_xlabel('Frames')
                ax_tt.set_ylabel('Intensity (counts)')
                ax_tt.set_title('Time trace ROI ' + str(roi_index + 1))

            x_positions = results[results[:, 1] == roi_index + 1, 2]
            y_positions = results[results[:, 1] == roi_index + 1, 3]

            ax_loc = fig.add_subplot(2, 2, 3)
            ax_loc.scatter(x_positions[invert(event_or_not_roi)], y_positions[invert(event_or_not_roi)],
                           label='non-events')
            ax_loc.scatter(x_positions[event_or_not_roi], y_positions[event_or_not_roi], label='events')
            if nm_or_pixels == 'nm':
                ax_loc.set_xlabel('x-position (nm)')
                ax_loc.set_ylabel('y-position (nm)')
            else:
                ax_loc.set_xlabel('x-position (pixels)')
                ax_loc.set_ylabel('y-position (pixels)')
            ax_loc.set_title('Scatter ROI ' + str(roi_index + 1))
            ax_loc.axis('equal')
            ax_loc.legend()

            x_positions = results_drift[results_drift[:, 1] == roi_index + 1, 2]
            y_positions = results_drift[results_drift[:, 1] == roi_index + 1, 3]

            ax_loc_drift = fig.add_subplot(2, 2, 4)
            ax_loc_drift.scatter(x_positions[invert(event_or_not_roi)], y_positions[invert(event_or_not_roi)],
                                 label='non-events')
            ax_loc_drift.scatter(x_positions[event_or_not_roi], y_positions[event_or_not_roi], label='events')
            if nm_or_pixels == 'nm':
                ax_loc_drift.set_xlabel('x-position (nm)')
                ax_loc_drift.set_ylabel('y-position (nm)')
            else:
                ax_loc_drift.set_xlabel('x-position (pixels)')
                ax_loc_drift.set_ylabel('y-position (pixels)')
            ax_loc_drift.set_title('Scatter ROI ' + str(roi_index + 1) + " post drift")
            ax_loc_drift.axis('equal')

            name = path + "/" + "ROI" + str(roi_index+1)+".png"
            fig.savefig(name, bbox_inches='tight')
