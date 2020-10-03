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
v1.2: minor improvement based on Sjoerd's feedback: 27/08/2020
v1.3: feedback of Peter meeting: 06/09/2020

"""

from os import mkdir
from numpy import asarray, invert, concatenate, linspace, nanmin, nanmax
from numpy import max as np_max
from numpy import min as np_min
from math import ceil, floor

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec

__self_made__ = True
DPI = 400
N_TICKS = 4


def plot_rois(ax, frame, roi_locations=None, roi_size=None, font_size=None):
    """

    Parameters
    ----------
    ax: axis to plot to
    frame : frame to plot
    roi_locations : locations to draw box around
    roi_size : Size of boxes to draw

    Returns
    -------
    None.

    """
    ax.imshow(frame, extent=[0, frame.shape[1], frame.shape[0], 0], aspect='auto')
    roi_size_1d = int((roi_size - 1) / 2)

    if roi_locations is not None and roi_size is not None:
        for roi in roi_locations:

            rect = patches.Rectangle((roi.x-roi_size_1d, roi.y-roi_size_1d), roi_size, roi_size,
                                     linewidth=0.5, edgecolor='r', facecolor='none')
            ax.add_patch(rect)
            if font_size is not None:
                ax.text(roi.x, roi.y, str(roi.index + 1), color='red', fontsize=font_size)


def find_range(results_drift, roi_locations):

    max_range = 0

    try:
        numerate_length = roi_locations.shape[0]
    except:
        numerate_length = len(roi_locations)

    for roi_index in range(numerate_length):
        x_positions = results_drift[results_drift[:, 1] == roi_index, 2]
        y_positions = results_drift[results_drift[:, 1] == roi_index, 3]

        if nanmax(x_positions) - nanmin(x_positions) > max_range:
            max_range = nanmax(x_positions) - nanmin(x_positions)
        if nanmax(y_positions) - nanmin(y_positions) > max_range:
            max_range = nanmax(y_positions) - nanmin(y_positions)

    if max_range > 999:
        max_range = ceil(max_range / 100) * 100
    elif max_range > 199:
        max_range = ceil(max_range / 50) * 50
    elif max_range > 99:
        max_range = ceil(max_range / 20) * 20
    else:
        max_range = ceil(max_range / 10) * 10

    return max_range


def save_overview(experiment):
    path = experiment.directory + "/Graphs"
    mkdir(path)

    figure_length = 6
    first_hsm = True
    for dataset in experiment.datasets:
        if dataset.type == "TT":
            figure_length += 1
        elif dataset.type == "HSM":
            if first_hsm:
                pass
            else:
                figure_length += 1

    fig = plt.figure(constrained_layout=True, figsize=(16, figure_length*4), dpi=DPI)
    widths = [1] * 4
    heights = [1] * figure_length
    gs = GridSpec(figure_length, 4, figure=fig, width_ratios=widths, height_ratios=heights)
    ax_overview = fig.add_subplot(gs[0:4, :])
    ax_sigma = fig.add_subplot(gs[4:6, :2])
    ax_int = fig.add_subplot(gs[4:6, 2:])

    plot_rois(ax_overview, experiment.frame_for_rois, roi_locations=experiment.rois, roi_size=7, font_size='large')
    ax_overview.set_xlabel('x (pixels)')
    ax_overview.set_ylabel('y (pixels)')
    ax_overview.set_title('Full first frame')

    ax_int.hist(experiment.roi_finder.int_list, bins=100)
    ax_int.set_xlabel('Integrated intensity (counts)')
    ax_int.set_ylabel('Occurrence')
    ax_int.set_title('Integrated intensity occurrence')

    ax_sigma.hist(experiment.roi_finder.sigma_list, bins=50)
    ax_sigma.set_xlabel('Sigma (pixels)')
    ax_sigma.set_ylabel('Occurrence')
    ax_sigma.set_title('Sigmas occurrence')

    name = path + "/" + "_Overview.png"
    plt.tight_layout()
    fig.savefig(name, bbox_inches='tight')


def save_graphs(frames, results, results_drift, roi_locations, method, nm_or_pixels, figures_option, path,
                event_or_not, settings, time_axis, hsm_result, hsm_raw, hsm_intensity, hsm_wavelengths):
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
    time_axis: time axis of experiment
    hsm_result: HSM results
    hsm_raw: Full HSM fit results
    hsm_intensity: HSM intensities found used to find results
    hsm_wavelengths: wavelengths used at HSM in nm

    Returns
    -------
    None really. Outputs graphs to disk
    """
    def lorentzian(params, x):
        return params[0] + params[1] / ((x - params[2]) ** 2 + (0.5 * params[3]) ** 2)

    if hsm_wavelengths is not None:
        hsm_wavelengths_ev = 1248 / linspace(np_min(hsm_wavelengths), np_max(hsm_wavelengths), num=50)
    else:
        hsm_wavelengths_ev = None

    path += "/Graphs"
    mkdir(path)
    plt.ioff()

    time_axis /= 1000

    linewidth = 1 / (2**(floor(len(time_axis) / 1500) - 1))

    if "Gaussian" in method:
        fig = plt.figure(constrained_layout=True, figsize=(16, 40), dpi=DPI)
        widths = [1] * 4
        heights = [1] * 10
        gs = GridSpec(10, 4, figure=fig, width_ratios=widths, height_ratios=heights)
        ax_overview = fig.add_subplot(gs[0:4, :])
        ax_sigma = fig.add_subplot(gs[4:6, :2])
        ax_int = fig.add_subplot(gs[4:6, 2:])

        overall_intensity = results[:, 4]
        ax_int.hist(overall_intensity, bins=100)
        ax_int.set_xlabel('Integrated intensity (counts)')
        ax_int.set_ylabel('Occurrence')
        ax_int.set_title('Integrated intensity occurrence')

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
        fig = plt.figure(constrained_layout=True, figsize=(16, 32), dpi=DPI)
        widths = [1] * 4
        heights = [1] * 8
        gs = GridSpec(8, 4, figure=fig, width_ratios=widths, height_ratios=heights)
        ax_overview = fig.add_subplot(gs[0:4, :])
        start_row = 4

    ax_overview.imshow(frames[0], extent=[0, frames[0].shape[1], frames[0].shape[0], 0], aspect='auto')
    roi_size = settings['roi_size']
    roi_size_1d = int((roi_size - 1) / 2)
    roi_locations_temp = roi_locations - roi_size_1d

    for roi_index, roi in enumerate(roi_locations_temp):
        rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                 linewidth=0.5, edgecolor='r', facecolor='none')
        ax_overview.add_patch(rect)
        ax_overview.text(roi[1], roi[0], str(roi_index + 1), color='red', fontsize='large')

    ax_overview.set_xlabel('x (pixels)')
    ax_overview.set_ylabel('y (pixels)')
    ax_overview.set_title('Full first frame')

    roi_list = []
    if roi_locations.shape[0] > 3:
        for roi_index in range(roi_locations.shape[0]):
            if roi_index % int(roi_locations.shape[0] / 3) == 0:
                roi_list.append(roi_index)
        if len(roi_list) < 4:  # sometimes modulo does give four due to small mismatch. This adds a fourth ROI.
            roi_list.append(range(roi_locations.shape[0])[-1])
        while len(roi_list) > 4:
            roi_list.pop()
    else:
        for i in range(4):
            value = i % roi_locations.shape[0]
            roi_list.append(value)

    max_range = find_range(results_drift, roi_locations)

    for roi_list_index, roi_index in enumerate(roi_list):
        row = start_row + int(int(roi_list_index / 2)*2)
        column = (roi_list_index % 2) * 2

        y = int(roi_locations[roi_index, 0])
        x = int(roi_locations[roi_index, 1])

        event_or_not_roi = event_or_not[:, roi_index]
        frame_index = [frame_index for frame_index, true_false in enumerate(event_or_not_roi) if not true_false][0]
        frame = asarray(frames[frame_index])
        my_roi = frame[y - roi_size_1d:y + roi_size_1d, x - roi_size_1d:x + roi_size_1d + 1]

        ax_frame = fig.add_subplot(gs[row, column])
        ax_frame.imshow(my_roi, extent=[0, my_roi.shape[1], my_roi.shape[0], 0], aspect='auto')
        ax_frame.set_xlabel('x (pixels)')
        ax_frame.set_ylabel('y (pixels)')
        ax_frame.set_title('Zoom-in ROI ' + str(roi_index + 1))

        if "Gaussian" in method:
            ax_tt = fig.add_subplot(gs[row, column + 1])
            intensities = results[results[:, 1] == roi_index, 4]
            ax_tt.plot(time_axis, intensities, linewidth=linewidth)
            ax_tt.set_xlabel('Time (s)')
            ax_tt.set_ylabel('Integrated intensity (counts)')
            ax_tt.set_title('Time trace ROI ' + str(roi_index + 1))
        elif "Sum" in method:
            ax_tt = fig.add_subplot(gs[row, column + 1])
            intensities = results[results[:, 1] == roi_index, 4]
            ax_tt.plot(time_axis, intensities, linewidth=linewidth)
            ax_tt.set_xlabel('Time (s)')
            ax_tt.set_ylabel('Summed intensity (counts)')
            ax_tt.set_title('Time trace ROI ' + str(roi_index + 1))

        x_positions = results_drift[results_drift[:, 1] == roi_index, 2]
        y_positions = results_drift[results_drift[:, 1] == roi_index, 3]

        ax_loc_drift = fig.add_subplot(gs[row + 1, column])
        ax_loc_drift.scatter(x_positions[invert(event_or_not_roi)], y_positions[invert(event_or_not_roi)],
                             label='non-events')
        ax_loc_drift.scatter(x_positions[event_or_not_roi], y_positions[event_or_not_roi], label='events')
        if nm_or_pixels == 'nm':
            ax_loc_drift.set_xlabel('x-position (nm)')
            ax_loc_drift.set_ylabel('y-position (nm)')
        else:
            ax_loc_drift.set_xlabel('x-position (pixels)')
            ax_loc_drift.set_ylabel('y-position (pixels)')
        ax_loc_drift.set_title('Scatter ROI ' + str(roi_index + 1) + " post drift corr")
        ax_loc_drift.axis('equal')
        ax_loc_drift.legend()

        y_min, y_max = ax_loc_drift.get_ylim()
        x_min, x_max = ax_loc_drift.get_xlim()
        y_center = (y_max + y_min) / 2
        x_center = (x_max + x_min) / 2
        ax_loc_drift.set_xlim(x_center - max_range / 2, x_center + max_range / 2)
        ax_loc_drift.set_ylim(y_center - max_range / 2, y_center + max_range / 2)
        ax_loc_drift.xaxis.set_major_locator(plt.MaxNLocator(N_TICKS))
        ax_loc_drift.yaxis.set_major_locator(plt.MaxNLocator(N_TICKS))

        if hsm_result is None:
            pass
        else:
            ax_hsm = fig.add_subplot(gs[row + 1, column + 1])
            hsm_intensities = hsm_intensity[hsm_intensity[:, 0] == roi_index, 1:]
            ax_hsm.scatter(hsm_wavelengths, hsm_intensities)
            hsm_params = hsm_raw[hsm_raw[:, 0] == roi_index, 1:]
            try:
                hsm_fit = lorentzian(hsm_params.flatten(), hsm_wavelengths_ev)
                ax_hsm.plot(1248 / hsm_wavelengths_ev, hsm_fit, 'r--')
            except:
                pass
            ax_hsm.set_xlabel('wavelength (nm)')
            ax_hsm.set_ylabel('intensity (arb. units)')
            ax_hsm.set_xlim(np_min(hsm_wavelengths) - 30, np_max(hsm_wavelengths) + 30)
            ax_hsm.set_title('HSM Result ROI {}'.format(str(roi_index + 1)))

    name = path + "/" + "_Overview.png"
    plt.tight_layout()
    fig.savefig(name, bbox_inches='tight')

    if figures_option == "All":

        max_range = find_range(results_drift, roi_locations)

        for roi_index in range(roi_locations.shape[0]):
            fig = plt.figure(figsize=(8, 8), dpi=DPI)

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
                intensities = results[results[:, 1] == roi_index, 4]
                ax_tt.plot(time_axis, intensities, linewidth=linewidth)
                ax_tt.set_xlabel('Time (s)')
                ax_tt.set_ylabel('Integrated intensity (counts)')
                ax_tt.set_title('Time trace ROI ' + str(roi_index + 1))
            elif "Sum" in method:
                ax_tt = fig.add_subplot(2, 2, 2)
                intensities = results[results[:, 1] == roi_index, 4]
                ax_tt.plot(time_axis, intensities, linewidth=linewidth)
                ax_tt.set_xlabel('Time (s)')
                ax_tt.set_ylabel('Summed intensity (counts)')
                ax_tt.set_title('Time trace ROI ' + str(roi_index + 1))

            x_positions = results_drift[results_drift[:, 1] == roi_index, 2]
            y_positions = results_drift[results_drift[:, 1] == roi_index, 3]

            ax_loc_drift = fig.add_subplot(2, 2, 3)
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
            ax_loc_drift.legend()

            y_min, y_max = ax_loc_drift.get_ylim()
            x_min, x_max = ax_loc_drift.get_xlim()
            y_center = (y_max + y_min) / 2
            x_center = (x_max + x_min) / 2
            ax_loc_drift.set_xlim(x_center - max_range / 2, x_center + max_range / 2)
            ax_loc_drift.set_ylim(y_center - max_range / 2, y_center + max_range / 2)
            ax_loc_drift.xaxis.set_major_locator(plt.MaxNLocator(N_TICKS))
            ax_loc_drift.yaxis.set_major_locator(plt.MaxNLocator(N_TICKS))

            if hsm_result is None:
                pass
            else:
                ax_hsm = fig.add_subplot(2, 2, 4)
                hsm_intensities = hsm_intensity[hsm_intensity[:, 0] == roi_index, 1:]
                ax_hsm.scatter(hsm_wavelengths, hsm_intensities)
                hsm_params = hsm_raw[hsm_raw[:, 0] == roi_index, 1:]
                try:
                    hsm_fit = lorentzian(hsm_params.flatten(), hsm_wavelengths_ev)
                    ax_hsm.plot(1248 / hsm_wavelengths_ev, hsm_fit, 'r--')
                except:
                    pass
                ax_hsm.set_xlabel('wavelength (nm)')
                ax_hsm.set_ylabel('intensity (arb. units)')
                ax_hsm.set_xlim(np_min(hsm_wavelengths) - 30, np_max(hsm_wavelengths) + 30)
                ax_hsm.set_title('HSM Result ROI {}'.format(str(roi_index + 1)))

            name = path + "/" + "ROI" + str(roi_index+1)+".png"
            plt.tight_layout()
            fig.savefig(name, bbox_inches='tight')
