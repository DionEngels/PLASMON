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
    font_size: size of font of roi number

    Returns
    -------
    None.

    """
    ax.imshow(frame, extent=[0, frame.shape[1], frame.shape[0], 0], aspect='auto')
    if roi_locations is not None and roi_size is not None:
        roi_size_1d = int((roi_size - 1) / 2)
        for roi in roi_locations:

            rect = patches.Rectangle((roi.x-roi_size_1d, roi.y-roi_size_1d), roi_size, roi_size,
                                     linewidth=0.5, edgecolor='r', facecolor='none')
            ax.add_patch(rect)
            if font_size is not None:
                ax.text(roi.x-roi_size_1d, roi.y-roi_size_1d, str(roi.index + 1), color='red', fontsize=font_size)


def find_range(dataset_name, rois):

    max_range = 0

    for roi in rois:
        results = roi.results[dataset_name]['result']
        y_positions = results[:, 1]
        x_positions = results[:, 2]

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


def plot_hsm(ax, result):
    def lorentzian(params, x):
        return params[0] + params[1] / ((x - params[2]) ** 2 + (0.5 * params[3]) ** 2)

    return

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


def make_tt_scatter(ax, result, event_or_not_boolean, dataset):
    ax.clear()
    y_pos = result[invert(event_or_not_boolean), 1]
    x_pos = result[invert(event_or_not_boolean), 2]
    y_pos_event = result[event_or_not_boolean, 1]
    x_pos_event = result[event_or_not_boolean, 2]

    if len(y_pos_event) != 0:
        ax.scatter(x_pos_event, y_pos_event, label='events')
    if len(y_pos) != 0:
        ax.scatter(x_pos, y_pos, label='non-events')

    ax.axis('equal')
    ax.legend()

    if dataset.settings['pixels_or_nm'] == 'nm':
        ax.set_xlabel('x-position (nm)')
        ax.set_ylabel('y-position (nm)')
    else:
        ax.set_xlabel('x-position (pixels)')
        ax.set_ylabel('y-position (pixels)')
    set_range_and_ticks(ax, dataset.figure_range)


def make_tt(ax, time_axis, result, method, roi_index):
    linewidth = 1 / (2 ** (floor(len(time_axis) / 1500) - 1))
    intensities = result[:, 3]
    ax.plot(time_axis, intensities, linewidth=linewidth)

    if "Gaussian" in method:
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Integrated intensity (counts)')
        ax.set_title('Time trace ROI ' + str(roi_index + 1))
    else:
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Summed intensity (counts)')
        ax.set_title('Time trace ROI ' + str(roi_index + 1))


def set_range_and_ticks(ax, max_range):
    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()
    y_center = (y_max + y_min) / 2
    x_center = (x_max + x_min) / 2
    ax.set_xlim(x_center - max_range / 2, x_center + max_range / 2)
    ax.set_ylim(y_center - max_range / 2, y_center + max_range / 2)
    ax.xaxis.set_major_locator(plt.MaxNLocator(N_TICKS))
    ax.yaxis.set_major_locator(plt.MaxNLocator(N_TICKS))


def save_overview(experiment):
    path = experiment.directory + "/Graphs"
    mkdir(path)

    figure_length_base = 6
    tt = []
    hsm = []
    for n_dataset, dataset in enumerate(experiment.datasets):
        if dataset.type == "TT":
            dataset.figure_range = find_range(dataset.name_result, experiment.rois)
            tt.append(n_dataset)
        elif dataset.type == "HSM":
            hsm.append(n_dataset)

    per_roi_length = len(tt) + max(len(hsm) - 1, 0)
    total_length = figure_length_base + 2 * per_roi_length

    fig = plt.figure(constrained_layout=True, figsize=(16, total_length*4), dpi=DPI)
    widths = [1] * 4
    heights = [1] * total_length
    gs = GridSpec(total_length, 4, figure=fig, width_ratios=widths, height_ratios=heights)
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

    roi_list = []
    if len(experiment.rois) > 3:
        for roi in experiment.rois:
            if roi.index % int(len(experiment.rois) / 3) == 0:
                roi_list.append(roi)
        if len(roi_list) < 4:  # sometimes modulo does give four due to small mismatch. This adds a fourth ROI.
            roi_list.append(experiment.rois[-1])
        while len(roi_list) > 4:
            roi_list.pop()
    else:
        for i in range(4):
            value = i % len(experiment.rois)
            roi_list.append(experiment.rois[value])

    for n_roi, roi in enumerate(roi_list):
        row = figure_length_base + int(floor(n_roi / 2)) * per_roi_length
        column = (n_roi % 2) * 2

        ax_frame = fig.add_subplot(gs[row, column])
        plot_rois(ax_frame, roi.get_roi(experiment.frame_for_rois, 7))
        ax_frame.set_xlabel('x (pixels)')
        ax_frame.set_ylabel('y (pixels)')
        ax_frame.set_title('Zoom-in ROI ' + str(roi.index + 1))

        for n_dataset in hsm:
            ax_hsm = fig.add_subplot(gs[row + n_dataset, column + 1])
            plot_hsm(ax_hsm, roi.results[experiment.datasets[n_dataset].name_result])
            ax_hsm.set_xlabel('wavelength (nm)')
            ax_hsm.set_ylabel('intensity (arb. units)')
            ax_hsm.set_title('HSM Result ROI {}'.format(str(roi.index + 1)))

        for n_dataset in tt:
            method = experiment.datasets[n_dataset].settings['method']
            ax_tt_scatter = fig.add_subplot(gs[row + len(hsm) + n_dataset, column])
            make_tt_scatter(ax_tt_scatter, roi.results[experiment.datasets[n_dataset].name_result]['result_post_drift'],
                            roi.results[experiment.datasets[n_dataset].name_result]['event_or_not'],
                            experiment.datasets[n_dataset])
            ax_tt_scatter.set_title('Scatter ROI ' + str(roi.index + 1) + " post drift corr")

            if "Gaussian" in method or "Sum" in method:
                ax_tt = fig.add_subplot(gs[row + len(hsm) + n_dataset, column + 1])
                make_tt(ax_tt, experiment.datasets[n_dataset].time_axis,
                        roi.results[experiment.datasets[n_dataset].name_result]['result'], method, roi.index)

    name = path + "/" + "_Overview.png"
    plt.tight_layout()
    fig.savefig(name, bbox_inches='tight')


def individual_figures(experiment):
    path = experiment.directory + "/Graphs"

    tt = []
    hsm = []
    for n_dataset, dataset in enumerate(experiment.datasets):
        if dataset.type == "TT":
            tt.append(n_dataset)
        elif dataset.type == "HSM":
            hsm.append(n_dataset)

    per_roi_length = len(tt) + max(len(hsm) - 1, 0)

    for roi in experiment.rois:
        fig = plt.figure(figsize=(8, 4*per_roi_length), dpi=DPI)

        ax_frame = fig.add_subplot(per_roi_length, 2 , 1)
        plot_rois(ax_frame, roi.get_roi(experiment.frame_for_rois, 7))
        ax_frame.set_xlabel('x (pixels)')
        ax_frame.set_ylabel('y (pixels)')
        ax_frame.set_title('Zoom-in ROI ' + str(roi.index + 1))

        for n_dataset in hsm:
            ax_hsm = fig.add_subplot(per_roi_length, 2, 1 + n_dataset * 2)
            plot_hsm(ax_hsm, roi.results[experiment.datasets[n_dataset].name_result])
            ax_hsm.set_xlabel('wavelength (nm)')
            ax_hsm.set_ylabel('intensity (arb. units)')
            ax_hsm.set_title('HSM Result ROI {}'.format(str(roi.index + 1)))

        for n_dataset in tt:
            method = experiment.datasets[n_dataset].settings['method']
            ax_tt_scatter = fig.add_subplot(per_roi_length, 2, 3 + n_dataset * 2 + (len(hsm) - 1) * 2)
            make_tt_scatter(ax_tt_scatter, roi.results[experiment.datasets[n_dataset].name_result]['result_post_drift'],
                            roi.results[experiment.datasets[n_dataset].name_result]['event_or_not'],
                            experiment.datasets[n_dataset])
            ax_tt_scatter.set_title('Scatter ROI ' + str(roi.index + 1) + " post drift corr")

            if "Gaussian" in method or "Sum" in method:
                ax_tt = fig.add_subplot(per_roi_length, 2, 4 + n_dataset * 2 + (len(hsm) - 1) * 2)
                make_tt(ax_tt, experiment.datasets[n_dataset].time_axis,
                        roi.results[experiment.datasets[n_dataset].name_result]['result'], method, roi.index)

        name = path + "/" + "ROI" + str(roi.index+1)+".png"
        plt.tight_layout()
        fig.savefig(name, bbox_inches='tight')
