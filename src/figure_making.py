# -*- coding: utf-8 -*-
"""
Created on Fri August 07 2020

@author: Dion Engels
PLASMON Data Analysis

figure_making

Everything related to making figures

----------------------------

v1.0: split from tools: 07/08/2020
v1.1: individual ROI figures
v1.2: minor improvement based on Sjoerd's feedback: 27/08/2020
v1.3: feedback of Peter meeting: 06/09/2020
v2.0: Figure making for v2.0 of GUI: 15/10/2020

"""

from os import mkdir
from numpy import invert, linspace, nanmin, nanmax
from numpy import max as np_max
from numpy import min as np_min
from math import ceil, floor

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import logging  # for logging warnings
logger = logging.getLogger('main')

__self_made__ = True
DPI = 400
N_TICKS = 4


def plot_rois(ax, frame, roi_locations=None, roi_size=None, roi_offset=None, font_size=None, overwrite=False):
    """
    Plots a microscope with rois if deisred
    ----------
    :param ax: axis to plot to
    :param frame : frame to plot
    :param roi_locations : locations to draw box around
    :param roi_size : Size of boxes to draw
    :param roi_offset: offset of the dataset compared to experiment ROIs
    :param font_size: size of font of roi number
    :param overwrite: Only called by ROI page, when the figure needs to be updated with new ROIs but same frame.
    :returns None. Edits ax object
    """
    if roi_offset is None:
        roi_offset = [0, 0]  # if no offset given, set to 0

    if not overwrite:
        # show frame
        ax.imshow(frame, extent=[0, frame.shape[1], frame.shape[0], 0], aspect='auto')

    # add ROIs if not None
    if roi_locations is not None and roi_size is not None:
        roi_size_1d = int((roi_size - 1) / 2)
        for roi in roi_locations:
            # red square for each ROI
            rect = patches.Rectangle((roi.x - roi_size_1d + roi_offset[1], roi.y + roi_offset[0] - roi_size_1d),
                                     roi_size, roi_size, linewidth=0.5, edgecolor='r', facecolor='none')
            ax.add_patch(rect)
            if font_size is not None:
                # ROI number
                ax.text(roi.x + roi_offset[1] - roi_size_1d, roi.y + roi_offset[0] - roi_size_1d, str(roi.index + 1),
                        color='red', fontsize=font_size)


def find_range(dataset_name, rois):
    """
    Finds the range of all scatter plots
    ----------------------
    :param dataset_name: Name of the dataset to get range for
    :param rois: All the ROIs
    :return: max_range: the highest range, will afterwards be used for all datasets
    """

    max_range = 0

    for roi in rois:
        # get results per ROI and find max_range of ROI. If larger than before, save.
        results = roi.results[dataset_name]['result']
        y_positions = results[:, 1]
        x_positions = results[:, 2]

        if nanmax(x_positions) - nanmin(x_positions) > max_range:
            max_range = nanmax(x_positions) - nanmin(x_positions)
        if nanmax(y_positions) - nanmin(y_positions) > max_range:
            max_range = nanmax(y_positions) - nanmin(y_positions)

    # settle on nice round number
    if max_range > 999:
        max_range = ceil(max_range / 100) * 100
    elif max_range > 199:
        max_range = ceil(max_range / 50) * 50
    elif max_range > 99:
        max_range = ceil(max_range / 20) * 20
    else:
        max_range = ceil(max_range / 10) * 10

    return max_range


def plot_hsm(ax, result, wavelengths):
    """
    Plotter of a HSM result
    ----------------------
    :param ax: ax object to edit
    :param result: result to plot
    :param wavelengths: wavelengths used to get result. Will become x-axis
    :return: None. Edits ax object
    """
    def lorentzian(params, x):
        """
        Lorentzian formula. Taken from SPectrA
        ----------------
        :param params: Parameters of lorentzian. Need to be four.
        :param x: x-axis. Wavelengths
        :return: array of values for current parameters and wavelengths
        """
        return params[0] + params[1] / ((x - params[2]) ** 2 + (0.5 * params[3]) ** 2)

    # scatter plot of actual results
    ax.scatter(wavelengths, result['intensity'])
    try:
        # try to fit and plot fit over
        wavelengths_ev = 1240 / linspace(np_min(wavelengths), np_max(wavelengths), num=50)
        hsm_fit = lorentzian(result['fit_parameters'], wavelengths_ev)
        ax.plot(1240 / wavelengths_ev, hsm_fit, 'r--')
    except:
        pass
    # set axis to same range every time
    ax.set_xlim(np_min(wavelengths) - 30, np_max(wavelengths) + 30)
    try:
        # add fit info if possible
        text = '\n'.join((r'SPR (nm)=%.1f' % (result['result'][0],),
                          r'linewidth (meV)=%.1f' % (result['result'][1],),
                          r'r^2=%.1f' % (result['result'][2],)))
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.05, 0.95, text, transform=ax.transAxes, verticalalignment='top', bbox=props)
    except:
        pass


def make_tt_scatter(ax, result, event_or_not_boolean, dataset):
    """
    Scatter plot of TT results
    -------------------------
    :param ax: ax object to edit
    :param result: results to plot
    :param event_or_not_boolean: boolean that tells whether or not each time step is event or not
    :param dataset: dataset object from which current results are obtained
    :return: None. Edits ax object
    """
    ax.clear()  # clear to be sure
    # plot events and non-events
    y_pos = result[invert(event_or_not_boolean), 1]
    x_pos = result[invert(event_or_not_boolean), 2]
    y_pos_event = result[event_or_not_boolean, 1]
    x_pos_event = result[event_or_not_boolean, 2]

    # scatter plot and set colours
    if len(y_pos_event) != 0:
        ax.scatter(x_pos_event, y_pos_event, color='red')
    if len(y_pos) != 0:
        ax.scatter(x_pos, y_pos, color='#1f77b4')

    # force same legend
    legend_events = Line2D([0], [0], marker='o', color='w', label='events',
                           markerfacecolor='red', markersize=10)
    legend_non_events = Line2D([0], [0], marker='o', color='w', label='non-events',
                               markerfacecolor='#1f77b4', markersize=10)

    # set square axis and put in legend
    ax.axis('equal')
    ax.legend(handles=[legend_events, legend_non_events])

    # set labels, range, and ticks
    if dataset.settings['pixels_or_nm'] == 'nm':
        ax.set_xlabel('x-position (nm)')
        ax.set_ylabel('y-position (nm)')
    else:
        ax.set_xlabel('x-position (pixels)')
        ax.set_ylabel('y-position (pixels)')
    set_range_and_ticks(ax, dataset.figure_range)


def set_range_and_ticks(ax, max_range):
    """
    Sets range and ticks for each scatter plot
    --------------------------------------
    :param ax: ax object to edit
    :param max_range: max range to put scatter too
    :return: None. Edits ax object
    """
    # get x and y range
    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()
    # set new x and y ranges
    y_center = (y_max + y_min) / 2
    x_center = (x_max + x_min) / 2
    ax.set_xlim(x_center - max_range / 2, x_center + max_range / 2)
    ax.set_ylim(y_center - max_range / 2, y_center + max_range / 2)
    # set ticks
    ax.xaxis.set_major_locator(plt.MaxNLocator(N_TICKS))
    ax.yaxis.set_major_locator(plt.MaxNLocator(N_TICKS))


def make_tt(ax, time_axis, time_axis_dim, result, method):
    """
    Makes line plot of TT
    --------------------
    :param ax: ax object to edit
    :param time_axis: time axis, the x-axis of the plot
    :param time_axis_dim: the dimension of the x-axis of the plot
    :param result: results to plot
    :param method: method used to get results. Impacts labels
    :return: None. Edits ax object
    """
    linewidth = 1 / (2 ** (floor(len(time_axis) / 1500) - 1))  # variable linewidth depending on length result
    # plot resulting intensities
    intensities = result[:, 3]
    ax.plot(time_axis, intensities, linewidth=linewidth)

    # put in label depending on method
    if "Gaussian" in method:
        ax.set_ylabel('Integrated intensity (counts)')
    else:
        ax.set_ylabel('Summed intensity (counts)')
    if time_axis_dim == 't':
        ax.set_xlabel('Time (s)')
    else:
        ax.set_xlabel('Frames (-)')


def save_overview(experiment):
    """
    Makes overview figure of experiment
    ---------------------------------
    :param experiment: Experiment to make overview figure of
    :return: None. Saves figure to disk
    """
    # force agg backend. Otherwise breaks due to threading
    mpl.use('agg', force=True)
    from matplotlib import pyplot as plt
    # set title size to small for all figures
    mpl.rcParams['axes.titlesize'] = 'small'

    # make graphs directory
    path = experiment.directory + "/Graphs"
    mkdir(path)
    name = path + "/" + "_Overview.png"

    # determine number of datasets in experiment and required figure size
    figure_length_base = 6
    tt = []
    hsm = []
    for n_dataset, dataset in enumerate(experiment.datasets):
        if dataset.type == "TT":
            dataset.figure_range = find_range(dataset.name_result, dataset.active_rois)
            tt.append(n_dataset)
        elif dataset.type == "HSM":
            hsm.append(n_dataset)

    per_roi_length = len(tt) + max(len(hsm), 1)
    total_length = figure_length_base + 2 * per_roi_length

    # make figure and grid in figure
    fig = plt.figure(constrained_layout=True, figsize=(16, total_length * 4), dpi=DPI)
    widths = [1] * 4
    heights = [1] * total_length
    gs = GridSpec(total_length, 4, figure=fig, width_ratios=widths, height_ratios=heights)
    try:
        # create base figures, overview, sigma histogram, and intensity histogram
        ax_overview = fig.add_subplot(gs[0:4, :])
        ax_sigma = fig.add_subplot(gs[4:6, :2])
        ax_int = fig.add_subplot(gs[4:6, 2:])

        # plot ROIs to overview
        plot_rois(ax_overview, experiment.frame_for_rois, roi_locations=experiment.rois, roi_size=7, font_size='large')
        ax_overview.set_xlabel('x (pixels)')
        ax_overview.set_ylabel('y (pixels)')
        ax_overview.set_title('Full first frame')

        # Histogram of intensities
        ax_int.hist(experiment.roi_finder.int_list, bins=100)
        ax_int.set_xlabel('Integrated intensity (counts)')
        ax_int.set_ylabel('Occurrence')
        ax_int.set_title('Integrated intensity occurrence')

        # Histogram of sigma
        ax_sigma.hist(experiment.roi_finder.sigma_list, bins=50)
        ax_sigma.set_xlabel('Sigma (pixels)')
        ax_sigma.set_ylabel('Occurrence')
        ax_sigma.set_title('Sigmas occurrence')

        # find the ROIs with all datasets if possible
        full_rois = []
        for roi in experiment.rois:
            if len(roi.results) == len(experiment.datasets):
                full_rois.append(roi)

        # if not enough full rois, just randomly sample all ROIs
        if len(full_rois) < 4:
            full_rois = experiment.rois
            logger.warning("Too few ROIs that are active in all datasets. Overview figure might have some empty slots")

        # sample ROIs to put in overview
        roi_list = []
        if len(full_rois) > 3:  # if 4 or more, just take those using modulo
            for roi in full_rois:
                if roi.index % int(len(full_rois) / 3) == 0:
                    roi_list.append(roi)
            if len(roi_list) < 4:  # sometimes modulo does give four due to small mismatch. This adds a fourth ROI.
                roi_list.append(full_rois[-1])
            while len(roi_list) > 4:
                roi_list.pop()
        else:
            for i in range(4):  # otherwise, double sample. Always get four
                value = i % len(full_rois)
                roi_list.append(full_rois[value])

        # for each ROI, plot all datasets
        for n_roi, roi in enumerate(roi_list):
            row = figure_length_base + int(floor(n_roi / 2)) * per_roi_length
            column = (n_roi % 2) * 2

            # first zoomed in frame
            ax_frame = fig.add_subplot(gs[row, column])
            plot_rois(ax_frame, roi.get_roi(experiment.frame_for_rois, 7, [0, 0]))
            ax_frame.set_xlabel('x (pixels)')
            ax_frame.set_ylabel('y (pixels)')
            ax_frame.set_title('Zoom-in ROI {}'.format(roi.index + 1))

            for index_dataset, n_dataset in enumerate(hsm):
                try:
                    # then each HSM
                    ax_hsm = fig.add_subplot(gs[row + index_dataset, column + 1])
                    plot_hsm(ax_hsm, roi.results[experiment.datasets[n_dataset].name_result],
                             experiment.datasets[n_dataset].wavelengths)
                    ax_hsm.set_xlabel('wavelength (nm)')
                    ax_hsm.set_ylabel('intensity (arb. units)')
                    ax_hsm.set_title('HSM {}\nROI {}'.format(experiment.datasets[n_dataset].name, roi.index + 1))
                except:
                    pass  # if this ROI does not have results for that dataset, skip

            for index_dataset, n_dataset in enumerate(tt):
                try:
                    # then each TT, first scatter
                    method = experiment.datasets[n_dataset].settings['method']
                    ax_tt_scatter = fig.add_subplot(gs[row + max(len(hsm), 1) + index_dataset, column])
                    make_tt_scatter(ax_tt_scatter,
                                    roi.results[experiment.datasets[n_dataset].name_result]['result_post_drift'],
                                    roi.results[experiment.datasets[n_dataset].name_result]['event_or_not'],
                                    experiment.datasets[n_dataset])
                    ax_tt_scatter.set_title('Scatter {}\nw/ drift corr ROI {}'.format(experiment.datasets[n_dataset].name,
                                                                                     roi.index + 1))
                    if "Gaussian" in method or "Sum" in method:
                        # and if possible, time trace
                        ax_tt = fig.add_subplot(gs[row + max(len(hsm), 1) + index_dataset, column + 1])
                        make_tt(ax_tt, experiment.datasets[n_dataset].time_axis,
                                experiment.datasets[n_dataset].time_axis_dim,
                                roi.results[experiment.datasets[n_dataset].name_result]['result'], method)
                        ax_tt.set_title('TT {}\nROI {}'.format(experiment.datasets[n_dataset].name, roi.index + 1))
                except:
                    pass  # if this ROI does not have results for that dataset, skip

        # save
        plt.tight_layout()
        fig.savefig(name, bbox_inches='tight')
        fig.clear()
        plt.close(fig)
    except Exception as e:
        # throw warning
        logger.error("Overview figure creation failed")
        logger.info("Info about overview figure creation failed", exc_info=e)
        # in case of crash, just save what you got
        plt.tight_layout()
        fig.savefig(name, bbox_inches='tight')
        fig.clear()
        plt.close(fig)


def individual_figures(experiment):
    """
    Makes individual figure of experiment
    ---------------------------------
    :param experiment: Experiment to make overview figure of
    :return: None. Saves figures to disk
    """
    # force agg backend. Otherwise breaks due to threading
    mpl.use('agg', force=True)
    from matplotlib import pyplot as plt
    # set title size to small for all figures
    mpl.rcParams['axes.titlesize'] = 'small'

    # set path, already created by overview
    path = experiment.directory + "/Graphs"

    # determine number of datasets in experiment and required figure size
    tt = []
    hsm = []
    for n_dataset, dataset in enumerate(experiment.datasets):
        if dataset.type == "TT":
            tt.append(n_dataset)
        elif dataset.type == "HSM":
            hsm.append(n_dataset)

    per_roi_length = len(tt) + max(len(hsm), 1)

    # for each ROI, plot all datasets
    for roi in experiment.rois:
        fig = plt.figure(figsize=(8, 4*per_roi_length), dpi=DPI)

        # first zoomed in frame
        ax_frame = fig.add_subplot(per_roi_length, 2, 1)
        plot_rois(ax_frame, roi.get_roi(experiment.frame_for_rois, 7, [0, 0]))
        ax_frame.set_xlabel('x (pixels)')
        ax_frame.set_ylabel('y (pixels)')
        ax_frame.set_title('Zoom-in ROI {}'.format(roi.index + 1))

        for index_dataset, n_dataset in enumerate(hsm):
            try:
                # then each HSM
                ax_hsm = fig.add_subplot(per_roi_length, 2, 2 + index_dataset * 2)
                plot_hsm(ax_hsm, roi.results[experiment.datasets[n_dataset].name_result],
                         experiment.datasets[n_dataset].wavelengths)
                ax_hsm.set_xlabel('wavelength (nm)')
                ax_hsm.set_ylabel('intensity (arb. units)')
                ax_hsm.set_title('HSM {}\nROI {}'.format(experiment.datasets[n_dataset].name, str(roi.index + 1)))
            except:
                pass  # if this ROI does not have results for that dataset, skip

        for index_dataset, n_dataset in enumerate(tt):
            try:
                # then each time trace scatter
                method = experiment.datasets[n_dataset].settings['method']
                ax_tt_scatter = fig.add_subplot(per_roi_length, 2, 1 + index_dataset * 2 + max(len(hsm), 1) * 2)
                make_tt_scatter(ax_tt_scatter,
                                roi.results[experiment.datasets[n_dataset].name_result]['result_post_drift'],
                                roi.results[experiment.datasets[n_dataset].name_result]['event_or_not'],
                                experiment.datasets[n_dataset])
                ax_tt_scatter.set_title('Scatter {}\nw/ drift corr ROI {}'.format(experiment.datasets[n_dataset].name,
                                                                                 roi.index + 1))

                if "Gaussian" in method or "Sum" in method:
                    # and if possible, time trace
                    ax_tt = fig.add_subplot(per_roi_length, 2, 2 + index_dataset * 2 + max(len(hsm), 1) * 2)
                    make_tt(ax_tt, experiment.datasets[n_dataset].time_axis,
                            experiment.datasets[n_dataset].time_axis_dim,
                            roi.results[experiment.datasets[n_dataset].name_result]['result'], method)
                    ax_tt.set_title('TT {}\nROI {}'.format(experiment.datasets[n_dataset].name, roi.index + 1))
            except:
                pass  # if this ROI does not have results for that dataset, skip

        # save
        name = path + "/" + "ROI" + str(roi.index+1)+".png"
        plt.tight_layout()
        fig.savefig(name, bbox_inches='tight')
        fig.clear()
        plt.close(fig)
