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
v0.5: own ND2Reader class to prevent warnings: 29/07/2020
v0.6: save drift and save figures: 03/08/2020
v0.6.1: better MATLAB ROI and Drift output

"""

from csv import DictWriter  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from numpy import zeros, savetxt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from datetime import datetime
from pims_nd2 import ND2_Reader
from nd2reader import ND2Reader
from nd2reader.parser import Parser
from os import mkdir

# translates the raw dictionary keys to user readable input
TRANSLATOR_DICT = {'int_max': 'Maximum Intensity', 'int_min': 'Minimum Intensity', 'sigma_min': "Minimum Sigma",
                   'sigma_max': "Maximum Sigma", 'corr_min': "Minimum Correlation",
                   'roi_size': "ROI size", 'filter_size': "Filter size",
                   'roi_side': "Side spacing", 'inter_roi': "ROI spacing"}

# %% Class to load in ND2s


class ND2ReaderForMetadata(ND2Reader):

    def __init__(self, filename):
        super(ND2Reader, self).__init__()
        self.filename = filename

        # first use the parser to parse the file
        self._fh = open(filename, "rb")
        self._parser = Parser(self._fh)

        # Setup metadata
        self.metadata = self._parser.metadata

        # Set data type
        self._dtype = self._parser.get_dtype_from_metadata()

        # Setup the axes
        self._setup_axes()

        # Other properties
        self._timesteps = None

    def get_metadata(self):
        metadata_dict = self.metadata

        metadata_dict.pop('rois')
        metadata_dict.pop('z_levels')
        metadata_dict.pop('frames')
        metadata_dict.pop('date')

        metadata_dict['pfs_status'] = self._parser._raw_metadata.pfs_status
        metadata_dict['pfs_offset'] = self._parser._raw_metadata.pfs_offset

        metadata_dict['timesteps'] = self.timesteps
        metadata_dict['frame_rate'] = self.frame_rate

        return metadata_dict


class ND2ReaderSelf(ND2_Reader):
    """
    Small class to read in ND2 using a prebuild ND2 Reader. Slightly edited to prevent it giving a warning
    """
    def __init__(self, filename, series=0, channel=0):
        self._clear_axes()
        self._get_frame_dict = dict()
        super().__init__(filename, series=series, channel=channel)

    def get_metadata(self):
        metadata_dict = self.metadata
        metadata_dict_filtered = {k: v for k, v in metadata_dict.items() if v is not None}
        del metadata_dict_filtered['time_start']
        del metadata_dict_filtered['time_start_utc']

        nd2_part_2 = ND2ReaderForMetadata(self.filename)
        metadata_dict_part2 = nd2_part_2.get_metadata()
        total_metadata = {**metadata_dict_filtered, **metadata_dict_part2}

        nd2_part_2.close()

        return total_metadata


# %% Saving stuff

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
        fieldnames = [k[0] for k in values.items()]
        writer = DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow(values)

        values_dict = {name: values}

        savemat(path + "/" + name + '.mat', values_dict)


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

    results_dict = {name: results}
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
    rois[:, 1] = height - rois[:, 1]  # MATLAB has origin in bottom left, Python top left. Switch y-axis to compensate
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


def save_graphs(results, results_drift, roi_locations, method, nm_or_pixels, figures_option, path):
    """
    Input results and drift-corrected results, puts out some example graphs for quick checking

    Parameters
    ----------
    results: results of fitting
    results_drift: drift-corrected results of fitting
    roi_locations : locations of ROIs
    method: method of fitting
    nm_or_pixels: nm or pixels, for labels
    figures_option: amount of figures to be plotted
    path: path in which to place graphs

    Returns
    -------
    None really. Outputs graphs to disk
    """
    path += "/Graphs"
    mkdir(path)

    fig = plt.figure()
    ax = plt.subplot(111)

    if "Gaussian" in method:
        overall_intensity = results[:, 4]
        ax.hist(overall_intensity, bins=100)
        ax.set_xlabel('Intensity (counts)')
        ax.set_ylabel('Occurrence')
        ax.set_title('Intensity occurrence')
        name = path + "/" + "_Histogram_intensity.png"
        fig.savefig(name, bbox_inches='tight')

        plt.cla()
        overall_sigma_x = results[:, 5]
        overall_sigma_y = results[:, 6]
        ax.hist(overall_sigma_x, bins=50)
        ax.hist(overall_sigma_y, bins=50)
        ax.set_xlabel('Sigmas (counts)')
        ax.set_ylabel('Occurrence')
        ax.set_title('Sigmas occurrence')
        name = path + "/" + "_Histogram_sigma.png"
        fig.savefig(name, bbox_inches='tight')

    for roi_index, roi in enumerate(roi_locations):
        if figures_option == "Few" and roi_index % int(roi_locations.shape[0]/5) != 0:
            continue

        if "Gaussian" in method:
            intensities = results[results[:, 1] == roi_index+1, 4]
            plt.cla()
            ax.plot(intensities)
            ax.set_xlabel('Frames')
            ax.set_ylabel('Intensity (counts)')
            ax.set_title('Time trace ROI ' + str(roi_index+1))
            name = path + "/" + str(roi_index+1) + "_time_trace.png"
            fig.savefig(name, bbox_inches='tight')

        x_positions = results[results[:, 1] == roi_index+1, 2]
        y_positions = results[results[:, 1] == roi_index + 1, 3]

        plt.cla()
        ax.scatter(x_positions, y_positions)
        if nm_or_pixels == 'nm':
            ax.set_xlabel('x-position (nm)')
            ax.set_ylabel('y-position (nm)')
        else:
            ax.set_xlabel('x-position (pixels)')
            ax.set_ylabel('y-position (pixels)')
        ax.set_title('Scatter ROI ' + str(roi_index + 1))
        name = path + "/" + str(roi_index + 1) + "_scatter.png"
        fig.savefig(name, bbox_inches='tight')

        x_positions = results_drift[results_drift[:, 1] == roi_index + 1, 2]
        y_positions = results_drift[results_drift[:, 1] == roi_index + 1, 3]

        plt.cla()
        ax.scatter(x_positions, y_positions)
        if nm_or_pixels == 'nm':
            ax.set_xlabel('x-position (nm)')
            ax.set_ylabel('y-position (nm)')
        else:
            ax.set_xlabel('x-position (pixels)')
            ax.set_ylabel('y-position (pixels)')
        ax.set_title('Scatter ROI ' + str(roi_index + 1) + " post drift")
        name = path + "/" + str(roi_index + 1) + "_scatter_post_drift.png"
        fig.savefig(name, bbox_inches='tight')
