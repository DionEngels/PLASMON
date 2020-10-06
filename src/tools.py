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
v0.7: metadata v2
v0.8: metadata v3
v1.0: bug fixes
v1.0.1: metadata goldmine found, to be implemented
v1.1: emptied out: 09/08/2020

"""
from numpy import zeros

__self_made__ = True


def change_to_nm(results, metadata, method):
    """
    Change the pixels to nm

    Parameters
    -----------------
    results: the results in pixels
    metadata: metadata to get the pixelsize from
    method: method, to check if sigmas also need to be changed

    Returns
    ------------------
    results: the results in nm

    """
    pixelsize_nm = metadata['pixel_microns'] * 1000
    results[:, 1] *= pixelsize_nm  # y position to nm
    results[:, 2] *= pixelsize_nm  # x position to nm
    if "Gaussian" in method:
        results[:, 4] *= pixelsize_nm  # sigma y to nm
        results[:, 5] *= pixelsize_nm  # sigma x to nm

    return results


def convert_to_matlab(experiment):

    for dataset in experiment.datasets:
        dataset.roi_offset = offset_to_matlab(dataset.roi_offset)
        for roi in dataset.active_rois:
            if dataset.type == "TT":
                roi.results[dataset.name_result]['result'] = \
                    result_to_matlab(roi.results[dataset.name_result]['result'],
                                     dataset.frame_for_rois.shape[0], dataset.settings['method'],
                                     dataset.settings['pixels_or_nm'], dataset.metadata)
            elif dataset.type == "HSM":
                pass  # no correction needed for HSM
            else:
                pass

    for roi in experiment.rois:
        roi.y, roi.index = roi_to_matlab(roi, experiment.frame_for_rois.shape[0])


def offset_to_matlab(offset):
    new = offset
    new[1] = offset[0]
    new[0] = offset[1]
    new[1] = -new[1]

    return new


def roi_to_matlab(roi, height):
    """
    Change from Python to MATLAB coordinates for ROI locations

    Parameters
    -----------------
    roi: a single ROI
    height: Height of microscope view

    Returns
    ------------------
    results: the roi.y coordinate switched to MATLAB coordinates
    """
    roi.y = height - roi.y

    return roi.y, roi.index + 1


def result_to_matlab(results, height, method, nm_or_pixels, metadata):
    """
    Change from Python to MATLAB coordinates for results

    Parameters
    -----------------
    results: the results in Python coordinates
    height: Height of microscope view
    method: to check whether or not sigmas also need to be changed
    nm_or_pixels: whether or not the results are in pixels or nm
    metadata: to get pixelsize

    Returns
    ------------------
    results: the results switched to MATLAB coordinates
    """
    results[:, 0] += 1  # add one to frame counting

    results[:, 1:3] = switch_axis_to_matlab_coordinates(results[:, 1:3], height, nm_or_pixels, metadata)  # switch x-y

    if "Gaussian" in method:
        results[:, 4:6] = switch_axis(results[:, 4:6])  # switch sigma x-y if Gaussian

    return results


def switch_axis_to_matlab_coordinates(array, height, nm_or_pixels="pixels", metadata=None):
    """
    Change from Python to MATLAB coordinates for an x,y system

    Parameters
    -----------------
    array: the array to be switched
    height: Height of microscope view
    nm_or_pixels: whether or not the results are in pixels or nm
    metadata: to get pixelsize

    Returns
    ------------------
    array: the array to MATLAB coordinates
    """
    array = switch_axis(array)
    if nm_or_pixels == "nm":
        pixelsize_nm = metadata['pixel_microns'] * 1000
        array[:, 1] = height*pixelsize_nm - array[:, 1]
    else:
        array[:, 1] = height - array[:, 1]

    return array


def switch_axis(array):
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

# OLD


def roi_to_python_coordinates(roi_locs, height):
    """
    Change from MATLAB to Python coordinates for ROI locations

    Parameters
    -----------------
    roi_locs: ROI locations
    height: Height of microscope view

    Returns
    ------------------
    results: the roi locations switched to MATLAB coordinates
    """
    roi_locs[:, 1] = height - roi_locs[:, 1]
    roi_locs = switch_axis(roi_locs)

    return roi_locs
