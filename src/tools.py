# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
PLASMON Data Analysis

Tools

Some additional tools used by PLASMON

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
v2.0: part of GUI v2: 15/10/2020

"""
from numpy import zeros, moveaxis

__self_made__ = True


def change_to_nm(results, metadata, method):
    """
    Change the pixels to nm
    -----------------
    :param: results: the results in pixels
    :param: metadata: metadata to get the pixelsize from
    :param: method: method, to check if sigmas also need to be changed
    :return: results: the results in nm
    """
    pixelsize_nm = metadata['pixel_microns'] * 1000
    results[:, 1] *= pixelsize_nm  # y position to nm
    results[:, 2] *= pixelsize_nm  # x position to nm
    if "Gaussian" in method:
        results[:, 4] *= pixelsize_nm  # sigma y to nm
        results[:, 5] *= pixelsize_nm  # sigma x to nm

    return results


def convert_to_matlab(experiment):
    """
    Converts an experiments results to MATLAB
    -----------------------
    :param experiment: experiment to convert
    :return: None. Changes experiment object
    """
    for dataset in experiment.datasets:
        # set offset to MATLAB coordinates
        dataset.roi_offset = offset_to_matlab(dataset.roi_offset)
        for roi in dataset.active_rois:
            if dataset.type == "TT":
                # Change results to matlab, both pre- and post-drift
                roi.results[dataset.name_result]['result'] = \
                    result_to_matlab(roi.results[dataset.name_result]['result'], dataset.settings['method'])
                roi.results[dataset.name_result]['result_post_drift'] = \
                    result_to_matlab(roi.results[dataset.name_result]['result_post_drift'], dataset.settings['method'])
                # change drift and raw to matlab
                roi.results[dataset.name_result]['drift'] = switch_axis(roi.results[dataset.name_result]['drift'])
                roi.results[dataset.name_result]['raw'] = raw_to_matlab(roi.results[dataset.name_result]['raw'])
            elif dataset.type == "HSM":
                # change HSM to matlab
                roi.results[dataset.name_result]['raw'] = raw_to_matlab(roi.results[ dataset.name_result]['raw'])
            else:
                pass

    # change ROIs to matlab
    experiment.rois = roi_to_matlab(experiment.rois)


def offset_to_matlab(offset):
    """
    Change offset to MATLAB
    ---------------------
    :param offset: base offet
    :return: MATLAB offset
    """
    new = offset
    new[1] = offset[0]
    new[0] = offset[1]
    return new


def roi_to_matlab(rois):
    """
    ROIs to MATLAB. Adds one to index
    --------------------------
    :param rois: All ROIs in experiment
    :return: New ROIs
    """
    for roi in rois:
        roi.index += 1
    return rois


def result_to_matlab(results, method):
    """
    Change from Python to MATLAB coordinates for results
    -----------------
    :param: results: the results in Python coordinates
    :param: method: to check whether or not sigmas also need to be changed
    :return: results: the results switched to MATLAB coordinates
    """
    results[:, 0] += 1  # add one to frame counting

    results[:, 1:3] = switch_axis(results[:, 1:3])  # switch x-y

    if "Gaussian" in method:
        results[:, 4:6] = switch_axis(results[:, 4:6])  # switch sigma x-y if Gaussian

    return results


def raw_to_matlab(raw):
    """
    Switch axis of raw for MATLAB
    -----------------------
    :param raw: raw data
    :return: new raw data
    """
    if raw.ndim > 2:
        raw = moveaxis(raw, 0, -1)

    return raw


def switch_axis(array):
    """
    Switches a single arrays values
    ----------
    :param: array : array to switch
    :return: new : switched array
    """
    new = zeros(array.shape)
    new[:, 1] = array[:, 0]
    new[:, 0] = array[:, 1]
    return new
