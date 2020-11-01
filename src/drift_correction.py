# -*- coding: utf-8 -*-
"""
Created on Thu 30-07-2020

@author: Dion Engels
PLASMON Data Analysis

drift_correction

This package is for the drift correction of PLASMON.

----------------------------

v0.1: drift correction v1: 31/07/2020
v0.1.1: bug fix and save drift: 03/08/2020
v1.0: more output just after initial release: 07/08/2020
v1.1: switch to Python coordinate system: 10/08/2020
v2.0: part of v2.0: 03/10/2020

"""

import numpy as np
from scipy.stats import norm

__self_made__ = True


class DriftCorrector:
    """
    Drift correction class of PLASMON. Takes results and corrects them for drift
    """
    def __init__(self, method):
        """
        Initialisation, does not do much
        ----------------------
        :param method: method used to get results
        """
        self.threshold_sigma = 5
        self.method = method

    def main(self, rois, name_dataset, n_frames):
        """
        Main, put in results and get out drift corrected results
        --------------------
        :param rois: all ROIs
        :param name_dataset: name of dataset that drift correction is to be done for
        :param n_frames: number of frames fitted
        :return: results_drift: drift corrected results
        """
        np.warnings.filterwarnings('ignore')  # ignore warnings of "nan" values to a real value

        # declare
        all_drift_x = np.zeros((n_frames, len(rois)))
        all_drift_y = np.zeros((n_frames, len(rois)))

        for roi_index, roi in enumerate(rois):
            # get drift for each ROI
            roi_drift_x, roi_drift_y, roi.results[name_dataset]['event_or_not'] = \
                self.find_drift(roi.results[name_dataset]['result'])
            all_drift_x[:, roi_index] = roi_drift_x
            all_drift_y[:, roi_index] = roi_drift_y

        # get mean drift
        mean_drift_x = np.nanmean(all_drift_x, axis=1)
        mean_drift_y = np.nanmean(all_drift_y, axis=1)

        # Set drift correct results per ROI
        for roi in rois:
            roi.results[name_dataset]['result_post_drift'], roi.results[name_dataset]['drift'] = \
                self.adjust_for_drift(mean_drift_x, mean_drift_y, roi.results[name_dataset]['result'])

    @staticmethod
    def find_first_non_nan(array):
        """
        Finds first non-NaN value in array
        ------------------
        :param array: array
        :return: index of first non-NaN value
        """
        for index, value in enumerate(array):
            if not np.isnan(value):
                return index

    def find_drift(self, roi_results):
        """
        Finds drift for a single ROI
        -----------------------
        :param roi_results: results of that ROI
        :return: roi_drift_x: drift in x-direction
        :return: roi_drift_y: drift in y-direction
        :return: event_or_not: boolean whether or not each frame is event or not
        """
        if "Gaussian" in self.method:
            # find cutoff if Gaussian method
            cutoff = self.find_cutoff(roi_results)
            event_or_not = roi_results[:, 3] > cutoff
            roi_results[event_or_not, 1:] = np.nan
        else:
            # otherwise all not
            event_or_not = [False]*roi_results.shape[0]

        # get drift

        if self.find_first_non_nan(roi_results[:, 1]) is not None:
            roi_drift_y = roi_results[:, 1] - roi_results[self.find_first_non_nan(roi_results[:, 1]), 1]
            roi_drift_x = roi_results[:, 2] - roi_results[self.find_first_non_nan(roi_results[:, 2]), 2]
        else:  # in case only NaNs at result, just pass normal results
            roi_drift_y = roi_results[:, 1]
            roi_drift_x = roi_results[:, 2]

        return roi_drift_x, roi_drift_y, event_or_not

    def find_cutoff(self, roi_results):
        """
        Find event cutoff intensity
        -----------------------
        :param roi_results: results for a single ROI
        :return: cutoff: the event cutoff intensity
        """
        int_ravel = roi_results[~np.isnan(roi_results[:, 3]), 3]
        mean = 0
        std = 0

        for _ in range(10):
            # for 10 times, fit norm to intensity and throw away outliers
            mean, std = norm.fit(int_ravel)
            int_ravel = int_ravel[int_ravel < mean + std * self.threshold_sigma]

        return mean + self.threshold_sigma * std

    @staticmethod
    def adjust_for_drift(mean_drift_x, mean_drift_y, results):
        """
        Adjust single ROI results for found drift
        -------------------------
        :param mean_drift_x: mean drift in x-direction
        :param mean_drift_y: mean drift in y-direction
        :param results: results of ROI
        :return: results_drift: drift corrected results of ROI
        :return: drift: drift of the ROI
        """
        results_drift = results.copy()

        # correct for drift
        results_drift[:, 1] -= mean_drift_y
        results_drift[:, 2] -= mean_drift_x

        # add drifts
        drift = np.stack((mean_drift_y, mean_drift_x), axis=1)

        return results_drift, drift
