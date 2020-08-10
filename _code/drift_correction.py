# -*- coding: utf-8 -*-
"""
Created on Thu 30-07-2020

@author: Dion Engels
MBx Python Data Analysis

drift_correction

This package is for the drift correction of MBx Python.

----------------------------

v0.1: drift correction v1: 31/07/2020
v0.1.1: bug fix and save drift: 03/08/2020
v1.0: more output just after initial release: 07/08/2020
v1.1: switch to Python coordinate system: 10/08/2020

"""

import numpy as np
from scipy.stats import norm


class DriftCorrector:
    """
    Drift correction class of MBx Python. Takes results and corrects them for drift
    """
    def __init__(self, method):
        """
        Initialisation, does not do much
        """
        self.threshold_sigma = 5
        self.method = method
        self.n_rois = 0
        self.n_frames = 0

    def main(self, results, rois, n_frames):
        """
        Main, put in results and get out drift corrected results

        Parameters
        ----------
        results: non-drift corrected results
        rois: list of ROIs
        n_frames: number of frames fitted

        Returns
        ----------
        results_drift: drift corrected results
        """
        np.warnings.filterwarnings('ignore')  # ignore warnings of "nan" values to a real value

        self.n_rois = rois.shape[0]
        self.n_frames = n_frames

        event_or_not_total = np.zeros((self.n_frames, self.n_rois), dtype=bool)

        all_drift_x = np.zeros((self.n_frames, self.n_rois))
        all_drift_y = np.zeros((self.n_frames, self.n_rois))

        for i in range(self.n_rois):  # find drift per ROI
            roi_results = results[results[:, 1] == i, :]
            roi_drift_x, roi_drift_y, event_or_not_roi = self.find_drift(roi_results)
            all_drift_x[:, i] = roi_drift_x
            all_drift_y[:, i] = roi_drift_y
            event_or_not_total[:, i] = event_or_not_roi
        mean_drift_x = np.nanmean(all_drift_x, axis=1)  # average drift of all ROIs
        mean_drift_y = np.nanmean(all_drift_y, axis=1)  # average drift of all ROIs

        results_drift, drift = self.adjust_for_drift(mean_drift_x, mean_drift_y, results)

        return results_drift, drift, event_or_not_total

    @staticmethod
    def find_first_non_nan(array):
        for index, value in enumerate(array):
            if not np.isnan(value):
                return index

    def find_drift(self, roi_results):

        if "Gaussian" in self.method:
            cutoff = self.find_cutoff(roi_results)
            event_or_not = roi_results[:, 4] > cutoff
            roi_results[event_or_not, 2:] = np.nan
        else:
            event_or_not = [False]*roi_results.shape[0]

        roi_drift_y = roi_results[:, 2] - roi_results[0, 2]
        roi_drift_x = roi_results[:, 3] - roi_results[0, 3]

        return roi_drift_x, roi_drift_y, event_or_not

    def find_cutoff(self, roi_results):
        int_ravel = roi_results[~np.isnan(roi_results[:, 4]), 4]
        mean = 0
        std = 0

        for _ in range(10):
            mean, std = norm.fit(int_ravel)
            int_ravel = int_ravel[int_ravel < mean + std * self.threshold_sigma]

        return mean + self.threshold_sigma * std

    def adjust_for_drift(self, mean_drift_x, mean_drift_y, results):

        results_drift = results.copy()

        mean_drift_x_repeat = np.repeat(mean_drift_x, self.n_rois)
        mean_drift_y_repeat = np.repeat(mean_drift_y, self.n_rois)

        results_drift[:, 2] -= mean_drift_y_repeat
        results_drift[:, 3] -= mean_drift_x_repeat

        drift = np.stack((mean_drift_y, mean_drift_x), axis=1)

        return results_drift, drift
