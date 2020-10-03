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
v2.0: part of v2.0: 03/10/2020

"""

import numpy as np
from scipy.stats import norm

__self_made__ = True


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

    def main(self, rois, name_dataset, n_frames):
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

        all_drift_x = np.zeros((n_frames, len(rois)))
        all_drift_y = np.zeros((n_frames, len(rois)))

        for roi_index, roi in enumerate(rois):
            roi_drift_x, roi_drift_y, roi.results[name_dataset]['event_or_not'] = \
                self.find_drift(roi.results[name_dataset]['result'])
            all_drift_x[:, roi_index] = roi_drift_x
            all_drift_y[:, roi_index] = roi_drift_y

        mean_drift_x = np.nanmean(all_drift_x, axis=1)  # average drift of all ROIs
        mean_drift_y = np.nanmean(all_drift_y, axis=1)

        for roi in rois:
            roi.results[name_dataset]['result_post_drift'], roi.results[name_dataset]['drift'] = \
                self.adjust_for_drift(mean_drift_x, mean_drift_y, roi.results[name_dataset]['result'])

    @staticmethod
    def find_first_non_nan(array):
        for index, value in enumerate(array):
            if not np.isnan(value):
                return index

    def find_drift(self, roi_results):

        if "Gaussian" in self.method:
            cutoff = self.find_cutoff(roi_results)
            event_or_not = roi_results[:, 3] > cutoff
            roi_results[event_or_not, 1:] = np.nan
        else:
            event_or_not = [False]*roi_results.shape[0]

        roi_drift_y = roi_results[:, 1] - roi_results[self.find_first_non_nan(roi_results[:, 1]), 1]
        roi_drift_x = roi_results[:, 2] - roi_results[self.find_first_non_nan(roi_results[:, 2]), 2]

        return roi_drift_x, roi_drift_y, event_or_not

    def find_cutoff(self, roi_results):
        int_ravel = roi_results[~np.isnan(roi_results[:, 3]), 3]
        mean = 0
        std = 0

        for _ in range(10):
            mean, std = norm.fit(int_ravel)
            int_ravel = int_ravel[int_ravel < mean + std * self.threshold_sigma]

        return mean + self.threshold_sigma * std

    @staticmethod
    def adjust_for_drift(mean_drift_x, mean_drift_y, results):

        results_drift = results.copy()

        results_drift[:, 1] -= mean_drift_y
        results_drift[:, 2] -= mean_drift_x

        drift = np.stack((mean_drift_y, mean_drift_x), axis=1)

        return results_drift, drift
