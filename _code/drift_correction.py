# -*- coding: utf-8 -*-
"""
Created on Thu 30-07-2020

@author: Dion Engels
MBx Python Data Analysis

drift_correction

This package is for the drift correction of MBx Python.

----------------------------

v0.1:

"""

import numpy as np


class DriftCorrector:
    """
    Drift correction class of MBx Python. Takes results and corrects them for drift
    """
    def __init__(self):
        """
        Initialisation, does not do much
        """
        self.threshold_sigma = 5

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
        n_rois = rois.shape[0]

        all_drift_x = np.zeros((n_frames, n_rois))
        all_drift_y = np.zeros((n_frames, n_rois))

        for i in range(1, n_rois+1):
            roi_results = results[results[:, 1] == i, :]
            roi_drift_x, roi_drift_y = self.find_drift(roi_results)
            all_drift_x[:, i-1] = roi_drift_x  # -1 since converted to MATLAB counting
            all_drift_y[:, i-1] = roi_drift_y  # -1 since converted to MATLAB counting
        mean_drift_x = np.nanmean(all_drift_x, axis=1)
        mean_drift_y = np.nanmean(all_drift_y, axis=1)

        results_drift = results.copy()

        for i in range(n_frames):
            results_drift[results_drift[:, 0] == i+1, 2] -= mean_drift_x[i]  # +1 since converted to MATLAB counting
            results_drift[results_drift[:, 0] == i+1, 3] -= mean_drift_y[i]  # +1 since converted to MATLAB counting

        return results_drift

    @staticmethod
    def find_drift(roi_results):

        roi_drift_x = roi_results[:, 2] - roi_results[0, 2]
        roi_drift_y = roi_results[:, 3] - roi_results[0, 3]

        return roi_drift_x, roi_drift_y
