# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:20:29 2020

@author: Dion Engels
MBx Python Data Analysis

ROI finding

This package holds all the MBx Python parts that are required to find the ROIs in the .nd2 files

----------------------------

v0.1.0, ROI detection:  31/05/2020
v0.1.1, conventional naming: 04/06/2020
v0.2.0: self-made ROI finding: 10/07/2020
v0.3.0: clean up and based on SPectrA; correlation, pixel_int, sigma and int
v0.4.0: bug fix ROI distance, ready for Peter review
v0.4.1: clean up
v0.4.2: changed "change_settings"
v0.5: removed pixel min
v1.0: bugfixes
v1.1: list creation bugfixes

"""
import numpy as np  # for linear algebra
from scipy.signal import convolve2d
from scipy.ndimage import median_filter
from scipy.ndimage.filters import maximum_filter

# %% Python ROI finder


class RoiFinder:
    """
    Class to find ROIs in the first frame of a microscope video
    """
    def __init__(self, frame, fitter, filter_size=9, corr_min=0.05,
                 sigma_min=0, sigma_max=None, int_min=None, int_max=None,
                 roi_size=7, settings=None):
        """
        Parameters
        ----------
        frame : Frame in which ROIs will be found
        fitter : Fitter to fit ROIs with. Is Gaussian fitter.
        filter_size : optional, size of filter in pixels. The default is 9.
        corr_min : optional, minimum correlation value threshold. The default is 0.05.
        sigma_min : optional, minimum sigma value threshold. The default is 0.
        sigma_max : optional, maximum sigma value threshold. The default is None.
        int_min : optional, minimum intensity value threshold. The default is None.
        int_max : optional, maximum intensity value threshold. The default is None.
        roi_size : optional, roi size.  The default is 7.
        settings : for resetting when frame has already been loaded before

        Returns
        -------
        None.

        """
        self.sigma_list = []
        self.int_list = []
        self.corr_list = []
        self.roi_locations = []

        self.base_frame = frame

        if settings is None:
            self.filter_size = int(filter_size)
            self.roi_size = roi_size
            self.roi_size_1d = int((self.roi_size - 1) / 2)
            self.side_distance = 11
            self.roi_distance = 6

            background = median_filter(frame, size=self.filter_size, mode='constant')
            background[background == 0] = np.min(background[background > 0])
            self.frame = frame.astype('float') - background

            self.corr_min = corr_min
            self.int_min = int_min
            self.sigma_min = sigma_min

            if sigma_max is None or int_max is None or int_min is None:
                self.sigma_max = 5
                self.int_max = np.inf
                self.int_min = 0
                self.roi_locations = self.main(fitter)
                if int_max is None or int_min is None:
                    self.int_sigma_limit(fitter, True, False)
                    if int_max is None:
                        self.int_max = np.max(self.int_list) * 1.5  # 50% margin
                    if int_min is None:
                        self.int_min = np.min(self.int_list) / 2
                if sigma_max is None:
                    self.int_sigma_limit(fitter, False, True)
                    self.sigma_max = np.max(self.sigma_list) * 1.5  # 50% margin
            else:
                self.sigma_max = sigma_max
                self.int_max = int_max
        else:
            self.sigma_min = settings['sigma_min']
            self.sigma_max = settings['sigma_max']

            self.corr_min = settings['corr_min']
            self.int_max = settings['int_max']
            self.int_min = settings['int_min']

            self.roi_size = settings['roi_size']
            self.roi_size_1d = int((self.roi_size - 1) / 2)

            self.side_distance = settings['roi_side']
            self.roi_distance = settings['inter_roi']
            self.filter_size = settings['filter_size']

            self.frame = settings['processed_frame']

    def change_settings(self, settings):
        """
        Changed the settings of the ROI finder. Only called by GUI.

        Parameters
        ----------
        settings: a dictionary of all settings

        Returns
        -------
        None.

        """
        self.sigma_min = settings['sigma_min']
        self.sigma_max = settings['sigma_max']

        self.corr_min = settings['corr_min']
        self.int_max = settings['int_max']
        self.int_min = settings['int_min']

        self.roi_size = settings['roi_size']
        self.roi_size_1d = int((self.roi_size - 1) / 2)

        self.side_distance = settings['roi_side']
        self.roi_distance = settings['inter_roi']

        if settings['filter_size'] != self.filter_size:
            self.filter_size = settings['filter_size']
            background = medfilt(self.base_frame, kernel_size=self.filter_size)
            background[background == 0] = np.min(background[background > 0])
            self.frame = self.base_frame.astype('float') - background

    def make_gaussian(self, size, fwhm=3, center=None):
        """
        Makes a 2D Gaussian

        Parameters
        ----------
        size : Size of Gaussian
        fwhm : FWHM. The default is 3.
        center : Center position of Gaussian. The default is None.

        Returns
        -------
        size by size array of 2D Gaussian

        """
        x = np.arange(0, size, 1, float)
        y = x[:, np.newaxis]

        if center is None:
            x0 = y0 = size // 2
        else:
            x0 = center[0]
            y0 = center[1]

        return np.exp(-4 * np.log(2) * ((x - x0) ** 2 + (y - y0) ** 2) / fwhm ** 2)

    def find_particles(self, return_corr):
        """
        Finds particles using correlation with 2D Gaussian

        Parameters
        ----------
        return_corr : If true, returns the values of the correlations for graphing GUI

        Returns
        -------
        beads : boolean array. One if ROI position, zero if not
        locations : List of ROI locations

        """

        compare = self.make_gaussian(self.roi_size)

        frame_convolution = convolve2d(self.frame, compare, mode='same')

        fp = np.ones((3, 3), dtype=bool)
        local_peaks = maximum_filter(frame_convolution, footprint=fp)
        local_peaks_bool = (frame_convolution == local_peaks)

        max_convolution = np.max(frame_convolution)

        beads = (frame_convolution * local_peaks_bool) > self.corr_min * max_convolution

        locations = np.transpose(np.where(beads == 1))

        if return_corr:
            corr = frame_convolution * local_peaks_bool / max_convolution

            for roi in locations:
                y = int(roi[0])
                x = int(roi[1])

                value = corr[y, x]
                self.corr_list.append(value)

        return beads, locations

    def adjacent_or_boundary_rois(self, roi_boolean):
        """
        Takes boolean of correlation with 2D Gaussian and checks if ROIs are too close to each other or side of frame

        Parameters
        ----------
        roi_boolean : Boolean value whether or not pixel is defined as ROI after correlation with 2D Gaussian

        Returns
        -------
        None officially. Adapts ROI locations

        """
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)

        for roi_index, roi in enumerate(self.roi_locations):

            y = int(roi[0])
            x = int(roi[1])

            my_roi = roi_boolean[y - self.side_distance:y + self.side_distance + 1,
                     x - self.side_distance:x + self.side_distance + 1]
            if my_roi.shape != (self.side_distance * 2 + 1, self.side_distance * 2 + 1):
                keep_boolean[roi_index] = False  # if this fails, the roi is on the boundary
                continue

            my_roi = roi_boolean[y - self.roi_distance:y + self.roi_distance + 1,
                     x - self.roi_distance:x + self.roi_distance + 1]

            trues_in_roi = np.transpose(np.where(my_roi == True))

            if trues_in_roi.shape[0] > 1:
                keep_boolean[roi_index] = False

        self.roi_locations = self.roi_locations[keep_boolean, :]

    def int_sigma_limit(self, fitter, return_int, return_sigmas):
        """
        Checks intensity and sigma of each defined ROI. Rejects them if outside thresholds

        Parameters
        ----------
        fitter : Gaussian fitter object
        return_int : Boolean whether or not intensity list should be returned
        return_sigmas : Boolean whether or not sigma list should be returned

        Returns
        -------
        None officially. Either adapts ROI_locations, sigma_list or int_list depending on the afforementioned booleans

        """
        keep_boolean = np.ones(self.roi_locations.shape[0], dtype=bool)

        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1])

            my_roi = self.frame[y - self.roi_size_1d:y + self.roi_size_1d + 1,
                     x - self.roi_size_1d:x + self.roi_size_1d + 1]

            result, its, success = fitter.fit_gaussian(my_roi, roi_index)

            if return_sigmas:
                self.sigma_list.append(result[4])
                self.sigma_list.append(result[3])
            elif return_int:
                self.int_list.append(result[0])

            if result[4] < self.sigma_min or result[3] < self.sigma_min:
                keep_boolean[roi_index] = False
            if result[4] > self.sigma_max or result[3] > self.sigma_max:
                keep_boolean[roi_index] = False
            if result[0] < self.int_min or result[0] > self.int_max:
                keep_boolean[roi_index] = False

        self.roi_locations = self.roi_locations[keep_boolean, :]

    def main(self, fitter, return_int=False, return_sigmas=False,
             return_corr=False):
        """
        Main of ROI finder. Calls all functions above.

        Parameters
        ----------
        fitter : Gaussian fitter used to find sigma and intensity.
        return_int : Boolean whether or not intensity list is to be returned. Used by GUI to plot. The default is False.
        return_sigmas : Boolean whether or not sigma list is to be returned. Used by GUI to plot. The default is False.
        return_corr : Boolean whether or not correlation list is to be returned. Used by GUI to plot.
        The default is False.

        Returns
        -------
        Depending on booleans returns either a list or simply the ROI locations

        """

        if return_corr and self.corr_list != []:
            return self.corr_list
        elif return_int and self.int_list != []:
            return self.int_list
        elif return_sigmas and self.sigma_list != []:
            return self.sigma_list

        if return_corr:
            corr_min = self.corr_min
            self.corr_min = 0.005

        roi_boolean, self.roi_locations = self.find_particles(return_corr)

        if return_corr:
            self.corr_min = corr_min
            return self.corr_list

        self.adjacent_or_boundary_rois(roi_boolean)

        self.int_sigma_limit(fitter, return_int, return_sigmas)

        if return_sigmas:
            return self.sigma_list
        elif return_int:
            return self.int_list
        else:
            return self.roi_locations
