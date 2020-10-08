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
v1.2: find max_its

"""
import src.tt as fitting
from src.class_dataset_and_class_roi import Roi

import numpy as np  # for linear algebra
from scipy.signal import convolve2d
from scipy.ndimage import median_filter
from scipy.ndimage.filters import maximum_filter
from math import ceil

from scipy.stats import norm

# %% Python ROI finder
__self_made__ = True


class RoiFinder:
    """
    Class to find ROIs in the first frame of a microscope video
    """
    def __init__(self, frame, settings=None):
        """
        Parameters
        ----------
        frame : Frame in which ROIs will be found
        settings : for resetting when frame has already been loaded before

        Returns
        -------
        None.

        """
        self.base_frame = frame
        fitter_settings = {'roi_size': 7, 'method': "Gaussian", 'rejection': "None"}
        self.fitter = fitting.Gaussian(fitter_settings, 300, 5, [0, 0])

        self.sigma_list = []
        self.int_list = []
        self.corr_list = []
        self.roi_locations = []

        if settings is None:
            self.filter_size = 9
            self.roi_size = 7
            self.roi_size_1d = int((self.roi_size - 1) / 2)
            self.side_distance = 11
            self.roi_distance = 6

            self.corr_min = 0.05
            self.sigma_min = 0
            self.sigma_max = np.inf
            self.int_min = 0
            self.int_max = np.inf

            background = median_filter(frame, size=self.filter_size)
            self.frame_bg = frame.astype('float') - background

            self.roi_locations = self.main()
            self.int_sigma_limit(return_int=True, return_sigmas=True)

            self.sigma_max = np.max(self.sigma_list) * 1.5  # 50% margin
            self.int_max = np.max(self.int_list) * 1.5  # 50% margin
            self.int_min = np.min(self.int_list) / 2
        else:
            self.filter_size = settings['filter_size']
            self.roi_size = settings['roi_size']
            self.roi_size_1d = int((self.roi_size - 1) / 2)
            self.side_distance = settings['roi_side']
            self.roi_distance = settings['inter_roi']

            self.corr_min = settings['corr_min']
            self.sigma_min = settings['sigma_min']
            self.sigma_max = settings['sigma_max']
            self.int_max = settings['int_max']
            self.int_min = settings['int_min']

            self.frame_bg = settings['processed_frame']

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
        self.roi_size = settings['roi_size']
        self.roi_size_1d = int((self.roi_size - 1) / 2)
        self.side_distance = settings['roi_side']
        self.roi_distance = settings['inter_roi']

        self.corr_min = settings['corr_min']
        self.sigma_min = settings['sigma_min']
        self.sigma_max = settings['sigma_max']
        self.int_max = settings['int_max']
        self.int_min = settings['int_min']

        if settings['filter_size'] != self.filter_size:
            self.filter_size = settings['filter_size']
            background = median_filter(self.base_frame, size=self.filter_size)
            self.frame_bg = self.base_frame.astype('float') - background
        else:
            processed_frame = settings.pop('processed_frame', None)
            if processed_frame is not None:
                self.frame_bg = processed_frame
            else:
                background = median_filter(self.base_frame, size=self.filter_size)
                self.frame_bg = self.base_frame.astype('float') - background

    def get_settings(self):
        return {'int_max': self.int_max, 'int_min': self.int_min,
                'sigma_min': self.sigma_min, 'sigma_max': self.sigma_max,
                'corr_min': self.corr_min, 'roi_size': self.roi_size, 'filter_size': self.filter_size,
                'roi_side': self.side_distance, 'inter_roi': self.roi_distance,
                'processed_frame': self.frame_bg}

    @staticmethod
    def make_gaussian(size, fwhm=3, center=None):
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
        """
        roi_locations = []
        compare = self.make_gaussian(self.roi_size)

        frame_convolution = convolve2d(self.frame_bg, compare, mode='same')

        fp = np.ones((3, 3), dtype=bool)
        local_peaks = maximum_filter(frame_convolution, footprint=fp)
        local_peaks_bool = (frame_convolution == local_peaks)

        max_convolution = np.max(frame_convolution)

        beads = (frame_convolution * local_peaks_bool) > self.corr_min * max_convolution

        locations = np.transpose(np.where(beads == 1))

        for roi in locations:
            roi_locations.append(Roi(roi[1], roi[0]))

        if return_corr:
            corr = frame_convolution * local_peaks_bool / max_convolution

            for roi in roi_locations:
                value = corr[roi.y, roi.x]
                self.corr_list.append(value)

        return beads, roi_locations

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
        remove_list = []

        for roi_index, roi in enumerate(self.roi_locations):
            my_roi = roi.get_roi(roi_boolean, self.side_distance, [0, 0])
            if my_roi.shape != (self.side_distance * 2 + 1, self.side_distance * 2 + 1):
                remove_list.append(roi_index)  # if this fails, the roi is on the boundary
                continue

            my_roi = roi.get_roi(roi_boolean, self.roi_distance, [0, 0])

            trues_in_roi = np.transpose(np.where(my_roi == True))

            if trues_in_roi.shape[0] > 1:
                remove_list.append(roi_index)

        self.roi_locations = [roi for roi_index, roi in enumerate(self.roi_locations) if roi_index not in remove_list]

    def int_sigma_limit(self, return_int=False, return_sigmas=False):
        """
        Checks intensity and sigma of each defined ROI. Rejects them if outside thresholds

        Parameters
        ----------
        return_int : Boolean whether or not intensity list should be returned
        return_sigmas : Boolean whether or not sigma list should be returned

        Returns
        -------
        None officially. Either adapts ROI_locations, sigma_list or int_list depending on the afforementioned booleans

        """
        remove_list = []

        for roi_index, roi in enumerate(self.roi_locations):
            my_roi = roi.get_roi(self.frame_bg, self.roi_size_1d, [0, 0])

            result, its, success = self.fitter.fit_gaussian(my_roi, roi_index)

            if return_sigmas:
                self.sigma_list.append(result[4])
                self.sigma_list.append(result[3])
            if return_int:
                self.int_list.append(result[0])

            if result[4] < self.sigma_min or result[3] < self.sigma_min:
                remove_list.append(roi_index)
            if result[4] > self.sigma_max or result[3] > self.sigma_max:
                remove_list.append(roi_index)
            if result[0] < self.int_min or result[0] > self.int_max:
                remove_list.append(roi_index)
        self.roi_locations = [roi for roi_index, roi in enumerate(self.roi_locations) if roi_index not in remove_list]

    def main(self, return_int=False, return_sigmas=False,
             return_corr=False):
        """
        Main of ROI finder. Calls all functions above.

        Parameters
        ----------
        return_int : Boolean whether or not intensity list is to be returned. Used by GUI to plot. The default is False.
        return_sigmas : Boolean whether or not sigma list is to be returned. Used by GUI to plot. The default is False.
        return_corr : Boolean whether or not correlation list is to be returned. Used by GUI to plot.
        The default is False.

        Returns
        -------
        Depending on booleans returns either a list or simply the ROI locations

        """

        saved_corr_min = None

        # directly return if already made lists
        if return_corr and self.corr_list != []:
            return self.corr_list
        elif return_int and self.int_list != []:
            return self.int_list
        elif return_sigmas and self.sigma_list != []:
            return self.sigma_list

        # for corr list, save actual corr and use a different one for now
        if return_corr:
            saved_corr_min = self.corr_min
            self.corr_min = 0.005
            self.find_particles(return_corr)
            self.corr_min = saved_corr_min
            return self.corr_list

        # find particles
        roi_boolean, self.roi_locations = self.find_particles(return_corr)

        # continue finding particles
        self.adjacent_or_boundary_rois(roi_boolean)

        self.int_sigma_limit(return_int=return_int, return_sigmas=return_sigmas)

        # set final ROI number
        for roi_index, roi in enumerate(self.roi_locations):
            roi.set_index(roi_index)

        # return correct thing
        if return_sigmas:
            return self.sigma_list
        elif return_int:
            return self.int_list
        else:
            return self.roi_locations

    def find_snr(self):
        intensity_list = self.main(return_int=True)

        mu, std = norm.fit(intensity_list)
        intensity_list2 = np.asarray(intensity_list)[intensity_list < mu]
        mu, std = norm.fit(intensity_list2)

        if mu >= 2000:
            max_its = 100
        else:
            int_under_which_more_its_are_needed = 2000
            max_its = ceil((int_under_which_more_its_are_needed - mu) / 1000) * 100 + 100

        return max_its

    def find_snr_load_from_other(self, roi_locations):

        int_list = []

        for roi_index, roi in enumerate(roi_locations):
            y = int(roi[0])
            x = int(roi[1])

            my_roi = roi.get_roi(self.frame_bg, self.roi_size_1d, [0, 0])
            result, its, success = self.fitter.fit_gaussian(my_roi, roi_index)

            int_list.append(result[0])

        mu, std = norm.fit(int_list)
        intensity_list2 = np.asarray(int_list)[int_list < mu]
        mu, std = norm.fit(intensity_list2)

        if mu >= 2000:
            max_its = 100
        else:
            int_under_which_more_its_are_needed = 2000
            max_its = ceil((int_under_which_more_its_are_needed - mu) / 1000) * 100 + 100

        return max_its
