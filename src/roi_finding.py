# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:20:29 2020

@author: Dion Engels
PLASMON Data Analysis

ROI finding

This package holds all the PLASMON parts that are required to find the ROIs in the .nd2 files

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
v2.0: part of GUI v2 release; 15/10/2020

"""
import src.tt as fitting
from src.class_dataset_and_class_roi import Roi

import numpy as np  # for linear algebra
from scipy.signal import convolve2d
from scipy.ndimage import median_filter
from scipy.ndimage.filters import maximum_filter

__self_made__ = True
# %% Python ROI finder


class RoiFinder:
    """
    Class to find ROIs
    """
    def __init__(self, frame, signed_data_type, settings=None):
        """
        Initialises ROI finder
        -----------------------
        :param frame: Frame in which ROIs will be found
        :param settings: for resetting when frame has already been loaded before
        :return: None. Sets up class
        """
        self.base_frame = frame
        # standard fitter settings
        fitter_settings = {'roi_size': 7, 'method': "Gaussian", 'rejection': False}
        self.fitter = fitting.Gaussian(fitter_settings, 300, 5, [0, 0])
        self.data_type = signed_data_type

        # setup lists
        self.sigma_list = []
        self.int_list = []
        self.corr_list = []
        self.roi_locations = []

        if settings is None:
            # set standard settings if not given
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

            # correct for background
            background = median_filter(frame, size=self.filter_size)
            self.frame_bg = frame.astype(self.data_type) - background

            # find ROI locations
            self.roi_locations = self.main()
            # get lists
            self.int_sigma_limit(return_int=True, return_sigmas=True)
            # set standard sigma/int min/max
            self.sigma_max = np.max(self.sigma_list) * 1.5  # 50% margin
            self.int_max = np.max(self.int_list) * 1.5  # 50% margin
            self.int_min = np.min(self.int_list) / 2
        else:
            # just take values from settings
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
        Changed the settings of the ROI finder.
        ----------------------------------
        :param: settings: a dictionary of all settings
        :return: None. Changes class
        """
        # copy settings
        self.roi_size = settings['roi_size']
        self.roi_size_1d = int((self.roi_size - 1) / 2)
        self.side_distance = settings['roi_side']
        self.roi_distance = settings['inter_roi']

        self.corr_min = settings['corr_min']
        self.sigma_min = settings['sigma_min']
        self.sigma_max = settings['sigma_max']
        self.int_max = settings['int_max']
        self.int_min = settings['int_min']

        # if filter size not changed, just take processed frame from dict to prevent processing it again
        if settings['filter_size'] != self.filter_size:
            self.filter_size = settings['filter_size']
            background = median_filter(self.base_frame, size=self.filter_size)
            self.frame_bg = self.base_frame.astype(self.data_type) - background
        else:
            processed_frame = settings.pop('processed_frame', None)
            if processed_frame is not None:
                self.frame_bg = processed_frame
            else:
                background = median_filter(self.base_frame, size=self.filter_size)
                self.frame_bg = self.base_frame.astype(self.data_type) - background

    def get_settings(self):
        """
        Get settings from ROI finder. Saves to dict
        --------------------------
        :return: settings dictionary
        """
        return {'int_max': self.int_max, 'int_min': self.int_min,
                'sigma_min': self.sigma_min, 'sigma_max': self.sigma_max,
                'corr_min': self.corr_min, 'roi_size': self.roi_size, 'filter_size': self.filter_size,
                'roi_side': self.side_distance, 'inter_roi': self.roi_distance,
                'processed_frame': self.frame_bg}

    @staticmethod
    def make_gaussian(size, fwhm=3, center=None):
        """
        Makes a 2D Gaussian
        ----------
        :param: size : Size of Gaussian
        :param: fwhm : FWHM. The default is 3.
        :param: center : Center position of Gaussian. The default is None.
        :return: size by size array of 2D Gaussian
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
        ----------
        :param: return_corr : If true, returns the values of the correlations for graphing GUI
        :return: beads : boolean array. One if ROI position, zero if not
        :return: roi_locations. List of ROIs
        """
        roi_locations = []
        # comparison gaussian
        compare = self.make_gaussian(self.roi_size)

        # compare with gaussian
        frame_convolution = convolve2d(self.frame_bg, compare, mode='same')

        # filter for maxima
        fp = np.ones((3, 3), dtype=bool)
        local_peaks = maximum_filter(frame_convolution, footprint=fp)
        local_peaks_bool = (frame_convolution == local_peaks)

        max_convolution = np.max(frame_convolution)

        # find beads
        beads = (frame_convolution * local_peaks_bool) > self.corr_min * max_convolution

        locations = np.transpose(np.where(beads == 1))

        # find ROIs
        for roi in locations:
            roi_locations.append(Roi(roi[1], roi[0]))

        # if return corr, return convolutions
        if return_corr:
            corr = frame_convolution * local_peaks_bool / max_convolution

            for roi in roi_locations:
                value = corr[roi.y, roi.x]
                self.corr_list.append(value)

        return beads, roi_locations

    def boundary_rois(self, roi_boolean):
        """
        Takes boolean of correlation with 2D Gaussian and checks if ROIs are too close to the side of the frame
        :param roi_boolean: Boolean value whether or not pixel is defined as ROI after correlation with 2D Gaussian
        :return: None officially. Adapts ROI locations
        """
        remove_list = []

        for roi_index, roi in enumerate(self.roi_locations):
            my_roi = roi.get_roi(roi_boolean, self.side_distance, [0, 0])
            if my_roi.shape != (self.side_distance * 2 + 1, self.side_distance * 2 + 1):
                remove_list.append(roi_index)  # if this fails, the roi is on the boundary
                continue

        # remove all bad ROIs from ROI list
        self.roi_locations = [roi for roi_index, roi in enumerate(self.roi_locations) if roi_index not in remove_list]

    def int_sigma_limit(self, return_int=False, return_sigmas=False):
        """
        Checks intensity and sigma of each defined ROI. Rejects them if outside thresholds
        ----------
        :param: return_int : Boolean whether or not intensity list should be returned
        :param: return_sigmas : Boolean whether or not sigma list should be returned
        :return: None officially.
        Either adapts ROI_locations, sigma_list or int_list depending on the aforementioned booleans

        """
        remove_list = []

        for roi_index, roi in enumerate(self.roi_locations):
            # gets ROI
            my_roi = roi.get_roi(self.frame_bg, self.roi_size_1d, [0, 0])

            # fits ROI
            result, its, success = self.fitter.fit_gaussian(my_roi)

            # return sigma or int if desired
            if return_sigmas:
                self.sigma_list.append(result[4])
                self.sigma_list.append(result[3])
            if return_int:
                self.int_list.append(result[0])

            # check if sigma and int without bounds, otherwise remove
            if result[4] < self.sigma_min or result[3] < self.sigma_min:
                remove_list.append(roi_index)
            if result[4] > self.sigma_max or result[3] > self.sigma_max:
                remove_list.append(roi_index)
            if result[0] < self.int_min or result[0] > self.int_max:
                remove_list.append(roi_index)
        self.roi_locations = [roi for roi_index, roi in enumerate(self.roi_locations) if roi_index not in remove_list]

    def make_new_boolean(self):
        """
        Checks the current self.roi_locations and makes this into a boolean matrix for adjacent ROIs to use
        :return: roi_boolean: matrix with ones at ROI locations
        """
        # make empty roi boolean
        roi_boolean = np.zeros(self.frame_bg.shape, dtype=bool)
        # for every ROI, set ROI location to true
        for roi in self.roi_locations:
            roi_boolean[roi.y, roi.x] = True

        return roi_boolean

    def adjacent_rois(self, roi_boolean):
        """
        Takes boolean of correlation with 2D Gaussian and checks if ROIs are too close to each other
        ----------
        :param: roi_boolean : Boolean value whether or not pixel is defined as ROI after correlation with 2D Gaussian
        :return: None officially. Adapts ROI locations

        """
        remove_list = []

        for roi_index, roi in enumerate(self.roi_locations):
            my_roi = roi.get_roi(roi_boolean, self.roi_distance, [0, 0])
            # if other ROIs in ROI, remove
            trues_in_roi = np.transpose(np.where(my_roi == True))

            if trues_in_roi.shape[0] > 1:
                remove_list.append(roi_index)

        # remove all bad ROIs from ROI list
        self.roi_locations = [roi for roi_index, roi in enumerate(self.roi_locations) if roi_index not in remove_list]



    def main(self, return_int=False, return_sigmas=False,
             return_corr=False):
        """
        Main of ROI finder. Calls all functions above.
        ----------
        :param: return_int : Boolean whether or not intensity list is to be returned. Used by GUI to plot.
        The default is False.
        :param: return_sigmas : Boolean whether or not sigma list is to be returned. Used by GUI to plot.
        The default is False.
        :param: return_corr : Boolean whether or not correlation list is to be returned. Used by GUI to plot.
        The default is False.
        :return: Depending on booleans returns either a list or simply the ROI locations

        """
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
            self.corr_min = 0.01
            self.find_particles(return_corr)
            self.corr_min = saved_corr_min
            return self.corr_list

        # find particles
        roi_boolean, self.roi_locations = self.find_particles(return_corr)

        # reject boundary ROIs
        self.boundary_rois(roi_boolean)

        # reject with intensity / sigma
        self.int_sigma_limit(return_int=return_int, return_sigmas=return_sigmas)

        # reject if too close to each other
        roi_boolean = self.make_new_boolean()
        self.adjacent_rois(roi_boolean)

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
