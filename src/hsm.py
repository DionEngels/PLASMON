# -*- coding: utf-8 -*-
"""
Created on Tue 04/08/2020

@author: Dion Engels
MBx Python Data Analysis

hsm

This package is for the HSM part of MBx Python.

----------------------------

v0.0.1: Loading in multiple nd2, finding .mats
v0.0.2: complete but not working
v0.0.3: continued development 31/08/2020
v0.0.4: correlations working: 07/09/2020
v0.1: working: 13/09/2020
v0.1.1: in GUI
v1.0: Working as desired and as in SPectrA: 29/09/2020

"""
# General
import os
import numpy as np

from warnings import warn
from src.warnings import InputWarning

# I/O
from scipy.io import loadmat
import mat73

# Scipy for signal processing
from scipy.ndimage import median_filter
from scipy.optimize import least_squares

# Own code
import src.tt as fitting
import src.figure_making as figuring
from src.class_dataset_and_class_roi import Dataset, normxcorr2

import matplotlib.pyplot as plt
__self_made__ = True

# %% HSM Dataset v2


class HSMDataset(Dataset):
    def __init__(self, experiment, nd2, name):
        super().__init__(experiment, name)
        self.type = "HSM"
        self.frames = np.asarray(nd2)
        self.metadata = nd2.get_metadata(verbose=False)
        self.wavelengths = None
        self.correction_file = None
        self.spec_wavelength = None
        self.spec_shape = None

        self.corrected, self.frame_for_rois = self.hsm_drift(verbose=False)

    def prepare_run(self, settings):
        def check_correct_chars(string):
            chars = set(string)

            for i in range(0, 10):
                chars.discard(str(i))

            chars.discard('[')
            chars.discard(']')
            chars.discard('.')
            chars.discard(',')
            chars.discard(':')
            chars.discard(' ')

            if len(chars) > 0:
                return False
            else:
                return True

        def parse_string_to_numpy_array(string):
            wavelength_list = []
            string = string.replace('[', '')
            string = string.replace(']', '')

            split_string = string.split(',')

            for split in split_string:
                split = split.replace(' ', '')
                try:
                    split_int = int(split)
                    wavelength_list.append(split_int)
                except ValueError:
                    range_split = split.split(':')
                    range_split = list(map(int, range_split))
                    wavelength_list.extend(np.arange(range_split[0], range_split[2] + range_split[1], range_split[1]))

            return np.asarray(wavelength_list)

        check = self.experiment.proceed_question("Are you sure?",
                                                 "You cannot change settings later. "
                                                 "Are you sure everything is set up correctly?")
        if not check:
            return False

        new_name = settings.pop('name', self.name)
        self.set_name(new_name)
        self.settings = settings

        self.correction_file = settings['correction_file']
        path = os.getcwd()
        path += ("/spectral_corrections/" + settings['correction_file'] + ".mat")
        try:  # for new MATLAB versions
            correction = loadmat(path)
            self.spec_wavelength = correction['SpectralCorrection'][0][0][0][0]
            self.spec_shape = correction['SpectralCorrection'][0][0][1][0]
        except NotImplementedError:  # for old MATLAB versions
            correction = mat73.loadmat(path)
            self.spec_wavelength = correction['SpectralCorrection']['Lambda']
            self.spec_shape = correction['SpectralCorrection']['SpecShape']

        if check_correct_chars(settings['wavelengths']):
            self.wavelengths = parse_string_to_numpy_array(settings['wavelengths'])
        else:
            warn("Wavelengths input not valid", InputWarning)
            return False

    # %% Correct for drift between frames
    def hsm_drift(self, verbose=False):
        """
        Corrects the drift between the HSM frames and adds them up for a merged frame to compare to the laser frame.

        :param verbose: If true, you get figures
        :return: data_output: all the frames aligned (without background correction)
        :return: data_merged: all the frames aligned and added up (with background correction)
        """
        offset = np.zeros((self.frames.shape[0], 2))
        offset_from_zero = np.zeros((self.frames.shape[0], 2))

        frame = self.frames[0, :, :]
        data_output = np.zeros(self.frames.shape, dtype=frame.dtype)
        data_merged_helper = np.zeros(self.frames.shape, dtype=np.int32)
        data_merged = np.zeros(frame.shape, dtype=np.int32)
        img_corrected_previous = 0

        for frame_index, frame in enumerate(self.frames):
            background = median_filter(frame, size=9, mode='constant')
            frame = frame.astype(np.int16) - background
            img_corrected = np.round((frame[int(frame.shape[0] * 0.25 - 1):int(frame.shape[0] * 0.75),
                                      int(frame.shape[1] * 0.25 - 1):int(frame.shape[1] * 0.75)]),
                                     0).astype(np.int16)
            data_merged_helper[frame_index, :, :] = frame.astype(np.int32)
            if frame_index > 0:
                frame_convolution = normxcorr2(img_corrected, img_corrected_previous)
                maxima = np.transpose(np.asarray(np.where(frame_convolution == np.amax(frame_convolution))))[0]
                offset[frame_index, :] = maxima - np.asarray(img_corrected.shape) + np.asarray([1, 1])

            img_corrected_previous = img_corrected

        for index in range(self.frames.shape[0]):
            if index == 0:
                offset_from_zero[index, :] = offset[index, :].copy()
            else:
                offset_from_zero[index, :] = offset_from_zero[index - 1, :] + offset[index, :]

        offset_from_center = offset_from_zero[int(round(self.frames.shape[0] / 2, 0)), :] - offset_from_zero
        max_offset = int(np.max(abs(offset_from_center)))
        size_frame = np.asarray(frame.shape, dtype=int)
        helper_size = tuple(size_frame + max_offset * 2)

        for frame_index, frame in enumerate(self.frames):
            bg = np.mean(frame)
            helper_image = np.ones(helper_size) * bg

            shift_dist = offset_from_center[frame_index, :].astype(int)
            helper_image[max_offset - shift_dist[-2]:size_frame[-2] + max_offset - shift_dist[-2],
                         max_offset - shift_dist[-1]:size_frame[-1] + max_offset - shift_dist[-1]] = frame

            data_output[frame_index, :, :] = helper_image[max_offset:size_frame[-2] + max_offset,
                                                          max_offset:size_frame[-1] + max_offset]

            helper_image = np.zeros(helper_size, dtype=np.int32)
            helper_image[max_offset - shift_dist[-2]:size_frame[-2] + max_offset - shift_dist[-2],
                         max_offset - shift_dist[-1]:size_frame[-1] + max_offset - shift_dist[-1]] = \
            data_merged_helper[frame_index, :, :]
            data_merged += helper_image[max_offset:size_frame[-2] + max_offset, max_offset:size_frame[-1] + max_offset]

            if verbose:
                fig, ax = plt.subplots(1)
                ax.imshow(data_output[frame_index, :, :],
                          extent=[0, data_output[frame_index, :, :].shape[1],
                                  data_output[frame_index, :, :].shape[0],
                                  0], aspect='auto')
                plt.title("Frame {} shifted final".format(frame_index))
                plt.show()

        if verbose:
            fig, ax = plt.subplots(1)
            ax.imshow(data_merged, extent=[0, data_merged.shape[1], data_merged.shape[0], 0], aspect='auto')
            plt.title("Result")
            plt.show()

        return data_output, data_merged

    # %% Run
    def run(self, verbose=False):
        """
        Main of HSM. Does all the work.

        :param verbose: True if you want figures
        :return: self.hsm_result: the actual result. An array with ROI index, Lorentzian results and r-squared
        :return: intensity_result: All the intensities per ROI used to fit the lorentzian
        """
        def find_nearest(match_array, value_array, match_index):
            """
            Finds and returns the nearest match in another array
            """
            idx = (np.abs(match_array - match_index)).argmin()
            return value_array[idx]
        # find correct shape for wavelength
        shape = np.asarray([find_nearest(self.spec_wavelength, self.spec_shape, nu) for nu in self.wavelengths])
        if verbose:
            fig, ax = plt.subplots(1)
            figuring.plot_rois(ax, self.frame_for_rois, roi_locations=self.active_rois, roi_size=9,
                               roi_offset=self.roi_offset)
            plt.show()

        # prep for fitting
        roi_size = 9
        roi_size_1d = int((roi_size - 1) / 2)
        fitter = fitting.GaussianBackground({'roi_size': roi_size, 'rejection': 'Loose', 'method': "Gaussian - Fit bg"},
                                            600, 6, self)

        pos_max = roi_size
        pos_min = 0
        sig_min = 0
        sig_max = roi_size_1d + 1

        # %% Fit every ROI for every frame

        for roi_index, roi in enumerate(self.active_rois):
            raw_intensity = np.zeros(self.frames.shape[0])
            intensity = np.zeros(self.frames.shape[0])
            hsm_result = np.zeros(3)
            frame_stack = roi.get_frame_stack(self.corrected, roi_size_1d, self.roi_offset)
            for frame_index, my_roi in enumerate(frame_stack):
                if verbose:
                    fig, ax = plt.subplots(1)
                    ax.imshow(my_roi, extent=[0, my_roi.shape[1], my_roi.shape[0], 0], aspect='auto')
                    ax.set_xlabel('x (pixels)')
                    ax.set_ylabel('y (pixels)')
                    ax.set_title('Zoom-in ROI #{}, frame #{}'.format(roi.index, frame_index))
                    plt.show()
                result, _, success = fitter.fit_gaussian(my_roi, roi.index)
                if success == 0 or \
                        result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max \
                        or result[3] < sig_min or result[3] > sig_max or result[4] < sig_min or result[4] > sig_max:
                    raw_intensity[frame_index] = np.nan
                    intensity[frame_index] = np.nan
                else:
                    raw_intensity[frame_index] = 2 * np.pi * result[0] * result[3] * result[4]
                    intensity[frame_index] = raw_intensity[frame_index] / shape[frame_index]

            # %% Fit the total intensity of a single ROI over all frames with Lorentzian
            if verbose:
                fig, ax = plt.subplots(1)
                ax.plot(self.wavelengths, intensity)
                ax.set_title('Result ROI #{}'.format(roi.index))
                plt.show()

            result, r_squared = self.fit_lorentzian(intensity, self.wavelengths, verbose=verbose)

            hsm_result[0] = 1248 / result[2]  # SP lambda
            hsm_result[1] = 1000 * result[3]  # linewidth
            hsm_result[2] = r_squared

            result_dict = {"type": self.type, "result": hsm_result, "fit_parameters": result,
                           "raw_intensity": raw_intensity, "intensity": intensity, "raw": frame_stack}
            roi.results[self.name_result] = result_dict

            if roi_index % (round(len(self.active_rois) / 5, 0)) == 0:
                self.experiment.progress_updater.status(roi_index, len(self.active_rois))

    def fit_lorentzian(self, scattering, wavelength, split=False, verbose=False):
        """
        Function to fit a lorentzian to the found intensities

        :param scattering: the scattering intensities found
        :param wavelength: the wavelengths of the found intensities
        :param split: if True, this function has been called recursively with only a part of the original array
        :param verbose: if True, you get a lot of images
        :return: result: resulting Lorentzian parameters
        :return r_squared: the r-squared of this fit
        """

        def lorentzian(params, x):
            return params[0] + params[1] / ((x - params[2]) ** 2 + (0.5 * params[3]) ** 2)

        def error_func(p, x, y):
            return lorentzian(p, x) - y

        def find_r_squared(f, p, x, y):
            res = y - f(p, x)
            ss_res = np.sum(res ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            return 1 - ss_res / ss_tot

        def compare_plot(x, y, p):
            fy = lorentzian(p, x)
            fig, ax = plt.subplots(1)
            ax.plot(x, y)
            ax.plot(x, fy)
            plt.show()

        # remove nans

        to_del = ~np.isnan(scattering)
        scattering = scattering[to_del]
        wavelength = wavelength[to_del]

        # return if not enough points

        if len(scattering) < 5:
            return [np.nan, np.nan, np.nan, np.nan], 0

        # convert to eV
        wavelength_ev = 1248 / wavelength

        # find max and min

        max_sca = np.nanmax(scattering[scattering < np.nanmax(scattering)])
        idx_max = np.nanargmax(scattering[scattering < np.nanmax(scattering)])
        min_sca = np.nanmin(scattering)
        idx_min = np.nanargmin(scattering)

        # init guess and first fit

        init_1w = abs(2 / (np.pi * max_sca) * np.trapz(scattering, wavelength_ev))
        init_guess = [min_sca, min_sca * init_1w / (2 * np.pi), wavelength_ev[idx_max], init_1w]

        result_full = least_squares(error_func, init_guess, args=(wavelength_ev, scattering))
        result = result_full.x
        result[3] = abs(result[3])
        r_squared = find_r_squared(lorentzian, result, wavelength_ev, scattering)

        # if bad fit, try standard values
        if r_squared < 0.9:
            result_full_std = least_squares(error_func, [-10, 100, 1248 / wavelength_ev[idx_max], 0.15],
                                            args=(wavelength_ev, scattering))
            result_std = result_full_std.x
            result_std[3] = abs(result_std[3])
            r_squared_std = find_r_squared(lorentzian, result_std, wavelength_ev, scattering)
            if r_squared_std > r_squared:
                result = result_std
                r_squared = r_squared_std

        # random try if bad fit
        if r_squared < 0.9:
            result_full_base = least_squares(error_func, [0, 1, 600, 1], args=(wavelength_ev, scattering))
            result_base = result_full_base.x
            result_base[3] = abs(result_base[3])
            r_squared_base = find_r_squared(lorentzian, result_base, wavelength_ev, scattering)
            if r_squared_base > r_squared:
                result = result_base
                r_squared = r_squared_base

        # if r_squared is too low, split
        if r_squared < 0.9 and split is False:
            wavelength_low = wavelength[:idx_min]
            wavelength_high = wavelength[idx_min:]
            scattering_low = scattering[:idx_min]
            scattering_high = scattering[idx_min:]
            result_low, r_squared_low = self.fit_lorentzian(scattering_low, wavelength_low, split=True)
            result_high, r_squared_high = self.fit_lorentzian(scattering_high, wavelength_high, split=True)

            if r_squared_high > r_squared and ~np.isnan(np.sum(result_high)):
                result = result_high
                r_squared = r_squared_high
            if r_squared_low > r_squared and ~np.isnan(np.sum(result_low)):
                result = result_low
            r_squared = find_r_squared(lorentzian, result, wavelength_ev, scattering)

        if verbose:
            compare_plot(wavelength_ev, scattering, result)

        return result, r_squared
