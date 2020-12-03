# -*- coding: utf-8 -*-
"""
Created on Tue 04/08/2020

@author: Dion Engels
PLASMON Data Analysis

hsm

This package is for the HSM part of PLASMON.

----------------------------

v0.0.1: Loading in multiple nd2, finding .mats
v0.0.2: complete but not working
v0.0.3: continued development 31/08/2020
v0.0.4: correlations working: 07/09/2020
v0.1: working: 13/09/2020
v0.1.1: in GUI
v1.0: Working as desired and as in SPectrA: 29/09/2020
v2.0: Completed for v2 of program: 15/10/2020

"""
# General
import os
import numpy as np

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

# %% HSM Fitter


class HSMFit(fitting.GaussianBackground):

    def __init__(self, roi_size_1d):
        super().__init__({'roi_size': int(roi_size_1d * 2 + 1), 'rejection': True, 'method': "Gaussian - Fit bg"},
                         1000, 6, [0, 0])
        self.init_sig = 0.8

    def fit_gaussian(self, data):
        """
        Gathers parameter estimate and calls fit.

        Parameters
        ----------
        data : ROI pixel values

        Returns
        -------
        p.x: solution of parameters
        p.nfev: number of iterations
        p.success: success or failure

        """
        def extract_data(p=None):
            if p is None:
                return np.ones(5), None, False
            else:
                return p.x, p.nfev, p.success
        # set bounds
        pos_max, pos_min, int_max, int_min, sig_max, sig_min = self.define_fitter_bounds()

        # set parameters
        background = self.fun_calc_bg(data)
        height = data[int(4.5), int(4.5)] - background
        if height < 0:
            height = 0
        params = np.array([height, 4.5, 4.5, self.init_sig, self.init_sig, background])

        # try first fit
        try:
            p = least_squares(lambda p: self.fun_gaussian(p, data), params, method='dogbox', max_nfev=self.max_its,
                              ftol=10e-12, gtol=10e-12, xtol=10e-12,
                              bounds=([int_min, pos_min, pos_min, sig_min, sig_min, 0],
                                      [int_max, pos_max, pos_max, sig_max, sig_max, np.inf]))
            res, its, success = extract_data(p)
        except Exception as e:
            # if exception, print and set comparison
            print(e)
            print(np.array([height, 4.5, 4.5, self.init_sig, self.init_sig, background]))
            # extract data
            res, its, success = extract_data(p=None)

        # if fit is bad, do new fit with Cauchy loss (better at low SNR)
        if success == 0 or \
                res[2] < pos_min or res[2] > pos_max or res[1] < pos_min or res[1] > pos_max or \
                res[0] <= int_min or res[0] > int_max or \
                res[3] <= sig_min or res[3] >= sig_max or res[4] <= sig_min or res[4] >= sig_max:
            params = np.array([height, 4.5, 4.5, self.init_sig, self.init_sig, background])
            try:
                p = least_squares(lambda p: self.fun_gaussian(p, data), params, method='dogbox', max_nfev=self.max_its,
                                  ftol=10e-12, gtol=10e-12, xtol=10e-12, loss='cauchy', f_scale=0.1,
                                  bounds=([int_min, pos_min, pos_min, sig_min, sig_min, 0],
                                          [int_max, pos_max, pos_max, sig_max, sig_max, np.inf]))
                # extract data
                res, its, success = extract_data(p)
            except Exception as e:
                # if exception, print
                print(e)
                print(np.array([height, 4.5, 4.5, self.init_sig, self.init_sig, background]))
        return [res, its, success]

    def define_fitter_bounds(self):
        """
        Defines fitter bounds, based on "Strict" bounds or not

        Returns
        -------
        pos_max : Max position of fit
        pos_min : Min position of fit
        int_max : Max intensity of fit
        int_min : Min intensity of fit
        sig_max : Max sigma of fit
        sig_min : Min sigma of fit

        """
        pos_max = self.roi_size
        pos_min = 0.0
        int_min = 0.0
        int_max = np.inf
        sig_min = 0.0
        sig_max = 2.0

        return pos_max, pos_min, int_max, int_min, sig_max, sig_min

    def fitter(self, frame_stack, shape, energy_width, *_):
        """
        Does Gaussian fitting for all frames for a single ROI
        --------------------------------------------------------
        :param frame_stack: HSM corrected frame stack for one ROI
        :return: HSM results for a single ROI
        """
        # predefine
        raw_intensity = np.zeros(frame_stack.shape[0])
        intensity = np.zeros(frame_stack.shape[0])
        raw_fits = np.zeros((frame_stack.shape[0], self.num_fit_params))

        pos_max, pos_min, int_max, int_min, sig_max, sig_min = self.define_fitter_bounds()

        # fit per frame
        for frame_index, my_roi in enumerate(frame_stack):  
            # fit
            result, _, success = self.fit_gaussian(my_roi)
            # save fittings if success
            if success == 0 or \
                    result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max or \
                    result[0] <= int_min or result[0] > int_max or \
                    result[3] <= sig_min or result[3] >= sig_max or result[4] <= sig_min or result[4] >= sig_max:
                raw_intensity[frame_index] = np.nan
                intensity[frame_index] = np.nan
                raw_fits[frame_index, :] = np.nan
            else:
                raw_intensity[frame_index] = 2 * np.pi * result[0] * result[3] * result[4]
                # for intensity, divide by shape correction and energy_width normalization
                intensity[frame_index] = raw_intensity[frame_index] / shape[frame_index] / energy_width[frame_index]
                raw_fits[frame_index, :] = result

        # reject very high sigma fits (50% above average)
        intensity[(raw_fits[:, 3] > np.nanmean(raw_fits[:, 3]) * 1.5) |
                  (raw_fits[:, 4] > np.nanmean(raw_fits[:, 4]) * 1.5)] = np.nan
        raw_intensity[(raw_fits[:, 3] > np.nanmean(raw_fits[:, 3]) * 1.5) |
                      (raw_fits[:, 4] > np.nanmean(raw_fits[:, 4]) * 1.5)] = np.nan
        raw_fits[(raw_fits[:, 3] > np.nanmean(raw_fits[:, 3]) * 1.5) |
                 (raw_fits[:, 4] > np.nanmean(raw_fits[:, 4]) * 1.5), :] = np.nan

        return raw_intensity, intensity, raw_fits

# %% HSM Dataset v2


class HSMDataset(Dataset):
    """
    HSM Dataset. Inherits from base dataset
    """
    def __init__(self, experiment, nd2, name, label=None):
        """
        Initialise HSM Dataset. Set base values and load nd2. Also create corrected frame.
        ------------------------
        :param experiment: parent experiment
        :param nd2: nd2 of HSM
        :param name: name of HSM
        :param label: a label that you can add. If added, update percentages will be placed there
        """
        super().__init__(experiment, nd2, name)
        self.type = "HSM"
        self.frames = np.asarray(nd2)
        self.metadata = nd2.get_metadata(verbose=False)
        self.wavelengths = None
        self.correction_file = None
        self.spec_wavelength = None
        self.spec_shape = None

        # create corrected merged frame and corrected frames
        self.corrected, self.frame_for_rois = self.hsm_drift(verbose=False, label=label)

    def prepare_run(self, settings):
        """
        Prepare the HSM dataset for run. Takes all settings and puts them in
        --------------------------
        :param settings: settings given by user
        :return: status: boolean whether or not success. Mostly edits class though.
        """
        def check_correct_chars(string):
            """
            Check if string for wavelengths given is correct
            ---------------------------
            :param string: wavelengths string
            :return: return if valid or not
            """
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
            """
            Parse wavelength string to numpy array
            ------------------------
            :param string: wavelength string
            :return: wavelength_list as array
            """
            wavelength_list = []
            string = string.replace('[', '')
            string = string.replace(']', '')

            # split on comma
            split_string = string.split(',')

            for split in split_string:
                # for each split, get rid of spaces
                split = split.replace(' ', '')
                try:
                    # try to make int, if successful, add to list
                    split_int = int(split)
                    wavelength_list.append(split_int)
                except ValueError:
                    # if not, it is a range. Add range to list
                    range_split = split.split(':')
                    range_split = list(map(int, range_split))
                    wavelength_list.extend(np.arange(range_split[0], range_split[2] + range_split[1], range_split[1]))

            return np.asarray(wavelength_list)

        # check with user
        if not self.experiment.proceed_question("Are you sure?", "You cannot change settings later. "
                                                "Are you sure everything is set up correctly?"):
            return False

        # set new name and settings
        new_name = settings.pop('name', self.name)
        if self.check_name_validity(new_name) is False:
            self.experiment.error_func("Invalid name", "MATLAB only accepts letters, numbers, and underscores. "
                                                       "I also accept - and spaces (these will become underscores in "
                                                       "MATLAB). Please only use these.")
            return False
        self.set_name(new_name)
        self.settings = settings

        # set correction file
        self.correction_file = settings['correction_file']
        path = os.getcwd()
        path += ("/spectral_corrections/" + settings['correction_file'] + ".mat")
        # load in correction file
        try:  # for new MATLAB versions
            correction = loadmat(path)
            self.spec_wavelength = correction['SpectralCorrection'][0][0][0][0]
            self.spec_shape = correction['SpectralCorrection'][0][0][1][0]
        except NotImplementedError:  # for old MATLAB versions
            correction = mat73.loadmat(path)
            self.spec_wavelength = correction['SpectralCorrection']['Lambda']
            self.spec_shape = correction['SpectralCorrection']['SpecShape']

        # Add wavelengths. Return false if fails
        try:
            if check_correct_chars(settings['wavelengths']):
                self.wavelengths = parse_string_to_numpy_array(settings['wavelengths'])
            else:
                self.experiment.error_func("Input error", "Wavelengths input not valid")
                return False
        except:
            self.experiment.error_func("Input error", "Wavelengths input not valid")
            return False

        # check wavelengths same size as array
        if len(self.wavelengths) != self.corrected.shape[0]:
            self.experiment.error_func("Input error", "Wavelengths not same length as the amount of frames loaded")
            return False

    # %% Correct for drift between frames
    def hsm_drift(self, verbose=False, label=None):
        """
        Corrects the drift between the HSM frames and adds them up for a merged frame to compare to the laser frame.
        ----------------------------------------
        :param verbose: If true, you get figures
        :param label: a label that you can add. If added, update percentages will be placed there
        :return: data_output: all the frames aligned (without background correction)
        :return: data_merged: all the frames aligned and added up (with background correction)
        """
        # pre-declare
        offset = np.zeros((self.frames.shape[0], 2))
        offset_from_zero = np.zeros((self.frames.shape[0], 2))

        frame = self.frames[0, :, :]
        data_output = np.zeros(self.frames.shape, dtype=self.data_type)
        data_merged_helper = np.zeros(self.frames.shape, dtype=self.data_type_signed)
        data_merged = np.zeros(frame.shape, dtype=self.data_type_signed)
        img_corrected_previous = 0

        # for each frame, correct for background and save corrected frame
        for frame_index, frame in enumerate(self.frames):
            background = median_filter(frame, size=9, mode='constant')
            frame = frame.astype(self.data_type_signed) - background
            # crop 5 pixels for background correction
            img_corrected = np.round((frame[5:-5, 5:-5]), 0).astype(self.data_type_signed)
            # save to data_merged_helper to prevent doing background correction again
            data_merged_helper[frame_index, :, :] = frame.astype(self.data_type_signed)
            # after first frame, correlate with previous frame
            if frame_index > 0:
                frame_convolution = normxcorr2(img_corrected, img_corrected_previous)
                maxima = np.transpose(np.asarray(np.where(frame_convolution == np.amax(frame_convolution))))[0]
                offset[frame_index, :] = maxima - np.asarray(img_corrected.shape) + np.asarray([1, 1])

            img_corrected_previous = img_corrected
            if label is not None:
                label.updater(text=f'HSM frames are being merged. '
                                   f'Progress {(frame_index + 1) / len(self.frames) * 100:.1f}%')

        # after all correlations, add up individual offsets to get offset from frame zero
        for index in range(self.frames.shape[0]):
            if index == 0:
                offset_from_zero[index, :] = offset[index, :].copy()
            else:
                offset_from_zero[index, :] = offset_from_zero[index - 1, :] + offset[index, :]

        # now convert to offset from center
        offset_from_center = offset_from_zero[int(round(self.frames.shape[0] / 2, 0)), :] - offset_from_zero
        # get max offset and determine size of helper frame
        max_offset = int(np.max(abs(offset_from_center)))
        size_frame = np.asarray(frame.shape, dtype=int)
        helper_size = tuple(size_frame + max_offset * 2)

        for frame_index, frame in enumerate(self.frames):
            # for each frame, make a helper frame (which is larger)
            bg = np.mean(frame)
            helper_image = np.ones(helper_size) * bg

            # and shift frame to correct position within the helper frame
            shift_dist = offset_from_center[frame_index, :].astype(int)
            helper_image[max_offset - shift_dist[-2]:size_frame[-2] + max_offset - shift_dist[-2],
                         max_offset - shift_dist[-1]:size_frame[-1] + max_offset - shift_dist[-1]] = frame

            # save to data_output
            data_output[frame_index, :, :] = helper_image[max_offset:size_frame[-2] + max_offset,
                                                          max_offset:size_frame[-1] + max_offset]
            # add to data merged by getting data_merged_helper
            helper_image = np.zeros(helper_size, dtype=self.data_type_signed)
            helper_image[max_offset - shift_dist[-2]:size_frame[-2] + max_offset - shift_dist[-2],
                         max_offset - shift_dist[-1]:size_frame[-1] + max_offset - shift_dist[-1]] = \
                data_merged_helper[frame_index, :, :]
            data_merged += helper_image[max_offset:size_frame[-2] + max_offset, max_offset:size_frame[-1] + max_offset]

            # if verbose, show result per frame
            if verbose:
                fig, ax = plt.subplots(1)
                ax.imshow(data_output[frame_index, :, :],
                          extent=[0, data_output[frame_index, :, :].shape[1],
                                  data_output[frame_index, :, :].shape[0],
                                  0], aspect='auto')
                plt.title("Frame {} shifted final".format(frame_index))
                plt.show()

        # if verbose, show overall result
        if verbose:
            fig, ax = plt.subplots(1)
            ax.imshow(data_merged, extent=[0, data_merged.shape[1], data_merged.shape[0], 0], aspect='auto')
            plt.title("Result")
            plt.show()

        return data_output, data_merged

    def find_energy_width(self):
        """
        Finds the bandwidth of the filters in eV. Assumed a 10 nm bandwidth, calculates the lower bound (lb)
        and upper bound (ub) and determines the bandwidth in eV.
        :return: energy_width: the bandwidth in eV for each filter.
        """
        energy_width = np.zeros(len(self.wavelengths))
        for index, wavelength in enumerate(self.wavelengths):
            wavelength_ev_lb = 1240 / (wavelength - 5)
            wavelength_ev_ub = 1240 / (wavelength + 5)
            energy_width[index] = wavelength_ev_lb - wavelength_ev_ub
        return energy_width

    # %% Run
    def run(self, verbose=False):
        """
        Main of HSM. Does all the work.
        ---------------------------------------
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
        energy_width = self.find_energy_width()
        # if verbose, show ROIs after correction
        if verbose:
            fig, ax = plt.subplots(1)
            figuring.plot_rois(ax, self.frame_for_rois, roi_locations=self.active_rois, roi_size=9,
                               roi_offset=self.roi_offset)
            plt.show()

        # prep for fitting
        roi_size_1d = 4
        fitter = HSMFit(roi_size_1d)

        # %% Fit every ROI for every frame
        for roi_index, roi in enumerate(self.active_rois):
            frame_stack = roi.get_frame_stack(self.corrected, roi_size_1d, self.roi_offset)
            # Fit with Gaussian fitter
            raw_intensity, intensity, raw_fits = fitter.fitter(frame_stack, shape, energy_width)

            # Fit the total intensity of a single ROI over all frames with Lorentzian
            if verbose:
                fig, ax = plt.subplots(1)
                ax.plot(self.wavelengths, intensity)
                ax.set_title('Result ROI #{}'.format(roi.index))
                plt.show()

            # fit lorentzian to individual fits
            result, r_squared = self.fit_lorentzian(intensity, self.wavelengths, verbose=verbose)

            result_dict = {"type": self.type, 'wavelengths': self.wavelengths, "lambda": 1240 / result[2],  # SP lambda
                           "linewidth": 1000 * result[3], 'R2': r_squared, "fit_parameters": result,  # linewidth
                           "raw_intensity": raw_intensity, "raw_fits": raw_fits,  # raw gaussian fits
                           "intensity": intensity, "raw": frame_stack}
            roi.results[self.name_result] = result_dict

            # progress update
            self.experiment.progress_updater.update_progress()

    def fit_lorentzian(self, scattering, wavelength, split=False, verbose=False):
        """
        Function to fit a lorentzian to the found intensities
        -----------------------------------------------
        :param scattering: the scattering intensities found
        :param wavelength: the wavelengths of the found intensities
        :param split: if True, this function has been called recursively with only a part of the original array
        :param verbose: if True, you get a lot of images
        :return: result: resulting Lorentzian parameters
        :return r_squared: the r-squared of this fit
        """

        def lorentzian(params, x):
            """
            Lorentzian formula. Taken from SPectrA
            ----------------
            :param params: Parameters of lorentzian. Need to be four.
            :param x: x-axis. Wavelengths
            :return: array of values for current parameters and wavelengths
            """
            return params[0] + params[1] / ((x - params[2]) ** 2 + (0.5 * params[3]) ** 2)

        def error_func(p, x, y):
            """
            Error function
            """
            return lorentzian(p, x) - y

        def find_r_squared(f, p, x, y):
            """
            Finds R^2 of a function and fitted result
            --------------------------
            :param f: function
            :param p: parameters
            :param x: x-axis
            :param y: true y-axis
            :return: R^2
            """
            res = y - f(p, x)
            ss_res = np.sum(res ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            return 1 - ss_res / ss_tot

        def compare_plot(x, y, p):
            """
            Automatically plot true x/y and fit
            ------------------------
            :param x: x-axis
            :param y: true y-axis
            :param p: parameters for fit
            :return: None. Shows fit
            """
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
        wavelength_ev = 1240 / wavelength

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
            result_full_std = least_squares(error_func, [-10, 100, 1240 / wavelength_ev[idx_max], 0.15],
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

        # if verbose, show comparison
        if verbose:
            compare_plot(wavelength_ev, scattering, result)

        return result, r_squared
