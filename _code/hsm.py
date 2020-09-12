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

"""
# General
import os
import numpy as np

# I/O
from scipy.io import loadmat
import mat73

# A lot of scipy for signal processing
from scipy.ndimage import median_filter, shift
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import scipy.fft as fft
from scipy.signal import fftconvolve

# Own code
from _code.nd2_reading import ND2ReaderSelf
import _code.fitters as fitting
import _code.figure_making as figuring

import matplotlib.pyplot as plt


class HSM:
    """
    Class for HSM of MBx Python
    """

    def __init__(self, directory, frame_zero, roi_locations, metadata, correction):

        def only_ints(list_to_check):
            new_list = []
            for string in list_to_check:
                try:
                    int(string[:-4])
                    new_list.append(string)
                except ValueError:
                    pass

            return new_list

        self.roi_locations = roi_locations
        self.hsm_result = np.zeros((roi_locations.shape[0], 4))
        self.metadata = metadata
        self.frame_zero = frame_zero
        self.frame_merge = None

        # load in nd2 frames

        self.nd2_dir = directory[0]
        nd2_dir_files = os.listdir(self.nd2_dir)
        self.nd2_files = [nd2 for nd2 in nd2_dir_files if nd2.endswith('.nd2')]

        if 'Merged Documents.nd2' in self.nd2_files:
            self.merged = [nd2 for nd2 in self.nd2_files if nd2.__contains__('Merged')]
            self.nd2_files = [nd2 for nd2 in self.nd2_files if not nd2.__contains__('Merged')]
            self.merged = ''.join(self.merged)
            nd2 = ND2ReaderSelf(self.nd2_dir + "/" + self.merged)
            self.frames = np.asarray(nd2, dtype=np.uint16).copy()
            nd2.close()
        else:
            self.nd2_files = only_ints(self.nd2_files)
            for nd2_index, nd2_name in enumerate(self.nd2_files):
                nd2 = ND2ReaderSelf(self.nd2_dir + "/" + nd2_name)
                frame = np.asarray(nd2, dtype=np.uint16).copy()
                if nd2_index == 0:
                    self.frames = np.zeros((len(self.nd2_files), frame.shape[1], frame.shape[2]), dtype=np.uint16)
                self.frames[nd2_index, :, :] = frame
                nd2.close()

        # load in spectral correction

        path = os.getcwd()
        path += ("/spectral_corrections/" + correction + ".mat")
        self.correction_path = path
        try:  # for new MATLAB versions
            correction = loadmat(path)
            self.spec_wavelength = correction['SpectralCorrection'][0][0][0][0]
            self.spec_shape = correction['SpectralCorrection'][0][0][1][0]
        except NotImplementedError:  # for old MATLAB versions
            correction = mat73.loadmat(path)
            self.spec_wavelength = correction['SpectralCorrection']['Lambda']
            self.spec_shape = correction['SpectralCorrection']['SpecShape']

    # %% Main
    def main(self, verbose=False):

        def find_nearest(match_array, value_array, match_index):
            idx = (np.abs(match_array - match_index)).argmin()
            return value_array[idx]

        def to_int(list_to_int):
            new_list = []
            for string in list_to_int:
                try:
                    int(string[:-4])
                    new_list.append(int(string[:-4]))
                except ValueError:
                    pass

            new_list = np.asarray(new_list, dtype=np.uint16)
            return new_list

        # find correct shape for wavelength

        wavelength = to_int(self.nd2_files)
        shape = np.asarray([find_nearest(self.spec_wavelength, self.spec_shape, nu) for nu in wavelength])

        # correct frames for drift

        self.frames, self.frame_merge = self.hsm_drift(verbose=False)

        # correct ROI locations for drift in HSM images

        size_hsm_frame = np.asarray(self.frame_merge.shape, dtype=int)
        size_laser_frame = np.asarray(self.frame_zero.shape, dtype=int)

        corr = normxcorr2_large(self.frame_zero, self.frame_merge)
        maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
        offset = maxima - np.asarray(self.frame_zero.shape) + np.asarray([1, 1])
        self.roi_locations += offset

        if verbose:
            fig, ax = plt.subplots(1)
            ax.imshow(self.frame_merge, extent=[0, self.frame_merge.shape[1], self.frame_merge.shape[0], 0],
                      aspect='auto')
            plt.title("Merged frame")
            plt.show()

            fig, ax = plt.subplots(1)
            ax.imshow(self.frame_zero, extent=[0, self.frame_zero.shape[1], self.frame_zero.shape[0], 0],
                      aspect='auto')
            plt.title("Frame zero")
            plt.show()

            fig, ax = plt.subplots(1)
            frame_670 = np.where(wavelength == 670)[0][0]
            cutout_670 = self.frames[frame_670, offset[0]:offset[0] + size_laser_frame[0],
                         offset[1]:offset[1] + size_laser_frame[1]]
            ax.imshow(cutout_670, extent=[0, cutout_670.shape[1], cutout_670.shape[0], 0], aspect='auto')
            plt.title("Frame 670 nm")
            plt.show()

            figuring.plot_rois(self.frame_merge, self.roi_locations, 9)

            figuring.plot_rois(self.frames[frame_670, :, :], self.roi_locations, 9)

        # prep for fitting

        raw_intensity = np.zeros((self.roi_locations.shape[0], self.frames.shape[0]))
        intensity = np.zeros((self.roi_locations.shape[0], self.frames.shape[0]))
        roi_size = 9
        roi_size_1d = int((roi_size - 1) / 2)
        fitter = fitting.Gaussian(roi_size, {}, "None", "Gaussian", 5, 500)

        pos_max = roi_size
        pos_min = 0
        sig_min = 0
        sig_max = roi_size_1d + 1

        def lorentzian(x, a, x0):
            return a / ((x - x0) ** 2 + a ** 2) / np.pi

        # %% Fit every ROI for every frame

        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1])

            if x < roi_size_1d or y < roi_size_1d or \
                    x > self.frames.shape[1] - roi_size_1d or y > self.frames.shape[2] - roi_size_1d:
                pass
            else:

                for frame_index, frame in enumerate(self.frames):
                    my_roi = frame[y - roi_size_1d:y + roi_size_1d + 1, x - roi_size_1d:x + roi_size_1d + 1]
                    if verbose:
                        fig, ax = plt.subplots(1)
                        ax.imshow(my_roi, extent=[0, my_roi.shape[1], my_roi.shape[0], 0], aspect='auto')
                        ax.set_xlabel('x (pixels)')
                        ax.set_ylabel('y (pixels)')
                        ax.set_title('Zoom-in ROI #{}, frame #{}'.format(roi_index, frame_index))
                        plt.show()
                    result, _, success = fitter.fit_gaussian(my_roi, roi_index)
                    if success == 0 or \
                            result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max or \
                            result[3] < sig_min or result[3] > sig_max or result[4] < sig_min or result[4] > sig_max:
                        raw_intensity[roi_index, frame_index] = np.nan
                        intensity[roi_index, frame_index] = np.nan
                    else:
                        raw_intensity[roi_index, frame_index] = 2 * np.pi * result[0] * result[3] * result[4]
                        intensity[roi_index, frame_index] = raw_intensity[roi_index, frame_index] / shape[frame_index]

                # %% Fit the total intensity of a single ROI over all frames with Lorentzian

                if verbose:
                    fig, ax = plt.subplots(1)
                    ax.plot(wavelength, raw_intensity[roi_index, :])
                    ax.set_title('Result ROI #{}'.format(roi_index))
                    plt.show()

                #res, cov = curve_fit(lorentzian, wavelength, intensity[roi_index, :])
                result, fit, r_squared = self.fit_lorentzian(intensity[roi_index, :], wavelength)
                #  perr = np.sqrt(np.diag(cov))  # NOT WORKING

                self.hsm_result[roi_index, 0] = roi_index
                self.hsm_result[roi_index, 1] = res[1]
                self.hsm_result[roi_index, 2] = res[0]
                #  self.hsm_result[roi_index, 3] = perr ** 2

        intensity_result = np.concatenate((np.array(range(self.roi_locations.shape[0]))[:, np.newaxis], intensity),
                                          axis=1)

        return self.hsm_result, intensity_result

    # %% Correct for drift between frames
    def hsm_drift(self, verbose=False):

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
            data_merged += helper_image[max_offset:size_frame[-2] + max_offset,
                                        max_offset:size_frame[-1] + max_offset]

            if verbose:
                fig, ax = plt.subplots(1)
                ax.imshow(data_output[frame_index, :, :],
                          extent=[0, data_output[frame_index, :, :].shape[1], data_output[frame_index, :, :].shape[0],
                                  0], aspect='auto')
                plt.title("Frame {} shifted final".format(frame_index))
                plt.show()

        if verbose:
            fig, ax = plt.subplots(1)
            ax.imshow(data_merged, extent=[0, data_merged.shape[1], data_merged.shape[0], 0], aspect='auto')
            plt.title("Result")
            plt.show()

        return data_output, data_merged

    @staticmethod
    def fit_lorentzian(scattering, wavelength):
        def lorentzian(p, x):
            return p[0] + p[1] / ((x - p[2])**2 + (0.5 * p[3])**2)

        def find_r_squared(f, p, x, y):
            res = y - f(p, x)
            ss_res = np.sum(res**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            return 1 - ss_res / ss_tot

        # remove nans

        bool = ~np.isnan(scattering)
        scattering = scattering[bool]
        wavelength = wavelength[bool]

        # prep

        wavelength_ev = 1248 / wavelength
        lorentz_ev = np.arange(1248 / np.min(wavelength), 1248 / np.max(wavelength),
                               (1248 / np.min(wavelength - 1248/np.max(wavelength))) / (len(wavelength) * 7))
        lorentz_ev = np.append(lorentz_ev, 1248 / np.max(wavelength))
        lorentz_w1 = 1248 / lorentz_ev
        fit = np.zeros((len(lorentz_w1), 2))

        # find max and min

        max_sca = np.max(scattering[scattering < np.max(scattering)])
        idx_max = np.argmax(scattering[scattering < np.max(scattering)])
        min_sca = np.min(scattering)

        # init guess and first fit

        init_1w = abs(2 / (np.pi * max_sca) * np.trapz(scattering, wavelength_ev))
        init_guess = [min_sca, min_sca * init_1w / (2 * np.pi), wavelength_ev[idx_max], init_1w]
        result, cov_x, dict, mesg, ier = leastsq(lorentzian, init_guess, args=(wavelength_ev, scattering))
        result[3] = abs(result[3])

        # if fit is bad, retry

        if result[3] < 0.02:
            lb = abs(idx_max - round(len(scattering)/4, 0)) - 1  # -1 for MATLAB
            hb = abs(round(3*len(scattering) /4, 0)) - 1  # -1 for MATLAB
            init_guess = [min_sca, min_sca * init_1w / (2 * np.pi), wavelength_ev[idx_max], init_1w]
            result, cov_x, dict, mesg, ier = leastsq(lorentzian, init_guess,
                                                     args=(wavelength_ev[lb:hb], scattering[lb:hb]))
            result[4] = abs(result[4])

        # calc r-sqaured. If bad, retry

        r_squared = find_r_squared(lambda p, x: lorentzian(p, x), result, wavelength, scattering)

        if r_squared < 0.9:
            result, cov_x, dict, mesg, ier = leastsq(lorentzian, [-10, 100, 1248 / wavelength[idx_max], 0.15],
                                                     args=(wavelength_ev, scattering))
            result[4] = abs(result[4])
            r_squared = find_r_squared(lambda p, x: lorentzian(p, x), result, wavelength, scattering)

        # return

        fit[:, 0] = lorentz_w1
        fit[:, 1] = lorentzian(result, lorentz_ev)

        return result, fit, r_squared


def normxcorr2(b, a):
    def conv2(a, b):
        ma, na = a.shape
        mb, nb = b.shape
        return fft.ifft2(fft.fft2(a, [2 * ma - 1, 2 * na - 1]) * fft.fft2(b, [2 * mb - 1, 2 * nb - 1]))

    c = conv2(a, np.flipud(np.fliplr(b)))
    a = conv2(a ** 2, np.ones(b.shape))
    b = sum(b.flatten() ** 2)
    c = c / np.sqrt(a * b)
    return c


def normxcorr2_large(template, image, mode="full"):
    """
    Input arrays should be floating point numbers.
    :param template: N-D array, of template or filter you are using for cross-correlation.
    Must be less or equal dimensions to image.
    Length of each dimension must be less than length of image.
    :param image: N-D array
    :param mode: Options, "full", "valid", "same"
    full (Default): The output of fftconvolve is the full discrete linear convolution of the inputs.
    Output size will be image size + 1/2 template size in each dimension.
    valid: The output consists only of those elements that do not rely on the zero-padding.
    same: The output is the same size as image, centered with respect to the ‘full’ output.
    :return: N-D array of same dimensions as image. Size depends on mode parameter.
    """

    # If this happens, it is probably a mistake
    if np.ndim(template) > np.ndim(image) or \
            len([i for i in range(np.ndim(template)) if template.shape[i] > image.shape[i]]) > 0:
        print("normxcorr2: TEMPLATE larger than IMG. Arguments may be swapped.")

    template = template - np.mean(template)
    image = image - np.mean(image)

    a1 = np.ones(template.shape)
    # Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
    ar = np.flipud(np.fliplr(template))
    out = fftconvolve(image, ar.conj(), mode=mode)

    image = fftconvolve(np.square(image), a1, mode=mode) - \
            np.square(fftconvolve(image, a1, mode=mode)) / (np.prod(template.shape))

    # Remove small machine precision errors after subtraction
    image[np.where(image < 0)] = 0

    template = np.sum(np.square(template))
    out = out / np.sqrt(image * template)

    # Remove any divisions by 0 or very close to 0
    out[np.where(np.logical_not(np.isfinite(out)))] = 0

    return out
