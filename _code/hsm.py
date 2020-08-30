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

"""
import os
from _code.nd2_reading import ND2ReaderSelf
import numpy as np
from scipy.io import loadmat
import mat73
from scipy.signal import convolve
from scipy.ndimage import median_filter
from scipy.ndimage import shift
import _code.fitters as fitting
from scipy.optimize import curve_fit


class HSM:
    """
    Class for HSM of MBx Python
    """
    def __init__(self, directory, frame_zero, roi_locations, metadata, correction):

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
    def main(self):

        def find_nearest(match_array, value_array, match_index):
            idx = (np.abs(match_array - match_index)).argmin()
            return value_array[idx]

        # find correct shape for wavelength

        wavelength = np.asarray([int(nd2[:-4]) for nd2 in self.nd2_files], dtype=np.uint16)
        shape = np.asarray([find_nearest(self.spec_wavelength, self.spec_shape, nu) for nu in wavelength])

        # correct frames for drift

        self.frames, self.frame_merge = self.hsm_drift()

        # correct ROI locations for drift in HSM images

        corr = convolve(self.frame_merge, self.frame_zero)
        maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
        offset = maxima - np.asarray(self.frame_merge.shape)
        self.roi_locations += offset

        # prep for fitting

        raw_intensity = np.zeros((self.roi_locations.shape[0], self.frames.shape[0]))
        intensity = np.zeros((self.roi_locations.shape[0], self.frames.shape[0]))
        roi_size = 9
        roi_size_1d = int((roi_size - 1) / 2)
        fitter = fitting.Gaussian(roi_size, {}, "None", "Gaussian", 5, "Yes")

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
                    result, _, success = fitter.fit_gaussian(my_roi, roi_index)
                    raw_intensity[roi_index, frame_index] = 2 * np.pi * result[0] * result[3] * result[4]
                    intensity[roi_index, frame_index] = raw_intensity[roi_index, frame_index] / shape[frame_index]

                # %% Fit the total intensity of a single ROI over all frames with Lorentzian

                res, cov = curve_fit(lorentzian, wavelength, intensity[roi_index, :])
                perr = np.sqrt(np.diag(cov))  # NOT WORKING

                self.hsm_result[roi_index, 0] = roi_index
                self.hsm_result[roi_index, 1] = res[1]
                self.hsm_result[roi_index, 2] = res[0]
                self.hsm_result[roi_index, 3] = perr ** 2

        intensity_result = np.concatenate((np.array(range(self.roi_locations.shape[0]))[:, np.newaxis], intensity),
                                          axis=1)

        return self.hsm_result, intensity_result

    # %% Correct for drift between frames
    def hsm_drift(self):

        frame = self.frames[0, :, :]
        data_output = np.zeros(self.frames.shape)
        data_merged = np.zeros(frame.shape)

        offset = np.zeros((self.frames.shape[0], 2))
        offset_from_zero = np.zeros((self.frames.shape[0], 2))

        img_corrected = np.zeros((self.frames.shape[0], int(frame.shape[0] * 0.5 + 1), int(frame.shape[1] * 0.5 + 1)),
                                 dtype=np.int16)

        for frame_index, frame in enumerate(self.frames):
            data_output[frame_index, :, :] = frame.copy()
            background = median_filter(frame, size=9, mode='constant')
            frame = frame.astype(np.int16) - background
            img_corrected_single = np.round((frame[int(frame.shape[0] * 0.25 - 1):int(frame.shape[0] * 0.75),
                                             int(frame.shape[1] * 0.25 - 1):int(frame.shape[1] * 0.75)]),
                                            0).astype(np.int16)
            img_corrected[frame_index, :, :] = img_corrected_single

            if frame_index > 0:
                frame_convolution = convolve(img_corrected[frame_index - 1, :, :], img_corrected[frame_index, :, :])
                maxima = np.transpose(np.asarray(np.where(frame_convolution == np.amax(frame_convolution))))[0]
                offset[frame_index, :] = maxima - np.asarray(img_corrected_single.shape)

        for index in range(self.frames.shape[0]):
            if index == 0:
                offset_from_zero[index, :] = offset[index, :].copy()
            else:
                offset_from_zero[index, :] = offset_from_zero[index - 1, :] + offset[index, :]

        offset_from_center = offset_from_zero[int(round(self.frames.shape[0] / 2, 0)), :] - offset_from_zero

        for frame_index, frame in enumerate(self.frames):
            shift_dist = tuple(offset_from_center[frame_index, :].tolist())
            bg = np.mean(frame)
            shifted_frame = shift(frame, shift_dist, cval=bg)
            data_output[frame_index, :, :] = shifted_frame
            data_merged += np.round(np.abs(shifted_frame - bg), 0)

        return data_output, data_merged
