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

"""
# General
import os
import numpy as np

# I/O
from scipy.io import loadmat
import mat73

# A lot of scipy for signal processing
from scipy.ndimage import median_filter, shift, correlate
from scipy.optimize import curve_fit
import scipy.fft as fft
from scipy.signal import fftconvolve

# Own code
from _code.nd2_reading import ND2ReaderSelf
import _code.fitters as fitting


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

        self.frames, self.frame_merge = self.hsm_drift()

        # correct ROI locations for drift in HSM images

        #  corr_norm = self.normxcorr2_large(self.frame_merge, self.frame_zero)
        #  corr_norm_2 = self.normxcorr2_large(self.frame_zero, self.frame_merge)
        #  corr_norm2 = correlate2d(self.frame_zero, self.frame_merge)
        corr = fftconvolve(self.frame_zero, self.frame_merge)

        maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
        offset = maxima - np.asarray(self.frame_zero.shape)
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
                #  perr = np.sqrt(np.diag(cov))  # NOT WORKING

                self.hsm_result[roi_index, 0] = roi_index
                self.hsm_result[roi_index, 1] = res[1]
                self.hsm_result[roi_index, 2] = res[0]
                #  self.hsm_result[roi_index, 3] = perr ** 2

        intensity_result = np.concatenate((np.array(range(self.roi_locations.shape[0]))[:, np.newaxis], intensity),
                                          axis=1)

        return self.hsm_result, intensity_result

    # %% Correct for drift between frames
    def hsm_drift(self):

        frame = self.frames[0, :, :]
        data_output = np.zeros(self.frames.shape, dtype=frame.dtype)
        data_merged = np.zeros(frame.shape, dtype=frame.dtype)

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
                frame_convolution = self.normxcorr2(img_corrected[frame_index - 1, :, :],
                                                    img_corrected[frame_index, :, :])
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
            data_merged += np.asarray(np.round(np.abs(shifted_frame - bg), 0), dtype=frame.dtype)

        return data_output, data_merged

    @staticmethod
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

    @staticmethod
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

