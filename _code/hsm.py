# -*- coding: utf-8 -*-
"""
Created on Tue 04/08/2020

@author: Dion Engels
MBx Python Data Analysis

hsm

This package is for the HSM part of MBx Python.

----------------------------

v0.0.1: Loading in multiple nd2, finding .mats

"""
import os
from _code.tools import ND2ReaderSelf
import numpy as np
from scipy.io import loadmat
import mat73
from scipy.signal import medfilt, fftconvolve
from scipy.ndimage import shift


class HSM:
    """
    Class for HSM of MBx Python
    """

    def __init__(self, directory, roi_locations, metadata, correction):

        self.roi_locations = roi_locations
        self.metadata = metadata

        # load in frame

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

    def main(self):

        def find_nearest(match_array, value_array, match_index):
            idx = (np.abs(match_array - match_index)).argmin()
            return value_array[idx]

        wavelength = np.asarray([int(nd2[:-4]) for nd2 in self.nd2_files], dtype=np.uint16)
        shape = np.asarray([find_nearest(self.spec_wavelength, self.spec_shape, nu) for nu in wavelength])

        self.frames, self.frame_merge = self.hsm_drift()

        return

    @staticmethod
    def normxcorr2(template, image, mode="full"):

        a1 = np.ones(template.shape)
        ar = np.flipud(np.fliplr(template))

        out = fftconvolve(image, ar.conj(), mode=mode)
        image = fftconvolve(np.square(image), a1, mode=mode) - \
                np.square(fftconvolve(image, a1, mode=mode)) / (np.prod(template.shape))

        # Remove small machine precision errors after subtraction
        image[np.where(image < 0)] = 0
        template = np.sum(np.square(template))
        out = out / np.sqrt(image * template)

        out[np.where(np.logical_not(np.isfinite(out)))] = 0
        return out

    def hsm_drift(self):

        frame = self.frames[0, :, :]

        data_output = np.zeros(self.frames.shape)
        data_merged = np.zeros(frame.shape)
        offset = np.zeros((self.frames.shape[0], 2))
        offset_from_zero = np.zeros((self.frames.shape[0], 2))
        img_corrected = np.zeros((self.frames.shape[0], int(frame.shape[0] * 0.5 + 1), int(frame.shape[1] * 0.5 + 1)),
                                 dtype=np.int16)
        img_corrected2 = np.zeros((self.frames.shape[0], int(frame.shape[0] * 0.5 + 1), int(frame.shape[1] * 0.5 + 1)),
                                  dtype=np.int16)

        for frame_index, frame in enumerate(self.frames):
            bg = np.mean(frame)
            data_output[frame_index, :, :] = frame.copy()
            img_corrected_single = np.round((frame[int(frame.shape[0] * 0.25 - 1):int(frame.shape[0] * 0.75),
                                             int(frame.shape[1] * 0.25 - 1):int(frame.shape[1] * 0.75)] - bg),
                                            0).astype(np.int16)

            background = medfilt(img_corrected_single, kernel_size=9)
            background[background == 0] = np.min(background[background > 0])
            img_corrected2[frame_index, :, :] = img_corrected_single - background
            img_corrected[frame_index, :, :] = img_corrected_single

            if frame_index > 0:
                frame_convolution = self.normxcorr2(img_corrected2[frame_index - 1, :, :],
                                                    img_corrected2[frame_index, :, :])
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
