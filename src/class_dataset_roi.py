# -*- coding: utf-8 -*-
"""
Created on Thu 01/10/2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

class dataset & roi

The dataset and ROI class of v2 of program. Dataset is one nd2 file, ROIs are region of interest.
"""
# GENERAL IMPORTS
import scipy.fft as fft
from scipy.signal import fftconvolve
import numpy as np

__self_made__ = True

# %% ROI


class Roi:

    def __init__(self, x, y):

        self.x = x
        self.y = y
        self.index = None

        self.results = {}

    def set_index(self, index):
        self.index = index

    def get_roi(self, frame, roi_size_1d):
        return frame[self.y - roi_size_1d:self.y + roi_size_1d + 1,
                     self.x - roi_size_1d:self.x + roi_size_1d + 1]

    def get_frame_stack(self, frames, roi_size_1d):
        return frames[:, self.y - roi_size_1d:self.y + roi_size_1d + 1,
                      self.x - roi_size_1d:self.x + roi_size_1d + 1]

    def in_frame(self, shape, offset):
        if self.x + offset[1] < 0 or self.x + offset[1] > shape[1]:
            in_frame_boolean = False
        elif self.y + offset[0] < 0 or self.y + offset[0] > shape[0]:
            in_frame_boolean = False
        else:
            in_frame_boolean = True

        return in_frame_boolean
# %% Dataset


class Dataset:
    def __init__(self, experiment, name):
        self.type = "Dataset"
        self.experiment = experiment
        self.name = name
        self.frames = None
        self.frame_for_rois = None
        self.metadata = None
        self.fitter = None
        self.drift_corrector = None
        self.roi_offset = None
        self.settings = None
        self.active_rois = []

    @staticmethod
    def parse_start_end(start, end):
        if start == "Leave empty for start" and end == "Leave empty for end":
            return slice(None), 0
        elif start == "Leave empty for start" and end != "Leave empty for end":
            return slice(0, int(end)), 0
        elif start != "Leave empty for start" and end == "Leave empty for end":
            return slice(int(start), None), start
        else:  # start != "Leave empty for start" and end != "Leave empty for end":
            return slice(int(start), int(end)), start

    @staticmethod
    def correlate_frames(frame_old, frame_new):
            if frame_old.shape == frame_new.shape:
                corr = normxcorr2(frame_old, frame_new)
                maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
                offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
            else:
                corr = normxcorr2_large(frame_old, frame_new)
                maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
                offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
            return offset

    def correlate(self, settings):
        x_slice, x_offset = self.parse_start_end(settings['x_min'], settings['x_max'])
        y_slice, y_offset = self.parse_start_end(settings['y_min'], settings['y_max'])
        offset_crop = np.asarray([y_offset, x_offset])

        experiment_frame_shape = self.experiment.frame_for_rois.shape
        frame_shape = self.frame_for_rois.shape

        # test offset crop
        if frame_shape[0] > experiment_frame_shape[0] and frame_shape[1] > experiment_frame_shape[1]:
            small_frame = self.experiment.frame_for_rois
            cropped_frame = self.frame_for_rois(y_slice, x_slice)
            offset = self.correlate_frames(small_frame, cropped_frame) - offset_crop
        elif frame_shape[0] < experiment_frame_shape[0] and frame_shape[1] < experiment_frame_shape[1]:
            small_frame = self.frame_for_rois
            cropped_frame = self.experiment.frame_for_rois(y_slice, x_slice)
            offset = self.correlate_frames(cropped_frame, small_frame) + offset_crop
        else:
            old_frame = self.experiment.frame_for_rois
            new_frame = self.frame_for_rois
            offset = self.correlate_frames(old_frame, new_frame)
        return offset

    def find_rois(self, settings, frame_for_rois, created_by):
        self.roi_offset = self.correlate(settings)
        self.active_rois = [roi for roi in self.experiment.rois if roi.in_frame(self.frame_for_rois.shape,
                                                                                self.roi_offset)]

# %% correlation


def normxcorr2(b, a):
    """
    Correlation of similar size frames
    """
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
    Correlation of frames with different sizes

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
