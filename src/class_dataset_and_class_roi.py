# -*- coding: utf-8 -*-
"""
Created on Thu 01/10/2020

----------------------------

@author: Dion Engels
MBx Python Data Analysis

class dataset & roi

The dataset and ROI class of v2 of program. Dataset is one nd2 file, ROIs are region of interest.

-----------------

v2.0: part of v2.0: 15/10/2020

"""
# GENERAL IMPORTS
import scipy.fft as fft
from scipy.signal import fftconvolve
import numpy as np

__self_made__ = True

# %% ROI


class Roi:
    """
    ROI class. Used to determine region of interest
    """
    def __init__(self, x, y):
        """
        Initialization of ROI class.
        ---------------------------
        :param x: x position
        :param y: y position
        """
        self.x = x
        self.y = y
        self.index = None

        self.results = {}

    def set_index(self, index):
        """
        Sets index of ROI
        ----------------
        :param index: index to set
        """
        self.index = index

    def get_roi(self, frame, roi_size_1d, offset):
        """
        Gets ROI for a certain frame, offset, and ROI size
        ------------------------------------
        :param frame: frame to get ROI of
        :param roi_size_1d: ROI size
        :param offset: offset of current ROI in that frame
        :return: Something by Something ROI around the x/y position of this ROI
        """
        return frame[self.y + offset[0] - roi_size_1d:self.y + offset[0] + roi_size_1d + 1,
                     self.x + offset[1] - roi_size_1d:self.x + offset[1] + roi_size_1d + 1]

    def get_frame_stack(self, frames, roi_size_1d, offset):
        """
        Gets ROI for a certain frame stack, offset, and ROI size
        ------------------------------------
        :param frames: frames to get ROI of
        :param roi_size_1d: ROI size
        :param offset: offset of current ROI in that frame
        :return: Something by Something ROI around the x/y position of this ROI, in time
        """
        return frames[:, self.y + offset[0] - roi_size_1d:self.y + offset[0] + roi_size_1d + 1,
                      self.x + offset[1] - roi_size_1d:self.x + offset[1] + roi_size_1d + 1]

    def in_frame(self, shape, offset, margin):
        """
        Checks whether or not this ROI is in the frame
        --------------------------
        :param shape: Shape of frame
        :param offset: offset of frame
        :param margin: margin required to edge
        :return: in_frame_boolean: whether or not in frame
        """
        if self.x + offset[1] < margin or self.x + offset[1] > shape[1] - margin:
            in_frame_boolean = False
        elif self.y + offset[0] < margin or self.y + offset[0] > shape[0] - margin:
            in_frame_boolean = False
        else:
            in_frame_boolean = True

        return in_frame_boolean

# %% Dataset


class Dataset:
    """
    Base dataset class. Each dataset type (HSM / TT) inherits from this
    """
    def __init__(self, experiment, name):
        """
        Init for dataset class. Sets name, type, name_result (for MATLAB) and some other base things to None
        -----------------------------------
        :param experiment: parent experiment
        :param name: name of dataset
        """
        self.type = "Dataset"
        self.experiment = experiment
        self.name = name.split(".")[0].split("/")[-1]
        self.filename = name
        self.name_result = 'res_{}'.format(self.name.replace(' ', '_'))
        self.frames = None
        self.frame_for_rois = None
        self.metadata = None
        self.fitter = None
        self.roi_offset = None
        self.settings = None
        self.active_rois = []

    @staticmethod
    def parse_start_end(start, end):
        """
        Parses start and end values often used in program to create a slice
        ---------------------------
        :param start: Start value
        :param end: End value
        :return: Slice from start to end value
        """
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
        """
        Correlates old frame with new frame and finds offset between the two
        -----------------------
        :param frame_old: Previous frame
        :param frame_new: New frame to correlate with
        :return: offset: the offset between the two frames
        """
        # if same matrix
        if np.array_equal(frame_new, frame_old):
            return [0, 0]

        # if same size
        if frame_old.shape == frame_new.shape:
            corr = normxcorr2(frame_old, frame_new)
            maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
            offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
        else:
            # if not same size
            corr = normxcorr2_large(frame_old, frame_new)
            maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
            offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
        return offset

    def correlate(self, settings):
        """
        The overall correlate function. Calls the above functions
        ----------------------------------
        :param settings: settings from user
        :return: offset: offset between new and old frame
        """
        # finds slices and offset
        x_slice, x_offset = self.parse_start_end(settings['x_min'], settings['x_max'])
        y_slice, y_offset = self.parse_start_end(settings['y_min'], settings['y_max'])
        offset_crop = np.asarray([y_offset, x_offset])

        experiment_frame_shape = self.experiment.frame_for_rois.shape
        frame_shape = self.frame_for_rois.shape

        # if frame is larger than experiment frame
        if frame_shape[0] > experiment_frame_shape[0] and frame_shape[1] > experiment_frame_shape[1]:
            small_frame = self.experiment.frame_for_rois
            cropped_frame = self.frame_for_rois[y_slice, x_slice]
            offset = self.correlate_frames(small_frame, cropped_frame) - offset_crop
        elif frame_shape[0] < experiment_frame_shape[0] and frame_shape[1] < experiment_frame_shape[1]:
            # if other way around
            small_frame = self.frame_for_rois
            cropped_frame = self.experiment.frame_for_rois[y_slice, x_slice]
            offset = self.correlate_frames(cropped_frame, small_frame) + offset_crop
        else:
            # if same size
            old_frame = self.experiment.frame_for_rois
            new_frame = self.frame_for_rois
            offset = self.correlate_frames(old_frame, new_frame)
        return offset

    def find_rois(self, settings):
        """
        Finds the ROIs within a new dataset. First correlates frames, then uses roi_in_frame func
        --------------------
        :param settings: settings from user
        :return: None. Edits dataset class by self.active_rois and self.roi_offset
        """
        self.roi_offset = self.correlate(settings)
        self.active_rois = [roi for roi in self.experiment.rois if roi.in_frame(self.frame_for_rois.shape,
                                                                                self.roi_offset,
                                                                                self.experiment.
                                                                                roi_finder.side_distance)]

    def set_name(self, new_name):
        """
        Set a name
        -----------------
        :param new_name: new name
        :return: None. Edits class
        """
        self.name = new_name
        self.name_result = 'res_{}'.format(self.name.replace(' ', '_'))

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
