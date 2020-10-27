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
from skimage.feature import match_template
import numpy as np
from pims.process import crop

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

    def get_frame_stack(self, frames, roi_size_1d, offset, shape):
        """
        Gets ROI for a certain frame stack, offset, and ROI size
        ------------------------------------
        :param frames: frames to get ROI of
        :param roi_size_1d: ROI size
        :param offset: offset of current ROI in that frame
        :param shape: shape of the frame
        :return: Something by Something ROI around the x/y position of this ROI, in time
        """
        return crop(frames, ((self.y + offset[0] - roi_size_1d, shape[0] + offset[0] - self.y - roi_size_1d - 1),
                             (self.x + offset[1] - roi_size_1d, shape[1] + offset[1] - self.x - roi_size_1d - 1)))

    def get_frame_stack_np(self, frames, roi_size_1d, offset):
        """
        Gets ROI for a certain frame stack, offset, and ROI size. This is for numpy implementation
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
    def __init__(self, experiment, nd2, name):
        """
        Init for dataset class. Sets name, type, name_result (for MATLAB) and some other base things to None
        -----------------------------------
        :param experiment: parent experiment
        :param name: name of dataset
        """
        self.type = "Dataset"
        self.experiment = experiment
        self.data_type = nd2.pixel_type
        bits = int("".join([s for s in str(self.data_type) if s.isdigit()]))  # complicated way to get #bits
        if bits < 16:
            self.data_type_signed = np.int16
        elif bits < 32:
            self.data_type_signed = np.int32
        else:
            self.data_type_signed = np.int64
        self.name = name.split(".")[0].split("/")[-1]
        self.filename = name
        self.name_result = self.set_result_name(self.name)
        self.frames = None
        self.frame_for_rois = None
        self.metadata = None
        self.fitter = None
        self.roi_offset = None
        self.settings = None
        self.active_rois = []

    @staticmethod
    def check_name_validity(new_name):
        """
        Check if name for a dataset is conform MATLAB requirements
        ---------------------------
        :param new_name: string of new name
        :return: return if valid or not
        """
        chars = set(new_name)

        for i in range(0, 10):
            chars.discard(str(i))

        # loop over all letters
        from string import ascii_letters
        for char in ascii_letters:
            chars.discard(char)

        chars.discard('-')
        chars.discard('_')
        chars.discard(' ')

        if len(chars) > 0:
            return False
        else:
            return True

    @staticmethod
    def set_result_name(new_name):
        """
        Set name for result struct
        --------------------
        :param new_name: new name to adapt to self.name
        :return:
        """
        # check matlab rules
        tmp_name = new_name.replace(' ', '_')  # no spaces
        tmp_name = tmp_name.replace('-', '_')  # no -
        tmp_name = tmp_name[-59:]  # take only last 59 characters
        return "res_" + tmp_name

    def set_name(self, new_name):
        """
        Set a name
        -----------------
        :param new_name: new name
        :return: None. Edits class
        """
        self.name = new_name
        self.name_result = self.set_result_name(new_name)

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
            return slice(0, None), 0
        elif start == "Leave empty for start" and end != "Leave empty for end":
            return slice(0, int(end)), 0
        elif start != "Leave empty for start" and end == "Leave empty for end":
            return slice(int(start), None), int(start)
        else:  # start != "Leave empty for start" and end != "Leave empty for end":
            return slice(int(start), int(end)), int(start)

    @staticmethod
    def correlate_frames_other_size(frame_small, frame_big):
        """
        Correlates other sizes small and big frame and finds offset between the two
        -----------------------
        :param frame_small: Previous frame
        :param frame_big: New frame to correlate with
        :return: offset: the offset between the two frames
        """
        corr = match_template(frame_big, frame_small)
        offset = -np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
        return offset

    @staticmethod
    def correlate_frames_same_size(frame_old, frame_new):
        """
        Correlates same sizes old frame and new frame and finds offset between the two
        -----------------------
        :param frame_old: Previous frame
        :param frame_new: New frame to correlate with
        :return: offset: the offset between the two frames
        """
        # if same matrix
        if np.array_equal(frame_new, frame_old):
            offset = np.asarray([0, 0])
        else:
            corr = normxcorr2(frame_old, frame_new)
            maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
            offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])

        return offset

    @staticmethod
    def check_slice_validity(total_frame, x_slice, y_slice):
        if x_slice.stop is None:
            x_slice = slice(x_slice.start, total_frame.shape[1])
        if y_slice.stop is None:
            y_slice = slice(y_slice.start, total_frame.shape[0])

        if x_slice.start > x_slice.stop or x_slice.start < 0 or x_slice.stop > total_frame.shape[1]:
            return False
        if y_slice.start > y_slice.stop or y_slice.start < 0 or y_slice.stop > total_frame.shape[0]:
            return False

        return True

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
        if frame_shape[0] > experiment_frame_shape[0] or frame_shape[1] > experiment_frame_shape[1]:
            small_frame = self.experiment.frame_for_rois
            cropped_frame = self.frame_for_rois[y_slice, x_slice]

            # if slices are present, check validity
            if self.check_slice_validity(self.frame_for_rois, x_slice, y_slice) is False:
                self.experiment.error_func("Slice invalid", "Slice size is not valid")
                return None

            # check if sliced frame is valid
            if small_frame.shape[0] > cropped_frame.shape[0] or small_frame.shape[1] > cropped_frame.shape[1]:
                self.experiment.error_func("Crop too tight", "Cropped frame now smaller than previously smaller frame")
                return None
            offset = -self.correlate_frames_other_size(small_frame, cropped_frame) + offset_crop
        elif frame_shape[0] < experiment_frame_shape[0] or frame_shape[1] < experiment_frame_shape[1]:
            # if other way around
            small_frame = self.frame_for_rois
            cropped_frame = self.experiment.frame_for_rois[y_slice, x_slice]

            # if slices are present, check validity
            if self.check_slice_validity(self.experiment.frame_for_rois, x_slice, y_slice) is False:
                self.experiment.error_func("Slice invalid", "Slice size is not valid")
                return None

            # check if sliced frame is valid
            if small_frame.shape[0] > cropped_frame.shape[0] or small_frame.shape[1] > cropped_frame.shape[1]:
                self.experiment.error_func("Crop too tight", "Cropped frame now smaller than previously smaller frame")
                return None
            offset = self.correlate_frames_other_size(small_frame, cropped_frame) - offset_crop
        else:
            # if same size
            old_frame = self.experiment.frame_for_rois
            new_frame = self.frame_for_rois
            offset = self.correlate_frames_same_size(old_frame, new_frame)
        return offset

    def find_rois(self, settings):
        """
        Finds the ROIs within a new dataset. First correlates frames, then uses roi_in_frame func
        --------------------
        :param settings: settings from user
        :return: None. Edits dataset class by self.active_rois and self.roi_offset
        """
        self.roi_offset = self.correlate(settings)
        if self.roi_offset is None:
            self.active_rois = []
            return
        self.active_rois = [roi for roi in self.experiment.rois if roi.in_frame(self.frame_for_rois.shape,
                                                                                self.roi_offset,
                                                                                self.experiment.
                                                                                roi_finder.side_distance)]

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
    b = int(sum(b.flatten().astype(np.int64) ** 2))
    c = c / np.sqrt(a * b)
    return c
