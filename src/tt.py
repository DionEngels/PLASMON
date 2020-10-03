# -*- coding: utf-8 -*-
"""
Created on Sun May 31 22:58:08 2020

@author: Dion Engels
MBx Python Data Analysis

Fitters

This package holds all the fitting algorithms of MBx Python.
This includes 2 Gaussian fitters and 3 Phasor fitters.

----------------------------

v0.0.1: Setup: 31/05/2020
v0.1: rainSTORM slow, but working: 04/06/2020
v0.2: optimizations after initial speed check: 17/06/2020
v0.3: optimization after second speed check: 27/06/2020
v0.4: removed all unneeded fitters
v0.5: Three versions of cached fitter
v0.6: FORTRAN enabled without cache, ST
v0.7: cleanup and declare_functions function: 14/07/2020
v0.8: PhasorSum, remove LastFit, rename, and clean up
v0.9: Gaussian fitting (with or without background) and three Phasor fitters: 24/07/2020
v0.9.1: NaN output instead of leaving out, MATLAB ready output
v0.9.2: removed minimum pixel requirement
v1.0: Integrated intensity output for Gaussians: 27/08/2020

"""
# %% Imports
from __future__ import division, print_function, absolute_import

import numpy as np  # general mathematics + arrays
from pims import Frame  # for converting ND2 generator to one numpy array

from scipy.optimize import _minpack, OptimizeResult  # for Gaussian fitter
from scipy.ndimage import median_filter  # for correlation with experiment

import src.mbx_fortran as fortran_linalg  # for fast self-made operations for Gaussian fitter
import src.mbx_fortran_tools as fortran_tools  # for fast self-made general operations
from src.class_dataset_roi import Dataset  # base dataset

from pyfftw import empty_aligned, FFTW  # for FFT for Phasor
from math import pi, atan2  # general mathematics
from cmath import phase  # general mathematics

__self_made__ = True

# %% Time trace class


class TimeTrace(Dataset):
    def __init__(self, experiment, nd2, name):
        super().__init__(experiment, name)
        self.frames = nd2
        background = median_filter(np.asarray(nd2[0]), size=9)
        self.frame_for_rois = np.asarray(nd2[0]).astype('float') - background
        self.metadata = nd2.get_metadata()
        self.pixels_or_nm = None
        self.slice = None
        self.n_cores = 1

    def prepare_run(self, settings):
        check = self.experiment.proceed_question("OK", "Cancel", "Are you sure?",
                                                 "Fitting may take a while. Are you sure everything is set up correctly?")
        if not check:
            return

        if settings['#cores'] > 1 and "Phasor" in settings['method']:
            cores_check = self.experiment.proceed_question("ok_cancel", "Just a heads up",
                                                           """Phasor will be used with one core since the
            overhead only slows it down""")
            if not cores_check:
                return
            settings['#cores'] = 1

        # TO WRITE FUNCTION FOR
        max_its = 300
        self.n_cores = settings['#cores']
        self.pixels_or_nm = settings['pixels_or_nm']
        self.slice, _ = self.parse_start_end(settings['frame_begin'], settings['frame_end'])

        if settings['method'] == "Phasor + Intensity":
            self.fitter = Phasor(settings)
        elif settings['method'] == "Phasor":
            self.fitter = PhasorDumb(settings)
        elif settings['method'] == "Gaussian - Fit bg":
            self.fitter = GaussianBackground(settings, max_its, 6)
        elif settings['method'] == "Gaussian - Estimate bg":
            self.fitter = Gaussian(settings, max_its, 5)
        else:
            self.fitter = PhasorSum(settings)

    def run(self):
        frames_to_fit = np.asarray(self.frames[self.slice])

        if "Phasor" in self.fitter.__name__:
            self.fitter.run(self, frames_to_fit, 0)
        else:
            results = np.array(1)
            if self.n_cores > 1:
                pass
            else:
                results = self.fitter.run(self, frames_to_fit, 0)

            for roi in self.active_rois:
                roi_results = results[results[:, 1] == roi.index, :]
                roi.results[self.name] = {"result": roi_results[:, 1:],
                                          "raw": roi.get_frame_stack(frames_to_fit, self.fitter.roi_size_1D)}


# %% Gaussian fitter with estimated background


TERMINATION_MESSAGES = {
    -1: "Improper input parameters status returned from `leastsq`",
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`ftol` termination condition is satisfied.",
    3: "`xtol` termination condition is satisfied.",
    4: "Both `ftol` and `xtol` termination conditions are satisfied."
}

FROM_MINPACK_TO_COMMON = {
    0: -1,  # Improper input parameters from MINPACK.
    1: 2,
    2: 3,
    3: 4,
    4: 1,
    5: 0
    # There are 6, 7, 8 for too small tolerance parameters,
    # but we guard against it by checking ftol, xtol, gtol beforehand.
}

EPS = np.finfo(float).eps


# noinspection DuplicatedCode
class Gaussian:
    """
    Gaussian fitter with estimated background, build upon Scipy Optimize Least-Squares
    """

    def __init__(self, settings, max_its, num_fit_params):
        """

        Parameters
        ----------
        :param settings: Fitting settings
        max_its: number of iterations limit
        :param num_fit_params: number of fitting parameters, 5 for without bg, 6 with bg

        Returns
        -------
        None.

        """
        self.result = []
        self.roi_locations = []
        self.roi_size = settings['roi_size']
        self.roi_size_1D = int((self.roi_size - 1) / 2)
        self.init_sig = 1.2  # Slightly on the high side probably
        self.threshold_method = settings['rejection']
        self.__name__ = settings['method']

        self.num_fit_params = num_fit_params
        self.params = np.zeros((1000, 2))

        self.x_scale = np.ones(num_fit_params)
        self.active_mask = np.zeros_like(self.x_scale, dtype=int)

        self.lb = np.resize(-np.inf, self.x_scale.shape)
        self.ub = np.resize(np.inf, self.x_scale.shape)
        self.use_one_sided = np.resize(False, self.x_scale.shape)
        self.empty_background = np.zeros(self.roi_size * 2 + (self.roi_size - 2) * 2, dtype=np.uint16)

        self.rel_step = EPS ** (1 / 3)
        self.comp = np.ones(num_fit_params)

        self.max_its = max_its

    def fun_find_max(self, roi):
        """
        Input ROI, returns max using FORTRAN

        Parameters
        ----------
        roi : region of interest

        Returns
        -------
        maximum

        """
        if self.roi_size == 9:
            return fortran_tools.max9(roi)
        elif self.roi_size == 7:
            return fortran_tools.max7(roi)

    def fun_find_min(self, roi):
        """
        Input ROI, returns minimum using FORTRAN

        Parameters
        ----------
        roi : region of interest

        Returns
        -------
        minimum

        """
        if self.roi_size == 9:
            return fortran_tools.min9(roi)
        elif self.roi_size == 7:
            return fortran_tools.min7(roi)

    def fun_calc_bg(self, roi):
        """
        Input ROI, returns background using FORTRAN

        Parameters
        ----------
        roi : region of interest

        Returns
        -------
        background

        """
        if self.roi_size == 9:
            return fortran_tools.calc_bg9(roi)
        elif self.roi_size == 7:
            return fortran_tools.calc_bg7(roi)

    def fun_norm(self, g):
        """
        Input 5 or 6 long array, returns norm using FORTRAN

        Parameters
        ----------
        g : vector of norm is to be found

        Returns
        -------
        norm of g

        """
        if self.num_fit_params == 5:
            return fortran_tools.norm5(g)
        elif self.num_fit_params == 6:
            return fortran_tools.norm6(g)

    def fun_gaussian(self, x, data):
        """
        Creates a gaussian with parameters x, returns difference between that and data using FORTRAN

        Parameters
        ----------
        x : Parameters of Gaussian
        data : ROI pixel values to be subtracted from created Gaussian

        Returns
        -------
        Difference between created 2D Gaussian and data

        """
        if self.num_fit_params == 5 and self.roi_size == 9:
            return fortran_linalg.gaussian(*x, self.roi_size, data)
        elif self.num_fit_params == 5 and self.roi_size == 7:
            return fortran_linalg.gaussian7(*x, self.roi_size, data)
        elif self.num_fit_params == 6 and self.roi_size == 9:
            return fortran_linalg.gs_bg(*x, self.roi_size, data)
        elif self.num_fit_params == 6 and self.roi_size == 7:
            return fortran_linalg.gs_bg7(*x, self.roi_size, data)

    def fun_jacobian(self, x0, data):
        """
        Generated jacobian using FORTRAN

        Parameters
        ----------
        x0 : Parameters of Gaussian
        data : ROI pixel values to be subtracted from created Gaussian

        Returns
        -------
        Jacobian

        """
        if self.num_fit_params == 5 and self.roi_size == 9:
            return fortran_linalg.dense_dif(x0, self.rel_step, self.comp,
                                            self.num_fit_params, self.roi_size, data)
        elif self.num_fit_params == 5 and self.roi_size == 7:
            return fortran_linalg.dense_dif7(x0, self.rel_step, self.comp,
                                             self.num_fit_params, self.roi_size, data)
        elif self.num_fit_params == 6 and self.roi_size == 9:
            return fortran_linalg.dense_dif_bg(x0, self.rel_step, self.comp,
                                               self.num_fit_params, self.roi_size, data)
        elif self.num_fit_params == 6 and self.roi_size == 7:
            return fortran_linalg.dense_dif_bg7(x0, self.rel_step, self.comp,
                                                self.num_fit_params, self.roi_size, data)

    def call_minpack(self, fun, x0, data, ftol, xtol, gtol, max_nfev, diag):
        """
        Caller of minpack, which is the iterative method of Python

        Parameters
        ----------
        fun : Function to be fitted
        x0 : Parameters
        data : Data to be compared to
        ftol : function tolerance
        xtol : parameter tolerance
        gtol : gradient tolerance
        max_nfev : max iterations
        diag : diagonal

        Returns
        -------
        Iterated result

        """
        n = x0.size
        epsfcn = EPS

        full_output = True
        factor = 100.0

        if max_nfev is None:
            # n squared to account for Jacobian evaluations.
            max_nfev = 100 * n * (n + 1)
        x, info, status = _minpack._lmdif(
            fun, x0, (), full_output, ftol, xtol, gtol,
            max_nfev, epsfcn, factor, diag)

        f = info['fvec']

        j = self.fun_jacobian(x0, data)

        cost = 0.5 * np.dot(f, f)
        g = j.T.dot(f)
        g_norm = self.fun_norm(g)

        nfev = info['nfev']
        njev = info.get('njev', None)

        status = FROM_MINPACK_TO_COMMON[status]
        active_mask = self.active_mask

        return OptimizeResult(
            x=x, cost=cost, fun=f, jac=j, grad=g, optimality=g_norm,
            active_mask=active_mask, nfev=nfev, njev=njev, status=status)

    def least_squares(self, x0, data, ftol=1e-8, xtol=1e-8, gtol=1e-8,
                      max_nfev=None):
        """
        The least-squares iterator. Does preparation before minpack iterator is called

        Parameters
        ----------
        x0 : parameters to be fitted
        data : data to be fitter to
        ftol : optional function tolerance. The default is 1e-8.
        xtol : optional parameter tolerance. The default is 1e-8.
        gtol : optional gradient tolerance. The default is 1e-8.
        max_nfev : optional max iterations. The default is None.

        Raises
        ------
        ValueError: in case you did anything wrong, should not happen if you just input a 2D Gaussian

        Returns
        -------
        Iterated result

        """

        def fun_wrapped(x):

            return self.fun_gaussian(x, data)

        if max_nfev is not None and max_nfev <= 0:
            raise ValueError("`max_nfev` must be None or positive integer.")

        if np.iscomplexobj(x0):
            raise ValueError("`x0` must be real.")

        if x0.ndim > 1:
            raise ValueError("`x0` must have at most 1 dimension.")

        f0 = fun_wrapped(x0)

        if f0.ndim != 1:
            raise ValueError("`fun` must return at most 1-d array_like.")

        # if not np.all(np.isfinite(f0)):
        #     raise ValueError("Residuals are not finite in the initial point.")

        result = self.call_minpack(fun_wrapped, x0, data, ftol, xtol, gtol,
                                   max_nfev, self.x_scale)

        result.message = TERMINATION_MESSAGES[result.status]
        result.success = result.status > 0

        return result

    def phasor_fit(self, data):
        """
        Does Phasor fit to estimate parameters using FORTRAN

        Parameters
        ----------
        data : ROI pixel values

        Returns
        -------
        pos_x : position in x-direction in ROI
        pos_y : position in y-direction in ROI

        """
        if self.roi_size == 9:
            x_re, x_im, y_re, y_im = fortran_tools.fft9(data)
        else:
            x_re, x_im, y_re, y_im = fortran_tools.fft7(data)

        roi_size = self.roi_size
        ang_x = atan2(x_im, x_re)
        if ang_x > 0:
            ang_x = ang_x - 2 * pi

        pos_x = abs(ang_x) / (2 * pi / roi_size) - 1

        ang_y = atan2(y_im, y_re)

        if ang_y > 0:
            ang_y = ang_y - 2 * pi

        pos_y = abs(ang_y) / (2 * pi / roi_size) - 1

        return pos_x, pos_y

    def phasor_guess(self, data):
        """
        Returns an initial guess based on phasor fitting

        Parameters
        ----------
        data : ROI pixel values

        Returns
        -------
        Parameters for Gaussian fitting. In order: height, position y, position x, y-sigma, x-sigma.

        """
        pos_x, pos_y = self.phasor_fit(data)
        height = data[int(pos_y), int(pos_x)] - self.fun_find_min(data)

        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig])

    def fit_gaussian(self, data, peak_index):
        """
        Gathers parameter estimate and calls fit.

        Parameters
        ----------
        data : ROI pixel values
        peak_index : Number of ROI that is being fitter.

        Returns
        -------
        p.x: solution of parameters
        p.nfev: number of iterations
        p.success: success or failure

        """
        if self.params[peak_index, 0] == 0:
            params = self.phasor_guess(data)
        else:
            params = self.phasor_guess(data)
            params[3:5] = self.params[peak_index, :]
        p = self.least_squares(params, data, max_nfev=self.max_its)  # , ftol=1e-10, xtol=1e-10, gtol=1e-10)

        return [p.x, p.nfev, p.success]

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
        pos_min = 0
        int_min = 0
        int_max = ((2 ** 16) - 1) * 1.5  # 50% margin over maximum pixel value possible
        sig_min = 0
        sig_max = self.roi_size_1D + 1

        return pos_max, pos_min, int_max, int_min, sig_max, sig_min

    def fitter(self, frame_index, frame, start_frame):
        """
        Fits all peaks by Gaussian least-squares fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame to be fitted
        start_frame: in case of multiprocessing, the number of the first frame passed to this thread/process

        Returns
        -------
        frame_result : fits of each frame

        """
        pos_max, pos_min, int_max, int_min, sig_max, sig_min = self.define_fitter_bounds()

        frame_result = np.zeros([len(self.roi_locations), 9])

        for roi in self.roi_locations:
            my_roi = roi.get_roi(frame, self.roi_size_1D)
            my_roi_bg = self.fun_calc_bg(my_roi)
            my_roi = my_roi - my_roi_bg

            result, its, success = self.fit_gaussian(my_roi, roi.index)

            if self.threshold_method == "None":
                if success == 0 or result[0] == 0:
                    self.params[roi.index, :] = [self.init_sig, self.init_sig]
                    success = 0
            else:
                if success == 0 or \
                        result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max or \
                        result[0] <= int_min or result[0] > int_max or \
                        result[3] < sig_min or result[3] > sig_max or result[4] < sig_min or result[4] > sig_max:
                    self.params[roi.index, :] = [self.init_sig, self.init_sig]
                    success = 0

            if success == 1:
                self.params[roi.index, :] = result[3:5]

                frame_result[roi.index, 0] = roi.index
                frame_result[roi.index, 1] = frame_index + start_frame
                # start position plus from center in ROI + half for indexing of pixels
                frame_result[roi.index, 2] = result[1] + roi.y - self.roi_size_1D + 0.5  # y
                frame_result[roi.index, 3] = result[2] + roi.x - self.roi_size_1D + 0.5  # x
                frame_result[roi.index, 4] = result[0]*result[3]*result[4]*2*pi  # Integrated intensity
                frame_result[roi.index, 5] = result[3]  # sigma y
                frame_result[roi.index, 6] = result[4]  # sigma x
                frame_result[roi.index, 7] = my_roi_bg
                frame_result[roi.index, 8] = its
            else:
                frame_result[roi.index, :] = np.nan
                frame_result[roi.index, 0] = roi.index
                frame_result[roi.index, 1] = frame_index + start_frame

        return frame_result

    def run(self, dataset, frames, start_frame):
        self.roi_locations = dataset.active_rois
        self.params = np.zeros((len(self.roi_locations), 2))
        try:
            n_frames = frames.shape[0]
        except AttributeError:
            n_frames = len(frames)

        tot_fits = 0

        for frame_index, frame in enumerate(frames):
            if frame_index == 0:
                frame_result = self.fitter(frame_index, frame, start_frame)
                n_fits = frame_result.shape[0]
                width = frame_result.shape[1]
                height = len(frames) * len(self.roi_locations)
                self.result = np.zeros((height, width))
                self.result[tot_fits:tot_fits + n_fits, :] = frame_result
                tot_fits += n_fits
            else:
                frame_result = self.fitter(frame_index, frame, start_frame)
                n_fits = frame_result.shape[0]
                self.result[tot_fits:tot_fits + n_fits, :] = frame_result
                tot_fits += n_fits

            if frame_index % (round(n_frames / 10, 0)) == 0:
                dataset.experiment.progress_function(frame_index + 1, n_frames + 1)

        dataset.experiment.progress_function(n_frames + 1, n_frames + 1)

        return self.result

# %% Gaussian fitter including background


# noinspection DuplicatedCode
class GaussianBackground(Gaussian):
    """
    Gaussian fitter with fitted background, build upon Scipy Optimize Least-Squares
    """

    def phasor_guess(self, data):
        """
        Returns an initial guess based on phasor fitting

        Parameters
        ----------
        data : ROI pixel values

        Returns
        -------
        Parameters for Gaussian fitting. In order: height, position y, position x, y-sigma, x-sigma, background.

        """
        pos_x, pos_y = self.phasor_fit(data)
        background = self.fun_calc_bg(data)
        height = data[int(pos_y), int(pos_x)] - background

        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig, background])

    def fitter(self, frame_index, frame, start_frame):
        """
        Fits all peaks by Gaussian least-squares fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame to be fitted
        start_frame: in case of multiprocessing, the number of the first frame passed to this thread/process

        Returns
        -------
        frame_result : fits of each frame

        """
        pos_max, pos_min, int_max, int_min, sig_max, sig_min = self.define_fitter_bounds()

        frame_result = np.zeros([len(self.roi_locations), 9])

        for roi in self.roi_locations:
            my_roi = roi.get_roi(frame, self.roi_size_1D)
            result, its, success = self.fit_gaussian(my_roi, roi.index)

            if self.threshold_method == "None":
                if success == 0 or result[0] == 0:
                    self.params[roi.index, :] = [self.init_sig, self.init_sig]
                    success = 0
            else:
                if success == 0 or \
                        result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max or \
                        result[0] <= int_min or result[0] > int_max or \
                        result[3] < sig_min or result[3] > sig_max or result[4] < sig_min or result[4] > sig_max:
                    self.params[roi.index, :] = [self.init_sig, self.init_sig]
                    success = 0

            if success == 1:
                self.params[roi.index, :] = result[3:5]

                frame_result[roi.index, 0] = roi.index
                frame_result[roi.index, 1] = frame_index + start_frame
                # start position plus from center in ROI + half for indexing of pixels
                frame_result[roi.index, 2] = result[1] + roi.y - self.roi_size_1D + 0.5  # y
                frame_result[roi.index, 3] = result[2] + roi.x - self.roi_size_1D + 0.5  # x
                frame_result[roi.index, 4] = result[0]*result[3]*result[4]*2*pi  # Integrated intensity
                frame_result[roi.index, 5] = result[3]  # sigma y
                frame_result[roi.index, 6] = result[4]  # sigma x
                frame_result[roi.index, 7] = result[5]  # background
                frame_result[roi.index, 8] = its
            else:
                frame_result[roi.index, :] = np.nan
                frame_result[roi.index, 0] = roi.index
                frame_result[roi.index, 1] = frame_index + start_frame

        return frame_result


# %% Phasor for ROI loops


class Phasor:
    """
    Phasor fitting using Fourier Transform. Also returns intensity of pixel in which Phasor position is found.
    """

    def __init__(self, settings):
        """
        Parameters
        ----------
        :param settings: input settings

        Returns
        -------
        None.

        """
        self.result = []
        self.roi_locations = []
        self.roi_size = settings['roi_size']
        self.roi_size_1D = int((self.roi_size - 1) / 2)
        self.threshold_method = settings['rejection']
        self.__name__ = settings['method']

    def fun_find_max(self, roi):
        """
        Input ROI, returns max using FORTRAN

        Parameters
        ----------
        roi : region of interest

        Returns
        -------
        maximum

        """
        if self.roi_size == 9:
            return fortran_tools.max9(roi)
        elif self.roi_size == 7:
            return fortran_tools.max7(roi)

    def fun_calc_bg(self, roi):
        """
        Input ROI, returns background using FORTRAN

        Parameters
        ----------
        roi : region of interest

        Returns
        -------
        background

        """
        if self.roi_size == 9:
            return fortran_tools.calc_bg9(roi)
        elif self.roi_size == 7:
            return fortran_tools.calc_bg7(roi)

    def fft_to_pos(self, fft_values):
        """
        Convert the found Fourier transform to a Phasor position

        Parameters
        ----------
        fft_values : Fourier transform of ROI

        Returns
        -------
        pos_x : x-position of Phasor
        pos_y : y-position of Phasor

        """
        ang_x = phase(fft_values[0, 1])
        if ang_x > 0:
            ang_x = ang_x - 2 * pi

        pos_x = abs(ang_x) / (2 * pi / self.roi_size) + 0.5

        ang_y = phase(fft_values[1, 0])

        if ang_y > 0:
            ang_y = ang_y - 2 * pi

        pos_y = abs(ang_y) / (2 * pi / self.roi_size) + 0.5

        return pos_x, pos_y

    def phasor_fit_stack(self, frame_stack, roi_index, y, x):
        """
        Applies phasor fitting to an entire stack of frames of one ROI

        Parameters
        ----------
        frame_stack : stack of frames of a single ROI
        roi_index : Index of ROI that is currently being fitted
        y : y-location with respect to total microscope frame of currently fitted ROI
        x : x-location with respect to total microscope frame of currently fitted ROI

        Returns
        -------
        roi_result : Result of all the frames of the current ROI

        """
        roi_result = np.zeros([frame_stack.shape[0], 5])

        roi_bb = empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = empty_aligned((frame_stack.shape[0], self.roi_size, self.roi_size_1D + 1), dtype='complex128')
        fft_values_list = FFTW(roi_bb, roi_bf, axes=(1, 2),
                               flags=('FFTW_MEASURE',),
                               direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)

        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):
            success = 1
            my_frame_bg = self.fun_calc_bg(frame)
            frame_max = self.fun_find_max(frame)

            pos_x, pos_y = self.fft_to_pos(fft_values)

            if self.threshold_method != "None" and (pos_x > self.roi_size or pos_x < 0):
                success = 0
            if self.threshold_method != "None" and (pos_y > self.roi_size or pos_y < 0):
                success = 0
            if frame_max == 0:
                success = 0

            if success == 1:
                roi_result[frame_index, 0] = frame_index
                # start position plus from center in ROI + half for indexing of pixels
                roi_result[frame_index, 1] = y + pos_y - self.roi_size_1D  # y
                roi_result[frame_index, 2] = x + pos_x - self.roi_size_1D  # x
                roi_result[frame_index, 3] = frame_max - my_frame_bg  # returns max peak
                roi_result[frame_index, 4] = my_frame_bg  # background
            else:
                roi_result[frame_index, :] = np.nan
                roi_result[frame_index, 0] = frame_index

        return roi_result

    def run(self, dataset, frames, _):
        self.roi_locations = dataset.active_rois

        tot_fits = 0

        for roi in self.roi_locations:
            frame_stack = roi.get_frame_stack(frames, self.roi_size_1D)
            roi_result = self.phasor_fit_stack(frame_stack, roi.index, roi.y, roi.x)

            result_dict = {"result": roi_result, "raw": frame_stack}
            roi.results[dataset.name] = result_dict
            tot_fits += roi_result.shape[0]

            if len(self.roi_locations) > 10:
                if roi.index % (round(len(self.roi_locations) / 10, 0)) == 0:
                    dataset.experiment.progress_function(roi.index + 1, len(self.roi_locations) + 1)

        dataset.experiment.progress_function(len(self.roi_locations) + 1, len(self.roi_locations) + 1)


# %% Dumb phasor ROI loop


class PhasorDumb(Phasor):
    """
    Dumb Phasor. Does not return intensity of found location.
    """

    def phasor_fit_stack(self, frame_stack, roi_index, y, x):
        """
        Applies phasor fitting to an entire stack of frames of one ROI

        Parameters
        ----------
        frame_stack : stack of frames of a single ROI
        roi_index : Index of ROI that is currently being fitted
        y : y-location with respect to total microscope frame of currently fitted ROI
        x : x-location with respect to total microscope frame of currently fitted ROI

        Returns
        -------
        roi_result : Result of all the frames of the current ROI

        """
        roi_result = np.zeros([frame_stack.shape[0], 3])

        roi_bb = empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = empty_aligned((frame_stack.shape[0], self.roi_size, self.roi_size_1D + 1), dtype='complex128')
        fft_values_list = FFTW(roi_bb, roi_bf, axes=(1, 2),
                               flags=('FFTW_MEASURE',),
                               direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)

        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):
            success = 1
            pos_x, pos_y = self.fft_to_pos(fft_values)

            if self.threshold_method != "None" and (pos_x > self.roi_size or pos_x < 0):
                success = 0
            if self.threshold_method != "None" and (pos_y > self.roi_size or pos_y < 0):
                success = 0

            if success == 1:
                roi_result[frame_index, 0] = frame_index
                # start position plus from center in ROI + half for indexing of pixels
                roi_result[frame_index, 1] = y + pos_y - self.roi_size_1D  # y
                roi_result[frame_index, 2] = x + pos_x - self.roi_size_1D  # x
            else:
                roi_result[frame_index, :] = np.nan
                roi_result[frame_index, 0] = frame_index

        return roi_result


# %% Phasor with sum


class PhasorSum(Phasor):
    """
    Phasor Sum. Also returns summation of entire ROI
    """

    def phasor_fit_stack(self, frame_stack, roi_index, y, x):
        """
        Applies phasor fitting to an entire stack of frames of one ROI

        Parameters
        ----------
        frame_stack : stack of frames of a single ROI
        roi_index : Index of ROI that is currently being fitted
        y : y-location with respect to total microscope frame of currently fitted ROI
        x : x-location with respect to total microscope frame of currently fitted ROI

        Returns
        -------
        roi_result : Result of all the frames of the current ROI

        """
        roi_result = np.zeros([frame_stack.shape[0], 4])

        roi_bb = empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = empty_aligned((frame_stack.shape[0], self.roi_size, self.roi_size_1D + 1), dtype='complex128')
        fft_values_list = FFTW(roi_bb, roi_bf, axes=(1, 2),
                               flags=('FFTW_MEASURE',),
                               direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)

        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):
            success = 1
            frame_sum = np.sum(frame)

            pos_x, pos_y = self.fft_to_pos(fft_values)

            if self.threshold_method != "None" and (pos_x > self.roi_size or pos_x < 0):
                success = 0
            if self.threshold_method != "None" and (pos_y > self.roi_size or pos_y < 0):
                success = 0
            if frame_sum == 0:
                success = 0

            if success == 1:
                roi_result[frame_index, 0] = frame_index
                # start position plus from center in ROI + half for indexing of pixels
                roi_result[frame_index, 1] = y + pos_y - self.roi_size_1D  # y
                roi_result[frame_index, 2] = x + pos_x - self.roi_size_1D  # x
                roi_result[frame_index, 3] = frame_sum  # returns summation
            else:
                roi_result[frame_index, :] = np.nan
                roi_result[frame_index, 0] = frame_index

        return roi_result
