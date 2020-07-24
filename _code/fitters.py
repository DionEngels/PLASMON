# -*- coding: utf-8 -*-
"""
Created on Sun May 31 22:58:08 2020

@author: Dion Engels
MBx Python Data Analysis

Gaussian fitting

----------------------------

v0.1: Setup: 31/05/2020
v0.2: Bugged rainSTORM inspired: 03/06/2020
v1.0: rainSTORM slow, but working: 04/06/2020
v1.1: first set of solvers: 15/06/2020
v2.0: optimilzations after initial speed check: 17/06/2020
v2.1: ran second set of solver performance: 21/06/2020
v2.2: phasor fitter over ROI, dumb phasor fitter,
      ROI mover, log gauss fitter, background fitting variable: 22/06/2020
v3.0: optimization after second speed check. New fitters:
      SettingsScipy, ScipyROIupdater, ScipyLog, ScipyFrameLastFit
      ScipyROIloopBackground, PhasorROILoop, PhasorROILoopDum: 27/06/2020
v3.1: Last fit frame loop with new method: 28/06/2020
v4.0: removed all unneeded fitters
v4.1: Dions fitter v1
v4.2: Dions fitter v2
v4.3: small optimization of Scipy
v4.4: attempt of stripped least squares function
v4.5: cached: 06/07/2020
v4.6: cached class
v4.7: added approx_derivative to class
v4.8: passing f0 through the class
v4.9: new method of background determination
v4.10: own dense difference
v4.11: only check intensity = 0 for save params
v4.12: also cache for derivative
v4.13: cleanup
v4.14: intensity invariant cache
v4.15: intensity back in cache
v4.16: background in fit
v4.17: further cleanup
v5.0: Three versions of cached fitter
v5.1: limited size of cache
v5.2: caching bugfix v1
v5.3: invariant caching test
v5.4: complete clear test
v5.5: complete clear FORTRAN
v5.6: invariant cache FORTRAN
v5.7: small improvemnet, remove derivate cache
v5.8: improvement fun_wrapped, deletion old gaussian
v5.9: without caching
v5.10: alternative last fits
v5.11: data to FORTRAN
v5.12: data to FORTAN v2
v5.13: dense_differnece in FORTRAN
v5.14: cleanup of FORTRAN and Python code
v6.0: FORTRAN enabled without cache, ST
v6.1: cutoff intensity
v6.2: stricter rejection and only save params if succesful
v6.3: threshold from ROI finder
v6.4: ROI size 7 implementation
v6.5: own norm
v6.6: FORTRAN tools: max, norm, and calc bg
v7.0: cleanup and declare_functions function: 14/07/2020
v7.1: revert lambda implementation for MT
v7.2: change in phasor guess
v7.3: further change in phasor bounds.
v7.4: changes for MP status updates
v7.5: removed any wavelength dependancy
v7.6: add FORTRAN fft and min
v7.7: variance as fit test
v7.8: back to sigma as fit
v7.9: max its of 100, working with new ROI finder
v8.0: PhasorSum, remove LastFit, rename, and clean up
v8.1: Rejection options
v8.2: change rejection options Phasor fits

"""
# %% Generic imports
from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.optimize import _minpack, OptimizeResult
from scipy.optimize._lsq.common import EPS
import _code.MBx_FORTRAN_v5 as fortran_linalg
import _code.MBx_FORTRAN_TOOLS_v2 as fortran_tools
from pyfftw import empty_aligned, FFTW  # for FFT
from pims import Frame  # for converting ND2 generator to one numpy array
from math import pi, atan2
from cmath import phase

# %% Scipy last fit guess exlcuding background

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


class Gaussian:
    """
    Gaussian fitter with estimated background, build upon Scipy Optimize Least Squares
    """

    def __init__(self, roi_size, thresholds, threshold_method, method, num_fit_params):
        """

        Parameters
        ----------
        roi_size : ROI size
        thresholds : minimum and maximum of fits
        threshold_method: way to apply threshold

        Returns
        -------
        None.

        """
        self.result = []
        self.roi_locations = np.ones((1, 2))
        self.roi_size = roi_size
        self.roi_size_1D = int((self.roi_size - 1) / 2)
        self.init_sig = 1.3  # Slightly on the high side probably
        self.thresholds = thresholds
        self.threshold_method = threshold_method
        self.__name__ = method

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

    def fun_find_max(self, roi):

        if self.roi_size == 9:
            return fortran_tools.max9(roi)
        elif self.roi_size == 7:
            return fortran_tools.max7(roi)

    def fun_find_min(self, roi):

        if self.roi_size == 9:
            return fortran_tools.min9(roi)
        elif self.roi_size == 7:
            return fortran_tools.min7(roi)

    def fun_calc_bg(self, roi):

        if self.roi_size == 9:
            return fortran_tools.calc_bg9(roi)
        elif self.roi_size == 7:
            return fortran_tools.calc_bg7(roi)

    def fun_norm(self, g):

        if self.num_fit_params == 5:
            return fortran_tools.norm5(g)
        elif self.num_fit_params == 6:
            return fortran_tools.norm6(g)

    def fun_gaussian(self, x, data):

        if self.num_fit_params == 5 and self.roi_size == 9:
            return fortran_linalg.gaussian(*x, self.roi_size, data)
        elif self.num_fit_params == 5 and self.roi_size == 7:
            return fortran_linalg.gaussian7(*x, self.roi_size, data)
        elif self.num_fit_params == 6 and self.roi_size == 9:
            return fortran_linalg.gs_bg(*x, self.roi_size, data)
        elif self.num_fit_params == 6 and self.roi_size == 7:
            return fortran_linalg.gs_bg7(*x, self.roi_size, data)

    def fun_jacobian(self, x0, data):

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

    def phasor_fit(self, data):
        """
        Parameters
        ----------
        data : input ROI
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

    def call_minpack(self, fun, x0, data, ftol, xtol, gtol, max_nfev, diag):

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

    def phasor_guess(self, data):
        """ Returns an initial guess based on phasor fitting"""
        pos_x, pos_y = self.phasor_fit(data)
        height = data[int(pos_y), int(pos_x)] - self.fun_find_min(data)

        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig])

    def fit_gaussian(self, data, peak_index):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""

        if self.params[peak_index, 0] == 0:
            params = self.phasor_guess(data)
        else:
            params = self.phasor_guess(data)
            params[3:5] = self.params[peak_index, :]
        p = self.least_squares(params, data, max_nfev=100)  # , gtol=1e-4, ftol=1e-4)

        return [p.x, p.nfev, p.success]

    def define_fitter_bounds(self):

        if self.threshold_method == "Strict":
            pos_max = self.roi_size
            pos_min = 0
            int_min = self.thresholds['int_min']
            int_max = self.thresholds['int_max']
            sig_min = self.thresholds['sigma_min']
            sig_max = self.thresholds['sigma_max']
        else:
            pos_max = self.roi_size
            pos_min = 0
            int_min = 0
            int_max = (2 ** 16) - 1
            sig_min = 0
            sig_max = self.roi_size_1D + 1

        return pos_max, pos_min, int_max, int_min, sig_max, sig_min

    def fitter(self, frame_index, frame, start_frame):
        """
        Fits all peaks by scipy least squares fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        start_frame:

        Returns
        -------
        frame_result : fits of each frame

        """

        pos_max, pos_min, int_max, int_min, sig_max, sig_min = self.define_fitter_bounds()

        frame_result = np.zeros([self.roi_locations.shape[0], 9])

        for peak_index, peak in enumerate(self.roi_locations):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y - self.roi_size_1D:y + self.roi_size_1D + 1, x - self.roi_size_1D:x + self.roi_size_1D + 1]
            my_roi_bg = self.fun_calc_bg(my_roi)
            my_roi = my_roi - my_roi_bg

            if self.threshold_method == "Strict" and self.fun_find_max(my_roi) < self.thresholds['pixel_min']:
                continue

            result, its, success = self.fit_gaussian(my_roi, peak_index)

            if self.threshold_method == "None":
                if success == 0:
                    continue
            else:
                if success == 0 or \
                        result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max or \
                        result[0] < int_min or result[0] > int_max or \
                        result[3] < sig_min or result[3] > sig_max or result[4] < sig_min or result[4] > sig_max:
                    continue

            self.params[peak_index, :] = result[3:5]

            frame_result[peak_index, 0] = frame_index + start_frame
            frame_result[peak_index, 1] = peak_index
            # start position plus from center in ROI + half for indexing of pixels
            frame_result[peak_index, 2] = result[1] + y - self.roi_size_1D + 0.5
            frame_result[peak_index, 3] = result[2] + x - self.roi_size_1D + 0.5
            frame_result[peak_index, 4] = result[0]
            frame_result[peak_index, 5] = result[3]
            frame_result[peak_index, 6] = result[4]
            frame_result[peak_index, 7] = my_roi_bg
            frame_result[peak_index, 8] = its

        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result

    def main(self, frames, metadata, roi_locations, gui=None, q=None,
             start_frame=0, n_frames=None):
        """
        Main for every fitter method, calls fitter function and returns fits

        Parameters
        ----------
        frames : All frames of .nd2 video
        metadata : Metadata of .nd2 video
        start_frame:
        n_frames:
        gui:
        q:
        roi_locations:

        Returns
        -------
        All localizations

        """
        if n_frames is None:
            n_frames = metadata['sequence_count']

        self.roi_locations = roi_locations
        self.params = np.zeros((self.roi_locations.shape[0], 2))

        tot_fits = 0

        for frame_index, frame in enumerate(frames):
            if frame_index == 0:
                frame_result = self.fitter(frame_index, frame, start_frame)
                n_fits = frame_result.shape[0]
                width = frame_result.shape[1]
                height = len(frames) * self.roi_locations.shape[0]
                self.result = np.zeros((height, width))
                self.result[tot_fits:tot_fits + n_fits, :] = frame_result
                tot_fits += n_fits
            else:
                frame_result = self.fitter(frame_index, frame, start_frame)
                n_fits = frame_result.shape[0]
                self.result[tot_fits:tot_fits + n_fits, :] = frame_result
                tot_fits += n_fits

            if frame_index % (round(n_frames / 10, 0)) == 0:
                if gui is not None:
                    gui.update_status(frame_index + 1, n_frames + 1)
                elif q is not None:
                    q.put([start_frame, frame_index + 1])
                else:
                    print('Done fitting frame ' + str(frame_index) + ' of ' + str(n_frames))

        self.result = np.delete(self.result, range(tot_fits, len(frames) * self.roi_locations.shape[0]), axis=0)

        if gui is not None:
            gui.update_status(n_frames, n_frames)
        if q is not None:
            q.put([start_frame, len(frames) + 1])

        return self.result


# %% Gaussian fitter including background


class GaussianBackground(Gaussian):

    def phasor_guess(self, data):
        """ Returns an initial guess based on phasor fitting"""
        pos_x, pos_y = self.phasor_fit(data)
        background = self.fun_calc_bg(data)
        height = data[int(pos_y), int(pos_x)] - background

        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig, background])

    def fitter(self, frame_index, frame, start_frame):
        """
        Fits all peaks by scipy least squares fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        start_frame :

        Returns
        -------
        frame_result : fits of each frame

        """

        pos_max, pos_min, int_max, int_min, sig_max, sig_min = self.define_fitter_bounds()

        frame_result = np.zeros([self.roi_locations.shape[0], 9])

        for peak_index, peak in enumerate(self.roi_locations):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y - self.roi_size_1D:y + self.roi_size_1D + 1, x - self.roi_size_1D:x + self.roi_size_1D + 1]

            if self.threshold_method == "Strict" and \
                    self.fun_find_max(my_roi) - self.fun_calc_bg(my_roi) < self.thresholds['pixel_min']:
                continue

            result, its, success = self.fit_gaussian(my_roi, peak_index)

            if self.threshold_method == "None":
                if success == 0:
                    continue
            else:
                if success == 0 or \
                        result[2] < pos_min or result[2] > pos_max or result[1] < pos_min or result[1] > pos_max or \
                        result[0] < int_min or result[0] > int_max or \
                        result[3] < sig_min or result[3] > sig_max or result[4] < sig_min or result[4] > sig_max:
                    continue

            self.params[peak_index, :] = result[3:5]

            frame_result[peak_index, 0] = frame_index + start_frame
            frame_result[peak_index, 1] = peak_index
            # start position plus from center in ROI + half for indexing of pixels
            frame_result[peak_index, 2] = result[1] + y - self.roi_size_1D + 0.5  # y
            frame_result[peak_index, 3] = result[2] + x - self.roi_size_1D + 0.5  # x
            frame_result[peak_index, 4] = result[0]
            frame_result[peak_index, 5] = result[3]
            frame_result[peak_index, 6] = result[4]
            frame_result[peak_index, 7] = result[5]
            frame_result[peak_index, 8] = its

        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result


# %% Phasor for ROI loops


class Phasor:
    """
    Phasor localizer over ROIs
    """

    def __init__(self, roi_size, thresholds, threshold_method, method):
        """
        Parameters
        ----------
        roi_size : ROI size
        thresholds : minimum and maximum of fits
        threshold_method: way to apply threshold
        method : Name of method

        Returns
        -------
        None.

        """

        self.result = []
        self.roi_locations = []
        self.roi_size = roi_size
        self.roi_size_1D = int((self.roi_size - 1) / 2)
        self.thresholds = thresholds
        self.threshold_method = threshold_method
        self.__name__ = method

    def fun_find_max(self, roi):

        if self.roi_size == 9:
            return fortran_tools.max9(roi)
        elif self.roi_size == 7:
            return fortran_tools.max7(roi)

    def fun_calc_bg(self, roi):

        if self.roi_size == 9:
            return fortran_tools.calc_bg9(roi)
        elif self.roi_size == 7:
            return fortran_tools.calc_bg7(roi)

    def fft_to_pos(self, fft_values):

        ang_x = phase(fft_values[0, 1])
        if ang_x > 0:
            ang_x = ang_x - 2 * pi

        pos_x = abs(ang_x) / (2 * pi / self.roi_size) + 0.5

        ang_y = phase(fft_values[1, 0])

        if ang_y > 0:
            ang_y = ang_y - 2 * pi

        pos_y = abs(ang_y) / (2 * pi / self.roi_size) + 0.5

        return pos_x, pos_y

    def phasor_fit_stack(self, frame_stack, roi, roi_index, y, x):

        roi_result = np.zeros([frame_stack.shape[0], 6])

        roi_bb = empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = empty_aligned((frame_stack.shape[0], self.roi_size, self.roi_size_1D + 1), dtype='complex128')
        fft_values_list = FFTW(roi_bb, roi_bf, axes=(1, 2),
                               flags=('FFTW_MEASURE',),
                               direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)

        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):

            my_frame_bg = self.fun_calc_bg(frame)
            frame_max = self.fun_find_max(frame)

            if self.threshold_method == "Strict" and frame_max < self.thresholds['pixel_min']:
                continue

            pos_x, pos_y = self.fft_to_pos(fft_values)

            if self.threshold_method != "None" and (pos_x > self.roi_size or pos_x < 0):
                continue
            if self.threshold_method != "None" and (pos_y > self.roi_size or pos_y < 0):
                continue

            roi_result[frame_index, 0] = frame_index
            roi_result[frame_index, 1] = roi_index
            # start position plus from center in ROI + half for indexing of pixels
            roi_result[frame_index, 2] = y + pos_y - self.roi_size_1D
            roi_result[frame_index, 3] = x + pos_x - self.roi_size_1D
            roi_result[frame_index, 4] = frame_max - my_frame_bg  # returns max peak
            roi_result[frame_index, 5] = my_frame_bg  # background

        roi_result = roi_result[roi_result[:, 2] > 0]

        return roi_result

    def main(self, frames, metadata, roi_locations, gui=None, n_frames=None):
        """
        Main for phasor over ROI loop, calls fitter function and returns fits.

        Parameters
        ----------
        frames : All frames of .nd2 video
        metadata : Metadata of .nd2 video
        roi_locations:
        gui:
        n_frames:

        Returns
        -------
        All localizations

        """
        self.roi_locations = roi_locations

        frame_stack_total = Frame(frames)

        tot_fits = 0
        created = 0

        for roi_index, roi in enumerate(self.roi_locations):
            y = int(roi[0])
            x = int(roi[1])
            frame_stack = frame_stack_total[:, y - self.roi_size_1D:y + self.roi_size_1D + 1,
                          x - self.roi_size_1D:x + self.roi_size_1D + 1]

            roi_result = self.phasor_fit_stack(frame_stack, roi, roi_index, y, x)

            if created == 0:
                n_fits = roi_result.shape[0]
                width = roi_result.shape[1]
                height = len(frames) * self.roi_locations.shape[0]
                self.result = np.zeros((height, width))
                created = 1
                self.result[tot_fits:tot_fits + n_fits, :] = roi_result
                tot_fits += n_fits
            else:
                n_fits = roi_result.shape[0]
                self.result[tot_fits:tot_fits + n_fits, :] = roi_result
                tot_fits += n_fits

            if self.roi_locations.shape[0] > 10:
                if roi_index % (round(self.roi_locations.shape[0] / 10, 0)) == 0:
                    if gui is None:
                        print('Done fitting ROI ' + str(roi_index) + ' of ' + str(self.roi_locations.shape[0]))
                    else:
                        gui.update_status(roi_index + 1, len(self.roi_locations) + 1)

        self.result = np.delete(self.result, range(tot_fits, len(frames) * self.roi_locations.shape[0]), axis=0)

        if gui is not None:
            gui.update_status(len(self.roi_locations), len(self.roi_locations))

        return self.result


# %% Dumb phasor ROI loop

class PhasorDumb(Phasor):

    def phasor_fit_stack(self, frame_stack, roi, roi_index, y, x):

        roi_result = np.zeros([frame_stack.shape[0], 4])

        roi_bb = empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = empty_aligned((frame_stack.shape[0], self.roi_size, self.roi_size_1D + 1), dtype='complex128')
        fft_values_list = FFTW(roi_bb, roi_bf, axes=(1, 2),
                               flags=('FFTW_MEASURE',),
                               direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)

        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):

            pos_x, pos_y = self.fft_to_pos(fft_values)

            if self.threshold_method != "None" and (pos_x > self.roi_size or pos_x < 0):
                continue
            if self.threshold_method != "None" and (pos_y > self.roi_size or pos_y < 0):
                continue

            roi_result[frame_index, 0] = frame_index
            roi_result[frame_index, 1] = roi_index
            # start position plus from center in ROI + half for indexing of pixels
            roi_result[frame_index, 2] = y + pos_y - self.roi_size_1D
            roi_result[frame_index, 3] = x + pos_x - self.roi_size_1D

        roi_result = roi_result[roi_result[:, 2] > 0]

        return roi_result


# %% Phasor with sum


class PhasorSum(Phasor):

    def phasor_fit_stack(self, frame_stack, roi, roi_index, y, x):

        roi_result = np.zeros([frame_stack.shape[0], 5])

        roi_bb = empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = empty_aligned((frame_stack.shape[0], self.roi_size, self.roi_size_1D + 1), dtype='complex128')
        fft_values_list = FFTW(roi_bb, roi_bf, axes=(1, 2),
                               flags=('FFTW_MEASURE',),
                               direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)

        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):

            frame_sum = np.sum(frame)

            pos_x, pos_y = self.fft_to_pos(fft_values)

            if self.threshold_method != "None" and (pos_x > self.roi_size or pos_x < 0):
                continue
            if self.threshold_method != "None" and (pos_y > self.roi_size or pos_y < 0):
                continue

            roi_result[frame_index, 0] = frame_index
            roi_result[frame_index, 1] = roi_index
            # start position plus from center in ROI + half for indexing of pixels
            roi_result[frame_index, 2] = y + pos_y - self.roi_size_1D
            roi_result[frame_index, 3] = x + pos_x - self.roi_size_1D
            roi_result[frame_index, 4] = frame_sum  # returns summation

        roi_result = roi_result[roi_result[:, 2] > 0]

        return roi_result
