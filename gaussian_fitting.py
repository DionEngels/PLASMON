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
"""
#%% Generic imports
from __future__ import division, print_function, absolute_import
import math
import numpy as np
import gauss_full4 as gauss2

#%% Base Phasor

from math import pi
import cmath

class base_phasor():
    
    def phasor_fit(self, data):
        """
        Parameters
        ----------
        data : input ROI
        
        return: pos x, pos y
        """
        
        fft_values = np.fft.fft2(data)

        roi_size = self.ROI_size
        ang_x = cmath.phase(fft_values[0, 1])
        if ang_x>0:
            ang_x=ang_x-2*pi

        pos_x = abs(ang_x)/(2*pi/roi_size)

        ang_y = cmath.phase(fft_values[1,0])

        if ang_y >0:
            ang_y = ang_y - 2*pi

        pos_y = abs(ang_y)/(2*pi/roi_size)

        if pos_x > 8.5:
            pos_x -= roi_size
        if pos_y > 8.5:
            pos_y -= roi_size
            
        return pos_x, pos_y
    
       
#%% Limited size cache
from collections import OrderedDict

class LimitedSizeDict(OrderedDict):
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        #self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)
                
#%% Cached scipy last fit guess exlcuding background 
from numpy.linalg import norm
from scipy.optimize import _minpack, OptimizeResult
from scipy.optimize._lsq.common import EPS

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

class scipy_last_fit_guess(base_phasor):
    """
    Build-in scipy least squares fitting, with last fit as initial guess
    """
    
    def __init__(self, metadata, ROI_size, wavelength, threshold, ROI_locations, METHOD, num_fit_params):
        """

        Parameters
        ----------
        metadata : metadata of ND2
        ROI_size : ROI size
        wavelength : laser wavelength
        threshold : # of standard deviations for threhsold
        ROI_locations : found ROIs

        Returns
        -------
        None.

        """
        self.result = []
        self.ROI_size = ROI_size
        self.ROI_size_1D = int((self.ROI_size-1)/2)
        self.ROI_locations = ROI_locations
        self.init_sig = wavelength/(2*metadata['NA']*math.sqrt(8*math.log(2)))/(metadata['calibration_um']*1000)*2 #*2 for better guess
        self.threshold_sigma = threshold
        self.__name__ = METHOD
        
        self.params = np.zeros((self.ROI_locations.shape[0],num_fit_params))
        
        self.indices = np.indices((ROI_size, ROI_size))
        self.cache = LimitedSizeDict(size_limit = 1000000) # 1000000 = 150 Mb
        self.x_scale = np.ones(num_fit_params)
        self.active_mask = np.zeros_like(self.x_scale, dtype=int)
        
        self.lb = np.resize(-np.inf, self.x_scale.shape)
        self.ub = np.resize(np.inf, self.x_scale.shape)
        self.use_one_sided = np.resize(False, self.x_scale.shape)
        self.empty_background = np.zeros(self.ROI_size*2+(self.ROI_size-2)*2, dtype=np.uint16)
                
        self.counter_cache = 0
        self.counter_calc = 0
        
    def dense_difference(self, fun, x0, f0, h):
        m = f0.size
        n = x0.size
        J_transposed = np.empty((n, m))
        h_vecs = np.diag(h)
        
        for i in range(h.size):
            x1 = x0 - h_vecs[i]
            x2 = x0 + h_vecs[i]
            dx = x2[i] - x1[i]
            f1 = fun(x1)
            f2 = fun(x2)
            df = f2 - f1
            
            J_transposed[i] = df / dx
        
        return J_transposed.T
        
    
    def compute_absolute_step(self, x0):
        rel_step = EPS**(1/3)
        return np.abs(rel_step * np.maximum(1.0, np.abs(x0)))
        
    def approx_derivative(self, fun, x0, f0, data):
         
        if x0.ndim > 1:
            raise ValueError("`x0` must have at most 1 dimension.")
    
        lb = self.lb
        ub = self.ub
    
        if lb.shape != x0.shape or ub.shape != x0.shape:
            raise ValueError("Inconsistent shapes between bounds and `x0`.")
    
        def fun_wrapped(x):
            
            key = tuple(x[1:5]) 
            if key in self.cache:
                self.counter_cache+=1
                return (self.cache[key]*x[0]+x[5]) - data
            else:
                result = gauss2.gaussian(*x, self.ROI_size)
                self.cache[key] = (result-x[5])/x[0]
                self.counter_calc+=1
                return result - data
       
        h = self.compute_absolute_step(x0)

        return self.dense_difference(fun_wrapped, x0, f0, h)
                
    def call_minpack(self, fun, x0, f0, data, jac, ftol, xtol, gtol, max_nfev, diag):

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
        
        j = self.approx_derivative(fun, x, f0, data)
        
        cost = 0.5 * np.dot(f, f)
        g = j.T.dot(f)
        g_norm = norm(g, ord=np.inf)
    
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
            
            key = tuple(x[1:5]) 
            if key in self.cache:
                self.counter_cache+=1
                return (self.cache[key]*x[0]+x[5]) - data
            else:
                result = gauss2.gaussian(*x, self.ROI_size)
                self.cache[key] = (result-x[5])/x[0]
                self.counter_calc+=1
                return result - data
        
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
    
        result = self.call_minpack(fun_wrapped, x0, f0, data, None, ftol, xtol, gtol,
                              max_nfev, self.x_scale)
    
        result.message = TERMINATION_MESSAGES[result.status]
        result.success = result.status > 0
    
        return result
    
    def determine_background(self, my_roi):
        roi_background = self.empty_background
        roi_background[0:self.ROI_size] = my_roi[:, 0]
        roi_background[self.ROI_size:self.ROI_size*2] = my_roi[:, -1]
        roi_background[self.ROI_size*2:self.ROI_size*2+self.ROI_size-2] = my_roi[0, 0:-2]
        roi_background[self.ROI_size*2+self.ROI_size-2:] = my_roi[-1, 0:-2]
        
        return np.mean(roi_background)
    
    def phasor_guess(self, data):
        """ Returns an initial guess based on phasor fitting"""
        pos_x, pos_y = self.phasor_fit(data)
        height = data[int(pos_y+0.5), int(pos_x+0.5)]-np.min(data)
        
        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig])
    
    def gaussian(self, height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
        
    def fitgaussian(self, data, peak_index):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
    
        if self.params[peak_index, 0] == 0:
            params = self.phasor_guess(data)
        else:
            params = self.params[peak_index, :]
        p = self.least_squares(params, np.ravel(data))#, gtol=1e-4, ftol=1e-4)
        
        self.params[peak_index, :] = p.x
        
        return [p.x, p.nfev, p.success]
           
    def fitter(self, frame_index, frame, peaks):
        """
        Fits all peaks by scipy least squares fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        peaks : all peaks

        Returns
        -------
        frame_result : fits of each frame

        """

        frame_result = np.zeros([peaks.shape[0], 9])
        
        self.cache._check_size_limit()
        
        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            my_roi_bg = self.determine_background(my_roi)
            my_roi = my_roi - my_roi_bg

            result, its, success = self.fitgaussian(my_roi, peak_index)

            if success == 0:
                continue

            frame_result[peak_index, 0] = frame_index
            frame_result[peak_index, 1] = peak_index
            #start position plus from center in ROI + half for indexing of pixels
            frame_result[peak_index, 2] = result[1]+y-self.ROI_size_1D+0.5
            frame_result[peak_index, 3] = result[2]+x-self.ROI_size_1D+0.5
            frame_result[peak_index, 4] = result[0]
            frame_result[peak_index, 5] = result[3]
            frame_result[peak_index, 6] = result[4]
            frame_result[peak_index, 7] = my_roi_bg
            frame_result[peak_index, 8] = its

        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result
    
    
    def main(self, frames, metadata):
        """
        Main for every fitter method, calls fitter function and returns fits

        Parameters
        ----------
        frames : All frames of .nd2 video
        metadata : Metadata of .nd2 video

        Returns
        -------
        All localizations

        """
        tot_fits = 0
        
        for frame_index, frame in enumerate(frames):
            if frame_index == 0:
                frame_result = self.fitter(frame_index, frame, self.ROI_locations)
                n_fits = frame_result.shape[0]
                width = frame_result.shape[1]
                height = len(frames)*self.ROI_locations.shape[0]
                self.result = np.zeros((height, width))
                self.result[tot_fits:tot_fits+n_fits,:] = frame_result
                tot_fits += n_fits
            else:
                frame_result = self.fitter(frame_index, frame, self.ROI_locations)
                n_fits = frame_result.shape[0]
                self.result[tot_fits:tot_fits+n_fits,:] = frame_result
                tot_fits += n_fits

            if frame_index % (round(metadata['sequence_count']/10,0)) == 0:
                print('Done fitting frame '+str(frame_index)+' of ' + str(metadata['sequence_count']))
            #print('Done fitting frame '+str(frame_index)+' of ' + str(metadata['sequence_count']))
        
        self.result = np.delete(self.result, range(tot_fits,len(frames)*self.ROI_locations.shape[0]), axis=0)

        return self.result
    
#%% Cached scipy last fit guess including background

class scipy_last_fit_guess_background(scipy_last_fit_guess):
    
    def phasor_guess(self, data):
        """ Returns an initial guess based on phasor fitting"""
        pos_x, pos_y = self.phasor_fit(data)
        background = self.determine_background(data)
        height = data[int(pos_y+0.5), int(pos_x+0.5)]-background
        
        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig, background])
    
    def gaussian(self, height, center_x, center_y, width_x, width_y, background):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2) + background
    
    def fitter(self, frame_index, frame, peaks):
        """
        Fits all peaks by scipy least squares fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        peaks : all peaks

        Returns
        -------
        frame_result : fits of each frame

        """

        frame_result = np.zeros([peaks.shape[0], 9])
        
        self.cache._check_size_limit()

        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]

            result, its, success = self.fitgaussian(my_roi, peak_index)

            if success == 0:
                continue

            frame_result[peak_index, 0] = frame_index
            frame_result[peak_index, 1] = peak_index
            #start position plus from center in ROI + half for indexing of pixels
            frame_result[peak_index, 2] = result[1]+y-self.ROI_size_1D+0.5 #y
            frame_result[peak_index, 3] = result[2]+x-self.ROI_size_1D+0.5 #x
            frame_result[peak_index, 4] = result[0]
            frame_result[peak_index, 5] = result[3]
            frame_result[peak_index, 6] = result[4]
            frame_result[peak_index, 7] = result[5]
            frame_result[peak_index, 8] = its

        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result 

#%% Cached scipy phasor guess

class scipy_phasor_fit_guess(scipy_last_fit_guess):
    
    def fitgaussian(self, data, peak_index):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
    
        if self.params[peak_index, 0] == 0:
            params = self.phasor_guess(data)
        else:
            params = self.params[peak_index, :]
            params[2], params[1]= self.phasor_fit(data)
        p = self.least_squares(params, np.ravel(data))#, gtol=1e-4, ftol=1e-4)
        
        self.params[peak_index, :] = p.x
        
        return [p.x, p.nfev, p.success]
    

#%% Phasor for ROI loops

import pyfftw 
from pims import Frame # for converting ND2 generator to one numpy array

class phasor_only_ROI_loop():
    """
    Phasor localizer over ROIs
    """
    def __init__(self, metadata, ROI_size, wavelength, threshold, ROI_locations, METHOD):
        """
        Parameters
        ----------
        metadata : metadata of ND2
        ROI_size : ROI size
        wavelength : laser wavelength
        threshold : # of standard deviations for threhsold
        ROI_locations : found ROIs

        Returns
        -------
        None.

        """

        self.result = []
        self.ROI_size = ROI_size
        self.ROI_size_1D = int((self.ROI_size-1)/2)
        self.ROI_locations = ROI_locations
        self.init_sig = wavelength/(2*metadata['NA']*math.sqrt(8*math.log(2)))/(metadata['calibration_um']*1000)*2 #*2 for better guess
        self.threshold_sigma = threshold
        self.__name__ = METHOD
    
    def phasor_fit_stack(self, frame_stack, roi, roi_index, y, x):
        
        roi_result = np.zeros([frame_stack.shape[0], 6])
        
        roi_bb = pyfftw.empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = pyfftw.empty_aligned((frame_stack.shape[0], 9, 5), dtype='complex128')
        fft_values_list = pyfftw.FFTW(roi_bb, roi_bf,axes=(1,2),flags=('FFTW_MEASURE',), direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)
        
        roi_size = frame_stack.shape[1]
        
        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):
            
            my_frame_bg = np.mean(np.append(np.append(np.append(frame[:, 0],frame[:, -1]), np.transpose(frame[0, 1:-2])), np.transpose(frame[-1, 1:-2])))
        
            if np.max(frame) < my_frame_bg + math.sqrt(my_frame_bg)*self.threshold_sigma:
                continue
            
            ang_x = cmath.phase(fft_values[0, 1])
            if ang_x>0:
                ang_x=ang_x-2*pi
        
            pos_x = abs(ang_x)/(2*pi/roi_size)
        
            ang_y = cmath.phase(fft_values[1,0])
        
            if ang_y >0:
                ang_y = ang_y - 2*pi
        
            pos_y = abs(ang_y)/(2*pi/roi_size)
        
            if pos_x > 8.5:
                pos_x -= roi_size
            if pos_y > 8.5:
                pos_y -= roi_size
                
            roi_result[frame_index, 0] = frame_index
            roi_result[frame_index, 1] = roi_index
            #start position plus from center in ROI + half for indexing of pixels
            roi_result[frame_index, 2] = y+pos_y-self.ROI_size_1D+0.5
            roi_result[frame_index, 3] = x+pos_x-self.ROI_size_1D+0.5
            roi_result[frame_index, 4] = np.max(frame) - my_frame_bg # returns max peak
            roi_result[frame_index, 5] = my_frame_bg # background
    
        roi_result = roi_result[roi_result[:, 2] > 0]    
    
        return roi_result
                
            
    def main(self, frames, metadata):
        """
        Main for phasor over ROI loop, calls fitter function and returns fits.

        Parameters
        ----------
        frames : All frames of .nd2 video
        metadata : Metadata of .nd2 video

        Returns
        -------
        All localizations

        """
        frame_stack_total = Frame(frames)
        
        tot_fits = 0
        created = 0

        for roi_index, roi in enumerate(self.ROI_locations):
            y = int(roi[0])
            x = int(roi[1])
            frame_stack = frame_stack_total[:,y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]

            roi_result = self.phasor_fit_stack(frame_stack, roi, roi_index, y, x)
            
            if created == 0:
                n_fits = roi_result.shape[0]
                width = roi_result.shape[1]
                height = len(frames)*self.ROI_locations.shape[0]
                self.result = np.zeros((height, width))
                created  = 1
                self.result[tot_fits:tot_fits+n_fits,:] = roi_result
                tot_fits += n_fits
            else:
                n_fits = roi_result.shape[0]
                self.result[tot_fits:tot_fits+n_fits,:] = roi_result
                tot_fits += n_fits
                
            if roi_index % (round(self.ROI_locations.shape[0]/10,0)) == 0:
                print('Done fitting ROI '+str(roi_index)+' of ' + str(self.ROI_locations.shape[0]))
            #print('Done fitting ROI '+str(roi_index)+' of ' + str(self.ROI_locations.shape[0]))
            
        self.result = np.delete(self.result, range(tot_fits,len(frames)*self.ROI_locations.shape[0]), axis=0)
                                
        return self.result
    
#%% Dumb phasor ROI loop

class phasor_only_ROI_loop_dumb(phasor_only_ROI_loop):
    
    def phasor_fit_stack(self, frame_stack, roi, roi_index, y, x):
        
        roi_result = np.zeros([frame_stack.shape[0], 4])
        
        roi_bb = pyfftw.empty_aligned(frame_stack.shape, dtype='float64')
        roi_bf = pyfftw.empty_aligned((frame_stack.shape[0], 9, 5), dtype='complex128')
        fft_values_list = pyfftw.FFTW(roi_bb, roi_bf,axes=(1,2),flags=('FFTW_MEASURE',), direction='FFTW_FORWARD')
        roi_bb = frame_stack
        fft_values_list = fft_values_list(roi_bb)
        
        roi_size = frame_stack.shape[1]
        
        for frame_index, (fft_values, frame) in enumerate(zip(fft_values_list, frame_stack)):
                        
            ang_x = cmath.phase(fft_values[0, 1])
            if ang_x>0:
                ang_x=ang_x-2*pi
        
            pos_x = abs(ang_x)/(2*pi/roi_size)
        
            ang_y = cmath.phase(fft_values[1,0])
        
            if ang_y >0:
                ang_y = ang_y - 2*pi
        
            pos_y = abs(ang_y)/(2*pi/roi_size)
        
            if pos_x > 8.5:
                pos_x -= roi_size
            if pos_y > 8.5:
                pos_y -= roi_size
                
            roi_result[frame_index, 0] = frame_index
            roi_result[frame_index, 1] = roi_index
            #start position plus from center in ROI + half for indexing of pixels
            roi_result[frame_index, 2] = y+pos_y-self.ROI_size_1D+0.5
            roi_result[frame_index, 3] = x+pos_x-self.ROI_size_1D+0.5
    
        return roi_result  
    