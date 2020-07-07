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
"""
#%% Generic imports
from __future__ import division, print_function, absolute_import
import math
import numpy as np

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
    
#%% Main

class main_localizer():
    """
    class that everything else inherits from.
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
    
#%% rainSTORM

class rainSTORM_Dion(main_localizer):
    """
    rainSTORM class
    """
    def __init__(self, metadata, ROI_size, wavelength, threshold, ROI_locations, METHOD):
        """
        Init of rainSTORM

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
        super().__init__(metadata, ROI_size, wavelength, threshold, ROI_locations, METHOD)
        self.allow_sig = [0.25, self.ROI_size_1D+1]
        self.max_its = 60
        self.maxtol = 0.2
        self.mintol = 0.06
        self.allow_x = 1

    def fitter(self, frame_index, frame, peaks):
        """
        Fits all peaks with 2D gaussian, iterative solver self-build

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        peaks : all peaks

        Returns
        -------
        frame_result : fits of each frame

        """
        init_x0 = 0
        n_fails = 0

        frame_result = np.zeros([peaks.shape[0], 12])

        for peak_index, peak in enumerate(peaks):
            y = int(peak[0])
            x = int(peak[1])
            myROI = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            #myROI_min = np.amin(myROI)
            myROI_bg = np.mean(np.append(np.append(np.append(myROI[:, 0],myROI[:, -1]), np.transpose(myROI[0, 1:-2])), np.transpose(myROI[-1, 1:-2])))
            myROI = myROI - myROI_bg
            
            if np.max(myROI) < math.sqrt(myROI_bg)*self.threshold_sigma:
                n_fails+= 1
                continue
            
            flag_y = False
            flag_x = False

            xx = np.transpose(range(-self.ROI_size_1D, self.ROI_size_1D+1))
            yy = np.transpose(range(-self.ROI_size_1D, self.ROI_size_1D+1))
            xSum = np.sum(myROI, axis=1)
            ySum = np.transpose(np.sum(myROI, axis=0))

            y0 = init_x0
            sig_y = self.init_sig
            C = xSum[self.ROI_size_1D]
            for i in range(0, self.max_its):
                fit = C*np.exp(-np.square(yy-y0)/(2*sig_y**2))
                beta = xSum - fit
                A = np.vstack((fit/C, fit* (yy-y0/sig_y**2), fit*(np.square(yy-y0)/sig_y**3)))
                b = np.matmul(A, beta)
                a = np.matmul(A, np.transpose(A))

                if abs(y0) > self.allow_x or sig_y < self.allow_sig[0] or sig_y > self.allow_sig[1]:
                    residue_y = 0
                    break

                residue_y = np.sum(np.square(beta))/np.sum(np.square(xSum))
                if residue_y < self.mintol and abs(y0) < self.allow_x and sig_y > self.allow_sig[0] and sig_y < self.allow_sig[1]:
                    break

                try:
                    dL = np.matmul(a, 1/b)
                    C = C+dL[0]
                    y_change = dL[1]
                    if abs(y_change) > 0.1:
                        y_change = y_change/abs(y_change)*0.1
                    y0 = y0+y_change
                    sig_change = dL[2]
                    if abs(sig_change) > 0.1:
                        sig_change = sig_change/abs(sig_change)*0.1
                    sig_y = sig_y +sig_change
                except:
                    residue_y = 0
                    break

            if residue_y < self.maxtol and abs(y0) < self.allow_x and sig_y > self.allow_sig[0] and sig_y < self.allow_sig[1]:
                fit_y = float(y) + y0 + 0.5 # add half for indexing of pixels
                flag_y = True

            if flag_y:
                x0 = init_x0
                sig_x = sig_y
                C = ySum[self.ROI_size_1D]
                for j in range(0, self.max_its):
                    fit = C*np.exp(-np.square(xx-x0)/(2*sig_x**2))
                    beta = ySum - fit
                    A = np.vstack((fit/C, fit*(xx-x0/sig_x**2), fit*(np.square(xx-x0)/sig_x**3)))
                    b = np.matmul(A, beta)
                    a = np.matmul(A, np.transpose(A))

                    if abs(x0) > self.allow_x or sig_x < self.allow_sig[0] or sig_x > self.allow_sig[1]:
                        residue_x = 0
                        break

                    residue_x = np.sum(np.square(beta))/np.sum(np.square(ySum))
                    if residue_x < self.mintol and abs(x0) < self.allow_x and sig_x > self.allow_sig[0] and sig_x < self.allow_sig[1]:
                        break

                    try:
                        dL = np.matmul(a, 1/b)
                        C = C+dL[0]/50
                        x0 = x0+dL[1]/50
                        sig_x = sig_x +dL[2]/50
                    except:
                        residue_x = 0
                        break

                if residue_x < self.maxtol and abs(x0) < self.allow_x and sig_x > self.allow_sig[0] and sig_x < self.allow_sig[1]:
                    fit_x = float(x) + x0 + 0.5 # add half for indexing of pixels
                    flag_x = True


            if flag_y and flag_x:
                frame_result[peak_index, 0] = frame_index
                frame_result[peak_index, 1] = peak_index
                frame_result[peak_index, 2] = fit_y
                frame_result[peak_index, 3] = fit_x
                frame_result[peak_index, 4] = np.max(myROI)
                frame_result[peak_index, 5] = sig_y
                frame_result[peak_index, 6] = sig_x
                frame_result[peak_index, 7] = (residue_y + residue_x) / 2
                frame_result[peak_index, 8] = residue_y
                frame_result[peak_index, 9] = residue_x
                frame_result[peak_index, 10] = i
                frame_result[peak_index, 11] = j

            else:
                n_fails += 1


        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result

#%% Summation

class summation(main_localizer):
    """
    Summation fitter

    """

    def fitter(self, frame_index, frame, peaks):
        """
        Fits all peaks by summing them

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        peaks : all peaks

        Returns
        -------
        frame_result : fits of each frame

        """

        frame_result = np.zeros([peaks.shape[0], 6])

        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            myROI = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            myROI_bg = np.mean(np.append(np.append(np.append(myROI[:, 0],myROI[:, -1]), np.transpose(myROI[0, 1:-2])), np.transpose(myROI[-1, 1:-2])))
            
            if np.max(myROI) < myROI_bg + math.sqrt(myROI_bg)*self.threshold_sigma:
                continue

            total = np.sum(myROI)

            frame_result[peak_index, 0] = frame_index
            frame_result[peak_index, 1] = peak_index
            frame_result[peak_index, 2] = y
            frame_result[peak_index, 3] = x
            frame_result[peak_index, 4] = np.max(myROI)-myROI_bg
            frame_result[peak_index, 5] = total


        frame_result = frame_result[frame_result[:, 4] > 0]
        
        return frame_result

#%% Build in using phasor as first estimate

from scipy import optimize

class scipy_phasor(main_localizer, base_phasor):

    """
    Build-in scipy least squares fitting, using phasor fitting as initial guess
    """

    def gaussian(self, height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
    
    def phasor_guess(self, data):
        """ Returns an initial guess based on phasor fitting"""
        pos_x, pos_y = self.phasor_fit(data)
        height = data[int(pos_y+0.5), int(pos_x+0.5)]-np.min(data)
        
        return height, pos_y, pos_x, self.init_sig, self.init_sig
        

    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = self.phasor_guess(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
         data)
        
        test_method ='lm' # different method test
        # test_tol =0.1 # tolerance test
        #test_bounds = ([-np.inf, params[1]*0.99, params[2]*0.99, -np.inf, -np.inf], [np.inf, params[1]*1.01, params[2]*1.01, np.inf, np.inf])
        
        
        p = optimize.least_squares(errorfunction, params, method=test_method)#, bounds=test_bounds) 

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

        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            my_roi_bg = np.mean(np.append(np.append(np.append(my_roi[:, 0],my_roi[:, -1]), np.transpose(my_roi[0, 1:-2])), np.transpose(my_roi[-1, 1:-2])))
            my_roi = my_roi - my_roi_bg

            result, its, success = self.fitgaussian(my_roi)

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

#%% Scipy phasor guess log scale 

# currently not used

class scipy_phasor_log(scipy_phasor):
    """
    Build-in scipy least squares fitting on log scale, using phasor fitting as initial guess
    """
    
    def gaussian(self, height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2
    
    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        pos_x, pos_y = self.phasor_fit(data)
        
        data[data< 1] = 1
        data = np.log(data)
        
        height = data[int(pos_y+0.5), int(pos_x+0.5)]-np.min(data)
        
        params = height, pos_y, pos_x, self.init_sig, self.init_sig
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
         data)
        p = optimize.least_squares(errorfunction, params)#, method='lm')

        return [p.x, p.nfev, p.success]
        
        
#%% Scipy last fit guess 
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

class scipy_last_fit_guess(scipy_phasor):
    """
    Build-in scipy least squares fitting, with last fit as initial guess
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
        super().__init__(metadata, ROI_size, wavelength, threshold, ROI_locations, METHOD)
        self.indices = np.indices((ROI_size, ROI_size))
        self.cache = {}  
        self.x_scale = np.ones(5)
        self.active_mask = np.zeros_like(self.x_scale, dtype=int)
        
        self.lb = np.resize(-np.inf, self.x_scale.shape)
        self.ub = np.resize(np.inf, self.x_scale.shape)
        self.use_one_sided = np.resize(False, self.x_scale.shape)
        self.empty_background = np.zeros(self.ROI_size*2+(self.ROI_size-2)*2, dtype=np.uint16)
        
        
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
        
    def approx_derivative(self, fun, x0, f0=None):
         
        if x0.ndim > 1:
            raise ValueError("`x0` must have at most 1 dimension.")
    
        lb = self.lb
        ub = self.ub
    
        if lb.shape != x0.shape or ub.shape != x0.shape:
            raise ValueError("Inconsistent shapes between bounds and `x0`.")
    
        def fun_wrapped(x):
            f = np.atleast_1d(fun(x))
            if f.ndim > 1:
                raise RuntimeError("`fun` return value has "
                                   "more than 1 dimension.")
            return f
       
        h = self.compute_absolute_step(x0)

        return self.dense_difference(fun_wrapped, x0, f0, h)
                
    def call_minpack(self, fun, x0, f0, jac, ftol, xtol, gtol, max_nfev, diag):

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
        
        j = self.approx_derivative(fun, x, f0=f0)
        
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
    
    def least_squares(self, fun, x0, ftol=1e-8, xtol=1e-8, gtol=1e-8, 
        max_nfev=None, args=(), kwargs={}):
        
        def fun_wrapped(x):
            
            key = tuple(x) 
            if key in self.cache:
                return self.cache[key]
            else:
                result = fun(x, *args, **kwargs)
                self.cache[key] = result
                return result
        
        if max_nfev is not None and max_nfev <= 0:
            raise ValueError("`max_nfev` must be None or positive integer.")
    
        if np.iscomplexobj(x0):
            raise ValueError("`x0` must be real.")
        
        if x0.ndim > 1:
            raise ValueError("`x0` must have at most 1 dimension.")
         
        f0 = fun_wrapped(x0)
    
        if f0.ndim != 1:
            raise ValueError("`fun` must return at most 1-d array_like.")
    
        if not np.all(np.isfinite(f0)):
            raise ValueError("Residuals are not finite in the initial point.")
    
        result = self.call_minpack(fun_wrapped, x0, f0, None, ftol, xtol, gtol,
                              max_nfev, self.x_scale)
    
        result.message = TERMINATION_MESSAGES[result.status]
        result.success = result.status > 0
    
        return result
    
    def phasor_guess(self, data):
        """ Returns an initial guess based on phasor fitting"""
        pos_x, pos_y = self.phasor_fit(data)
        height = data[int(pos_y+0.5), int(pos_x+0.5)]-np.min(data)
        
        return np.array([height, pos_y, pos_x, self.init_sig, self.init_sig])
    
    def fitgaussian(self, data, peak_index):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
    
        if np.sum(self.params[peak_index, :]) == 0:
            params = self.phasor_guess(data)
        else:
            params = self.params[peak_index, :]
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*self.indices) -
        data)  
        p = self.least_squares(errorfunction, params)#, gtol=1e-4, ftol=1e-4)
        
        self.params[peak_index, :] = p.x
        
        return [p.x, p.nfev, p.success]
    
    
    def determine_background(self, my_roi):
        roi_background = self.empty_background
        roi_background[0:self.ROI_size] = my_roi[:, 0]
        roi_background[self.ROI_size:self.ROI_size*2] = my_roi[:, -1]
        roi_background[self.ROI_size*2:self.ROI_size*2+self.ROI_size-2] = np.transpose(my_roi[0, 0:-2])
        roi_background[self.ROI_size*2+self.ROI_size-2:] = np.transpose(my_roi[-1, 0:-2])
        
        return np.mean(roi_background)
        
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

        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            my_roi_bg = self.determine_background(my_roi)
            #my_roi_bg = np.mean(np.append(np.append(np.append(my_roi[:, 0], 
            #my_roi[:, -1]), np.transpose(my_roi[0, 1:-2])), np.transpose(my_roi[-1, 1:-2])))
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
        
        self.params = np.zeros((self.ROI_locations.shape[0],5))
        
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

#%% ROI loops

from pims import Frame # for converting ND2 generator to one numpy array

class main_localizer_roi_loop(main_localizer):
    """
    Class that every looping over ROI fitter inherits from
    """
    
    def main(self, frames, metadata):
        """
        Main for every fitter method, calls fitter function and returns fits. Loops over ROIs

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
            self.params = []

            for frame_index, frame in enumerate(frame_stack):
                my_roi_bg = np.mean(np.append(np.append(np.append(frame[:, 0],frame[:, -1]), np.transpose(frame[0, 1:-2])), np.transpose(frame[-1, 1:-2])))
                my_roi = frame - my_roi_bg

                threshold = math.sqrt(my_roi_bg)*self.threshold_sigma
                maximum = np.max(my_roi)

                if maximum < threshold:
                    continue
                frame_result = self.fitter(frame_index, my_roi, roi, roi_index, [])
                
                if frame_result == []:
                    continue
                
                if created == 0:
                    n_fits = frame_result.shape[0]
                    width = frame_result.shape[1]
                    height = len(frames)*self.ROI_locations.shape[0]
                    self.result = np.zeros((height, width))
                    created = 1
                    self.result[tot_fits:tot_fits+n_fits,:] = frame_result
                    tot_fits += n_fits
                else:
                    n_fits = frame_result.shape[0]
                    self.result[tot_fits:tot_fits+n_fits,:] = frame_result
                    tot_fits += n_fits

            print('Done fitting ROI '+str(roi_index)+' of ' + str(self.ROI_locations.shape[0]))
            
        self.result = np.delete(self.result, range(tot_fits,len(frames)*self.ROI_locations.shape[0]), axis=0)
                                
        return self.result
    
#%% Phasor for ROI loops

import pyfftw 

class phasor_only_ROI_loop(main_localizer):
    """
    Phasor localizer over ROIs
    """
    
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
    
#%% Dion fitter method
from numpy import inner, max, diag, eye, Inf
from numpy.linalg import solve

class dions_fitter(base_phasor, main_localizer):
    
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
        self.init_sig = wavelength/(2*metadata['NA']*math.sqrt(8*math.log(2)))/(metadata['calibration_um']*1000) #*2 for better guess
        self.threshold_sigma = threshold
        self.__name__ = METHOD
        
        self.sigmas = np.zeros((self.ROI_locations.shape[0]))
        self.max_its = 60
        self.maxtol = 0.2
        self.mintol = 0.06
    
    def gaussian(self, height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
    
    def fJ(self, pars, x, y = 0):
        "Calculation of function and Jacobian for one-dimensional Gaussian."
        A, m, s, offs = pars[0:4]
        f = A*np.exp( - (x-m)**2 / (2*s**2))
        if 1:
            J = np.empty(shape = (4,)+x.shape, dtype = np.float_)
            J[0] = 1.0/A * f
            J[1] = f*(x-m)/s**2
            J[2] = f*(x-m)**2/s**3
            J[3] = 1
            return f + (offs - y), J
        return f + (offs - y) 
    
    def LM(self, fun, pars, args, 
           tau = 1e-2, eps1 = 1e-6, eps2 = 1e-6, kmax = 20):
        """Implementation of the Levenberg-Marquardt algorithm in pure
        Python. Solves the normal equations."""
        
        p = pars
        f, J = fun(p, *args)
    
        A = inner(J,J)
        g = inner(J,f)
    
        I = eye(len(p))
        
        success = True
    
        k = 0; nu = 2
        mu = tau * max(diag(A))
        stop = norm(g, Inf) < eps1
        while not stop and k < kmax:
            k += 1
    
            try:
                d = solve( A + mu*I, -g)
            except np.linalg.LinAlgError:
                stop = True
                success = False
                break
    
            if norm(d) < eps2*(norm(p) + eps2):
                stop = True
                break
    
            pnew = p + d
    
            fnew, Jnew = fun(pnew, *args)
            rho = (norm(f)**2 - norm(fnew)**2)/inner(d, mu*d - g)
            
            if rho > 0:
                p = pnew
                A = inner(Jnew, Jnew)
                g = inner(Jnew, fnew)
                f = fnew
                J = Jnew
                if (norm(g, Inf) < eps1):
                    stop = True
                    break
                mu = mu * max([1.0/3, 1.0 - (2*rho - 1)**3])
                nu = 2.0
            else:
                mu = mu * nu
                nu = 2*nu
    
        else:
            success = False
        
        return p, k, success
        
    def find_sigma(self, roi, pos_x, pos_y, init_sigma):
            
        xx = np.transpose(range(0, self.ROI_size_1D*2+1))
        yy = np.transpose(range(0, self.ROI_size_1D*2+1))
        x_sum = np.sum(roi, axis=1)
        y_sum = np.transpose(np.sum(roi, axis=0))
        
        startpars_x = [np.max(x_sum), pos_x, init_sigma, np.min(x_sum)]
        startpars_y = [np.max(y_sum), pos_y, init_sigma, np.min(y_sum)]
        
        result_x, its_x, success_x = self.LM(self.fJ, startpars_x, (xx, x_sum), tau = 1e-4)
        
        result_y, its_y, success_y = self.LM(self.fJ, startpars_y, (yy, y_sum), tau = 1e-4)
        
        return (result_y[2]+result_x[2])/2, its_x+its_y, success_y and success_y
        
        
    
    def fit(self, roi, peak_index):
        
        pos_x, pos_y = self.phasor_fit(roi)
        
        if self.sigmas[peak_index] == 0:
            init_sigma = self.init_sig
        else:
            init_sigma = self.sigmas[peak_index]
        
        sigma, its, success = self.find_sigma(roi, pos_x, pos_y, init_sigma)
        
        if success == False:
            return 0, 0, success
        
        self.sigmas[peak_index] = sigma
        
        height = np.max(roi)
                
        roi_fit = self.gaussian(height, pos_x, pos_y, sigma, sigma)(*np.indices(roi.shape))
        
        factor = height/np.max(roi_fit)
        res_height = factor*height
        
        result = res_height, pos_y, pos_x, sigma, sigma
        
        return result, its, success
    
    def fitter(self, frame_index, frame, peaks):
        frame_result = np.zeros([peaks.shape[0], 9])

        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            my_roi_bg = np.mean(np.append(np.append(np.append(my_roi[:, 0],my_roi[:, -1]), np.transpose(my_roi[0, 1:-2])), np.transpose(my_roi[-1, 1:-2])))
            my_roi = my_roi - my_roi_bg

            result, its, success = self.fit(my_roi, peak_index)

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
    