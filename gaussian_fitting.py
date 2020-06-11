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

"""
#%% Generic imports
import math
import numpy as np

#%% Main

class main_localizer():
    """
    class that everything else inherits from, has ROI finders in it.
    """

    def __init__(self, metadata, ROI_size, wavelength, threshold, ROI_locations):
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
        self.init_sig = wavelength/(2*metadata['NA']*math.sqrt(8*math.log(2)))/(metadata['calibration_um']*1000)
        self.threshold_sigma = threshold

    def find_peaks_v2(self, frame):
        """
        Find peaks in a single frame

        Parameters
        ----------
        frame : Frame

        Returns
        -------
        peaks : All peaks in frame

        """
        peaks = []
        for ROI_index, ROI in enumerate(self.ROI_locations):
            y = int(ROI[0])
            x = int(ROI[1])
            myROI = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]
            myROI_bg = np.mean(np.append(np.append(np.append(myROI[:, 0], myROI[:, -1]), np.transpose(myROI[0, 1:-2])), np.transpose(myROI[-1, 1:-2])))
            threshold = myROI_bg + math.sqrt(myROI_bg)*self.threshold_sigma

            maximum = np.max(myROI)
            if maximum < threshold:
                continue
            indices = np.where(myROI == maximum)

            indices = np.asarray([x for x in zip(indices[0], indices[1])])
            indices = indices[0, :]
            indices[0] = indices[0]-self.ROI_size_1D+y
            indices[1] = indices[1]-self.ROI_size_1D+x

            peak = np.append(indices, [maximum, ROI_index])
            if peaks == []:
                peaks = peak
            else:
                try:
                    peaks = np.vstack((peaks, peak))
                except:
                    pass

        return peaks

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
        for frame_index, frame in enumerate(frames):
            peaks = self.find_peaks_v2(frame)
            frame_result = self.fitter(frame_index, frame, peaks)
            if frame_index == 0:
                self.result = frame_result
            else:
                self.result = np.vstack((self.result, frame_result))

            print('Done fitting frame '+str(frame_index)+' of ' + str(metadata['sequence_count']))

        return self.result


#%% rainSTORM


class rainSTORM_Dion(main_localizer):
    """
    rainSTORM class

    """

    def __init__(self, metadata, ROI_size, wavelength, threshold, ROI_locations):
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
        super().__init__(metadata, ROI_size, wavelength, threshold, ROI_locations)
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
                frame_result[peak_index, 1] = peak[3]
                frame_result[peak_index, 2] = fit_y
                frame_result[peak_index, 3] = fit_x
                frame_result[peak_index, 4] = peak[2]
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

            total = np.sum(myROI)

            frame_result[peak_index, 0] = frame_index
            frame_result[peak_index, 1] = peak[3]
            frame_result[peak_index, 2] = y
            frame_result[peak_index, 3] = x
            frame_result[peak_index, 4] = peak[2]
            frame_result[peak_index, 5] = total


        return frame_result


#%% Phasor

from scipy.fftpack import ifftn
from math import pi
import cmath

class phasor_only(main_localizer):
    """
    Max Bergkamp inspired phasor fitting

    Returns
    -------
    frame_result : fits of each frame

    """
    def fitter(self, frame_index, frame, peaks):
        """
        Fits all peaks by phasor fitting

        Parameters
        ----------
        frame_index : Index of frame
        frame : frame
        peaks : all peaks

        Returns
        -------
        frame_result : fits of each frame

        """

        frame_result = np.zeros([peaks.shape[0], 5])

        for peak_index, peak in enumerate(peaks):

            y = int(peak[0])
            x = int(peak[1])

            my_roi = frame[y-self.ROI_size_1D:y+self.ROI_size_1D+1, x-self.ROI_size_1D:x+self.ROI_size_1D+1]

            fft_values = ifftn(my_roi)

            roi_size = float(self.ROI_size)
            ang_x = cmath.phase(fft_values[0, 1])
            if ang_x>0:
                ang_x=ang_x-2*pi

            pos_x = roi_size - abs(ang_x)/(2*pi/roi_size)

            ang_y = cmath.phase(fft_values[1,0])

            if ang_y >0:
                ang_y = ang_y - 2*pi

            pos_y = roi_size - abs(ang_y)/(2*pi/roi_size)

            if pos_x > 8.5:
                pos_x -= roi_size
            if pos_y > 8.5:
                pos_y -= roi_size

            frame_result[peak_index, 0] = frame_index
            frame_result[peak_index, 1] = peak[3]
            #start position plus from center in ROI + half for indexing of pixels
            frame_result[peak_index, 2] = y+pos_y-self.ROI_size_1D+0.5
            frame_result[peak_index, 3] = x+pos_x-self.ROI_size_1D+0.5
            frame_result[peak_index, 4] = peak[2]

        return frame_result

#%% Fully build-in option

from scipy import optimize

class scipy_least_squares(main_localizer):

    """
    Build-in scipy least squares fitting, using centroid fitting as initial guess

    Returns
    -------
    frame_result : fits of each frame
    """

    def gaussian(self, height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

    def moments(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution by calculating its
        moments """
        success = 1
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        if y > 9 or x > 9 or x < 0 or y < 0:
            success = 0
            return 0, 0, 0, 0, 0, success

        col = data[:, int(y)]
        width_x = self.init_sig
        row = data[int(x), :]
        width_y = self.init_sig
        height = data.max()
        return height, x, y, width_x, width_y, success

    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = self.moments(data)
        success = params[-1]
        params = params[0:-1]
        if success == 0:
            return [0, 0, success]
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
                                     data)
        p = optimize.least_squares(errorfunction, params)

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

        frame_result = np.zeros([peaks.shape[0], 8])

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
            frame_result[peak_index, 1] = peak[3]
            #start position plus from center in ROI + half for indexing of pixels
            frame_result[peak_index, 2] = result[1]+y-self.ROI_size_1D+0.5
            frame_result[peak_index, 3] = result[2]+x-self.ROI_size_1D+0.5
            frame_result[peak_index, 4] = result[0]
            frame_result[peak_index, 5] = result[3]
            frame_result[peak_index, 6] = result[4]
            frame_result[peak_index, 7] = its

        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result

#%% Build in using phasor as first estimate

class scipy_phasor_guess(scipy_least_squares):
    """
    Build-in scipy least squares fitting, using phasor fitting as initial guess
    """

    def phasor_guess(self, data):
        """
        Parameters
        ----------
        data : input ROI

        Returns
        -------
        height of gaussian, x position, y position, sigma_x and sigma_y
        """
        fft_values = ifftn(data)

        roi_size = self.ROI_size
        ang_x = cmath.phase(fft_values[0, 1])
        if ang_x>0:
            ang_x=ang_x-2*pi

        pos_x = roi_size - abs(ang_x)/(2*pi/roi_size)

        ang_y = cmath.phase(fft_values[1,0])

        if ang_y >0:
            ang_y = ang_y - 2*pi

        pos_y = roi_size - abs(ang_y)/(2*pi/roi_size)

        if pos_x > 8.5:
            pos_x -= roi_size
        if pos_y > 8.5:
            pos_y -= roi_size


        height = data[int(pos_y+0.5), int(pos_x+0.5)]-np.min(data)

        return height, pos_y, pos_x, self.init_sig, self.init_sig


    def fitgaussian(self, data):
        params = self.phasor_guess(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
         data)
        p = optimize.least_squares(errorfunction, params)

        return [p.x, p.nfev, p.success]

#%% ROI loops

from pims import Frame # for converting ND2 generator to one numpy array

class scipy_phashor_guess_roi_loop(scipy_phasor_guess):
    """
    Build-in scipy least squares fitting, using phasor fitting as initial guess, now looping first of ROIs, then frames
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
                if self.result == []:
                    self.result = frame_result
                else:
                    self.result = np.vstack((self.result, frame_result))

            print('Done fitting ROI '+str(roi_index)+' of ' + str(self.ROI_locations.shape[0]))

        return self.result

    def fitter(self, frame_index, my_roi, roi, roi_index, params):
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

        frame_result = np.zeros([1, 8])

        y = int(roi[0])
        x = int(roi[1])

        result, its, success = self.fitgaussian(my_roi, params)

        if success == 0:
            return []

        frame_result[0, 0] = frame_index
        frame_result[0, 1] = roi_index
        #start position plus from center in ROI + half for indexing of pixels
        frame_result[0, 2] = result[1]+y-self.ROI_size_1D+0.5
        frame_result[0, 3] = result[2]+x-self.ROI_size_1D+0.5
        frame_result[0, 4] = result[0]
        frame_result[0, 5] = result[3]
        frame_result[0, 6] = result[4]
        frame_result[0, 7] = its

        frame_result = frame_result[frame_result[:, 4] > 0]

        return frame_result

#%% Fitter based on last fit

class scipy_last_fit_guess_roi_loop(scipy_phashor_guess_roi_loop):
    """
    Build-in scipy least squares fitting, using last fit as initial guess, now looping first of ROIs, then frames
    """
    def fitgaussian(self, data, params):
        """
        Fits gausian to data and returns fit parameters, number of its and success parameter
        """
        if self.params == []:
            params = self.phasor_guess(data)
        else:
            params = self.params
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
         data)
        p = optimize.least_squares(errorfunction, params)

        self.params = p.x

        return [p.x, p.nfev, p.success]

#%% Parallel proccessed version
import time

class scipy_last_fit_guess_roi_loop_parallel(scipy_last_fit_guess_roi_loop):
    """
    Build-in scipy least squares fitting, using last fit as initial guess, now looping first of ROIs, then frames, parallel processed
    """

    def main(self, frames, metadata, q):
        q.put(super().main(frames, metadata))
        time.sleep(0.1)