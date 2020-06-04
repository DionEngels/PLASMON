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
#%% Imports
import math
import numpy as np

#%% Code

def doastorm_3D():
    """
    3D Daostorm inspired fitting

    https://github.com/ZhuangLab/storm-analysis/tree/master/storm_analysis/daostorm_3d

    Returns
    -------
    None.

    """
    pass

class rainSTORM_Dion():
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
        self.result = []
        self.ROI_size = ROI_size
        self.ROI_size_1D = int((self.ROI_size-1)/2)
        self.allow_sig = [0.25, self.ROI_size_1D+1]
        self.max_its = 60
        self.maxtol = 0.2
        self.mintol = 0.06
        self.allow_x = 1
        self.init_sig = wavelength/(2*metadata['NA']*math.sqrt(8*math.log(2)))/(metadata['calibration_um']*1000)
        self.threshold_sigma = threshold
        self.ROI_locations = ROI_locations


    def main(self, frames, metadata):
        """
        Main for rainSTORM, calls other functions and returns all fits

        Parameters
        ----------
        frames : All frames of .nd2 video
        metadata : Metadata of .nd2 video

        Returns
        -------
        All localizations

        """
        for frame_index, frame in enumerate(frames):
            frame_result = self.main_loop(frame, frame_index)
            if frame_index == 0:
                self.result = frame_result
            else:
                self.result = np.vstack((self.result, frame_result))

            print('Done fitting frame '+str(frame_index)+' of ' + str(metadata['sequence_count']))

        return self.result


    def main_loop(self, frame, frame_index):
        """
        Loop for every frame

        Parameters
        ----------
        frame : The frame to be fitted
        frame_index : Index of frame

        Returns
        -------
        Fits of this frame

        """
        peaks = self.find_peaks_v2(frame)
        return self.fitter(frame_index, frame, peaks)


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
        for ROI in self.ROI_locations:
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

            peak = np.append(indices, [maximum])
            if peaks == []:
                peaks = peak
            else:
                try:
                    peaks = np.vstack((peaks, peak))
                except:
                    pass

        return peaks

    def fitter(self, frame_index, frame, peaks):
        """
        Fits all peaks

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
                fit_y = float(y) + y0 + 0.5
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
                    fit_x = float(x) + x0 + 0.5
                    flag_x = True


            if flag_y and flag_x:
                frame_result[peak_index, 0] = frame_index
                frame_result[peak_index, 1] = peak_index
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





def max_bergkamp():
    """
    Max Bergkamp inspired fitting

    Returns
    -------
    None.

    """
    pass
