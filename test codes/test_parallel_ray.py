# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 00:38:19 2020

@author: s150127
"""
# import sys
import numpy as np
from pims import ND2_Reader  # reader of ND2 files
from pims import Frame  # for converting ND2 generator to one numpy array
import time
import multiprocessing as mp

import sys

sys.path.append('../')
# sys.path.insert(1, 'C:\Users\s150127\OneDrive - TU Eindhoven\Werk\MBx Python analysis\Code')

import gaussian_fitting
import tools
import cProfile
import ray

pr = cProfile.Profile()

# %% Init
ROI_SIZE = 9
WAVELENGTH = 637  # nm
THRESHOLD = 5  # X*sigma
filenames = ("C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2",)

n_processes = int(mp.cpu_count() / 4)


# %% Main

def main(fitter, frames, metadata, q):
    q.put(fitter.main(frames, metadata))


# %% Main
if __name__ == '__main__':

    pr.enable()

    ROI_locations = np.load('ROI_locations.npy')
    ROI_locations = ROI_locations - 1
    ROI_locations = ROI_locations[0:8, :]

    ## switch array columns since MATLAB gives x,y. Python likes y,x
    ROI_locations = tools.switch(ROI_locations)

    for name in filenames:
        with ND2_Reader(name) as ND2:

            ## parse ND2 info
            frames = ND2
            metadata = ND2.metadata
            frames = frames[0:20]

            start = time.time()
            processes = []
            q = mp.Queue()
            for i in range(0, n_processes):
                roi_locations_split = np.array_split(ROI_locations, n_processes)
                scipy_last_fit_guess_roi_loop = gaussian_fitting.scipy_last_fit_guess_roi_loop(metadata, ROI_SIZE,
                                                                                               WAVELENGTH, THRESHOLD,
                                                                                               roi_locations_split[i])
                processes.append(mp.Process(target=main, args=(scipy_last_fit_guess_roi_loop, frames, metadata, q)))

            for p in processes:
                p.start()

            for p in processes:
                p.join()

            result = [q.get() for p in processes]

            first = True
            for res in result:
                if first == True:
                    results = res
                    first = False
                else:
                    results = np.vstack((results, res))
            print('Time taken: ' + str(round(time.time() - start, 3)) + ' s. Fits done: ' + str(results.shape[0]))
            pr.disable()
            pr.print_stats(sort='time')
