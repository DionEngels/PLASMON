# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 00:38:19 2020

@author: Dion Engels
MBx Python Data Analysis

Multiprocessing main

----------------------------

v1.0: ROI split: 15/07/2020
v1.1: Bugfix saving correct info
v1.2: Frame split: 16/07/2020
v1.3: cleanup
v1.4: final for now, checked frame splitting

"""

import numpy as np
from pims import ND2_Reader # reader of ND2 files
import time
import multiprocessing as mp
import os

import _code.fitters as fitting
import _code.roi_finding as roi_finding
import _code.tools as tools
#import cProfile

#pr = cProfile.Profile()

#%% Init
ROI_SIZE = 7
WAVELENGTH = 637 #nm
filenames = ("C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2",)

n_processes = int(mp.cpu_count())

#%% Main

def split_frames(length, n_processes):

    split_result = [None]*n_processes
    start = 0
    for j in range(n_processes):
        if j+1 == n_processes:
            arr_length = length-start
        else:
            arr_length = int(length/n_processes*(1+0.05*(n_processes/2-j)))
        split_result[j] = np.arange(start, start+arr_length)
        start += arr_length

    return split_result


def main(name, fitter, frames_split, roi_locations, shared):
    
    ND2 = ND2_Reader(name)
    frames = ND2
    metadata = ND2.metadata
    metadata['sequence_count'] = len(frames_split)
    frames = frames[frames_split]

    local_result = fitter.main(frames, metadata, roi_locations)
    
    for result_index, result in enumerate(local_result):
        #print(result)
        shared[9*result_index:9*(result_index+1)] = result[:]
        # for index, value in enumerate(result):
        #     shared[9*result_index+index] = value

    print("Done fiting")

#%% Main
if __name__ == '__main__':

    #pr.enable()

    for name in filenames:
        with ND2_Reader(name) as ND2:

            basedir = os.getcwd()
            directory = name.split(".")[0].split("/")[-1]
            path = os.path.join(basedir, directory)
            try:
                os.mkdir(path)
            except:
                pass

            ## parse ND2 info
            frames = ND2
            metadata = ND2.metadata
            #frames = frames[0:2]
            
            roi_finder = roi_finding.roi_finder(ROI_SIZE, frames[0])#, intensity_min = 800)
            fitter = fitting.scipy_last_fit_guess(metadata, ROI_SIZE,
                                                  WAVELENGTH, roi_finder.intensity_min, 
                                                  "ScipyLastFitGuess", 5)
            roi_locations = roi_finder.main(frames[0], fitter)
            
            #roi_locations = roi_locations[0:4, :]
            
            #fitter = fitting.phasor_only_ROI_loop(metadata, ROI_SIZE, WAVELENGTH, roi_finder.intensity_min, "Phasor")
            
            start = time.time()

            shared = [None]*n_processes
            result = np.zeros((len(roi_locations)*len(frames), 9))
            
            frames_split = np.array_split(list(range(metadata['sequence_count'])), n_processes)
            #frames_split = split_frames(metadata['sequence_count'], n_processes)
            processes = [None]*n_processes
            
            for i in range(0, n_processes):
                shared[i] = mp.Array('d', int(9*len(roi_locations)*len(frames_split[i])))
                processes[i] = (mp.Process(target=main, args=(name, fitter, frames_split[i], roi_locations, shared[i])))

            for p in processes:
                p.start()

            for p in processes:
                p.join()

            counter = 0
            for i, share in enumerate(shared):
                arr = np.frombuffer(share.get_obj()) # mp_arr and arr share the same memor
                results = arr.reshape((len(roi_locations)*len(frames_split[i]), 9))
                #print(results)
                result[counter:counter+len(roi_locations)*len(frames_split[i])] = results
                counter+=len(roi_locations)*len(frames_split[i])

            #print(result)

            result = result[result[:, 3] != 0]
            #print(result)

            print('Time taken: ' + str(round(time.time() - start, 3)) + ' s. Fits done: ' + str(result.shape[0]))

            metadata_filtered = {k: v for k, v in metadata.items() if v is not None}
            del metadata_filtered['time_start']
            del metadata_filtered['time_start_utc']

            ## ROI_locations dict
            ROI_locations_dict = dict(zip(['x', 'y'], roi_locations.T))

            ## Localization dict
            results_dict = {'Localizations': result}

            # %% save everything
            tools.save_to_csv_mat('metadata', metadata_filtered, directory)
            tools.save_to_csv_mat('ROI_locations', ROI_locations_dict, directory)
            tools.save_to_csv_mat('Localizations', results_dict, directory)

    #pr.disable()
    #pr.print_stats(sort='time')
