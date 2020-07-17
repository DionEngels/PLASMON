# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 00:38:19 2020

@author: s150127
"""
# import sys
import numpy as np
from pims import ND2_Reader # reader of ND2 files
import time
import multiprocessing as mp
import os

import _code.fitters as fitting
import _code.roi_finding as roi_finding
import _code.tools as tools

# pr = cProfile.Profile()

# %% Init
ROI_SIZE = 7
WAVELENGTH = 637 #nm
filenames = ("C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2",)

n_processes = int(mp.cpu_count()/4)


#%% Main

def mt_main(name, fitter, frames_split, roi_locations, shared, q):
    nd2 = ND2_Reader(name)
    frames = nd2
    metadata = nd2.metadata
    metadata['sequence_count'] = len(frames_split)
    frames = frames[frames_split]

    local_result = fitter.main(frames, metadata, roi_locations, q=q, start_frame=frames_split[0])

    #q.put([frames_split[0], 5])

    for result_index, result in enumerate(local_result):
        shared[9 * result_index:9 * (result_index + 1)] = result[:]

    print("Done fiting")
    
    # q.put([frames_split[0], 25])


# %% Main
if __name__ == '__main__':
    
    # pr.enable()

    for name in filenames:
        with ND2_Reader(name) as ND2:

            basedir = os.getcwd()
            directory = name.split(".")[0].split("/")[-1]
            path = os.path.join(basedir, directory)
            try:
                os.mkdir(path)
            except:
                pass

            # parse ND2 info
            frames = ND2
            metadata = ND2.metadata
            # frames = frames[0:2]
            
            roi_finder = roi_finding.roi_finder(ROI_SIZE, frames[0])#, intensity_min = 800)
            fitter = fitting.scipy_last_fit_guess(metadata, ROI_SIZE,
                                                  WAVELENGTH, roi_finder.intensity_min, 
                                                  "ScipyLastFitGuess", 5)
            roi_locations = roi_finder.main(frames[0], fitter)
            
            #roi_locations = roi_locations[0:4, :]
            
            #fitter = fitting.phasor_only_ROI_loop(metadata, ROI_SIZE, WAVELENGTH, roi_finder.intensity_min, "Phasor")
            
            start = time.time()
            q = mp.Queue()

            shared = [None] * n_processes
            results = np.zeros((len(roi_locations) * len(frames), 9))

            end_frame = 100 #metadata['sequence_count']

            frames_split = np.array_split(list(range(end_frame)), n_processes)
            processes = [None] * n_processes

            for i in range(0, n_processes):
                shared[i] = mp.Array('d', int(9 * len(roi_locations) * len(frames_split[i])))
                processes[i] = (mp.Process(target=mt_main, args=(name, fitter, frames_split[i], roi_locations, shared[i], q)))

            for p in processes:
                p.start()

            queue_dict = {}
            while sum(queue_dict.values()) < 100:
                for p in processes:
                    queue = q.get()
                    queue_dict[queue[0]] = queue[1]
                    print(sum(queue_dict.values()))

            for p in processes:
                p.join()

            counter = 0
            for i, share in enumerate(shared):
                arr = np.frombuffer(share.get_obj())  # mp_arr and arr share the same memory
                result = arr.reshape((len(roi_locations) * len(frames_split[i]), 9))
                results[counter:counter + len(roi_locations) * len(frames_split[i])] = result
                counter += len(roi_locations) * len(frames_split[i])

            results = results[results[:, 3] != 0]
            
            #print(result)
            
            print('Time taken: ' + str(round(time.time() - start, 3)) + ' s. Fits done: ' + str(results.shape[0]))
                
            #pr.disable()
            #pr.print_stats(sort='time')

            metadata_filtered = {k: v for k, v in metadata.items() if v is not None}
            del metadata_filtered['time_start']
            del metadata_filtered['time_start_utc']

            ## ROI_locations dict
            ROI_locations_dict = dict(zip(['x', 'y'], roi_locations.T))

            ## Localization dict
            results_dict = {'Localizations': results}

            # %% save everything
            tools.save_to_csv_mat('metadata', metadata_filtered, directory)
            tools.save_to_csv_mat('ROI_locations', ROI_locations_dict, directory)
            tools.save_to_csv_mat('Localizations', results_dict, directory)
