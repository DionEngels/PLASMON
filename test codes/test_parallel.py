# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 00:38:19 2020

@author: s150127
"""

import numpy as np
from pims import ND2_Reader # reader of ND2 files
from pims import Frame # for converting ND2 generator to one numpy array
import concurrent.futures
import time

def main(frames, metadata, roi):

        frame_stack_total = Frame(frames)
        ROI_size_1D = 4

        y = int(roi[0])
        x = int(roi[1])
        frame_stack = frame_stack_total[:,y-ROI_size_1D:y+ROI_size_1D+1, x-ROI_size_1D:x+ROI_size_1D+1]

        return frame_stack


def do_something(seconds):
    print(f'Sleeping {seconds} second(s)...')
    time.sleep(seconds)
    return f'Done Sleeping...{seconds}'

if __name__ == '__main__':
    ROI_locations = np.load('ROI_locations.npy')
    ROI_locations = ROI_locations[0:2,:]
    name = "C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2"

    with ND2_Reader(name) as ND2:
        frames = ND2
        metadata = ND2.metadata
        frames = frames[0:10]

        with concurrent.futures.ProcessPoolExecutor() as executor:
            #secs = [5, 4, 3, 2, 1]
            #result = [executor.submit(do_something, sec) for sec in secs]
            result = [executor.submit(main, frames, metadata, roi) for roi in ROI_locations]

            for f in concurrent.futures.as_completed(result):
                print(f.result())
