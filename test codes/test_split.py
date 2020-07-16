# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 13:54:01 2020

@author: s150127
"""

import numpy as np
import sys
from pims import ND2_Reader # reader of ND2 files
sys.path.append('../')
name = "C:/Users/s150127/Downloads/_MBx dataset/1nMimager_newGNRs_100mW.nd2"


ROI_locations = np.load('ROI_locations.npy')
ROI_locations = ROI_locations - 1
#ROI_locations = ROI_locations[0:8, :]

test_ROI_split = np.array_split(ROI_locations, 4)

ND2 = ND2_Reader(name)

metadata = ND2.metadata

frames_list = list(range(metadata['sequence_count']))
frames_list_split = np.array_split(frames_list, 4)

frames = ND2

for i, splitter in enumerate(frames_list_split):
    
    frames_split = frames[splitter]