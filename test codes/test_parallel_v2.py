# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 00:38:19 2020

@author: s150127
"""

import numpy as np
from pims import ND2_Reader # reader of ND2 files
from pims import Frame # for converting ND2 generator to one numpy array
import time
import multiprocessing as mp
#%% Main

def main(i, q):
    q.put(i*i)

if __name__ == '__main__':
    q = mp.Queue()
    processes = [mp.Process(target=main, args=(i, q)) for i in range(2,5)]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    result = [q.get() for p in processes]
    print(result)