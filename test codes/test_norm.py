# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 13:24:49 2020

@author: s150127
"""
import numpy as np
from numpy.linalg import norm
import time

num_loop = 1000000

loops = list(range(0, num_loop))

def make_array():
    return (np.random.rand(5)-0.5)*100

def own_checker(g):
    
    maximum = 0
    
    for number in g:
        if number > maximum:
            maximum = number
        elif -number > maximum:
            maximum = -number
            
    return maximum

start = time.time()

for loop in loops:
    g1 = make_array()
    g_norm1 = norm(g1, ord=np.inf)

print('Time taken np: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    g2 = make_array()
    g_norm2 = np.max(np.abs(g2))

print('Time taken abs/max: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    g3 = make_array()
    g_norm3 = abs(g3).max()

print('Time taken abs/max v2: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    g3 = make_array()
    g_norm3 = np.abs(g3).max()

print('Time taken abs/max v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    g4 = make_array()
    g_norm4 = own_checker(g4)

print('Time taken self: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

import MBx_FORTRAN_TEST_v2 as fortran

start = time.time()

for loop in loops:
    g5 = make_array()
    g_norm5 = fortran.norm5(g5)

print('Time taken fortran: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

import MBx_FORTRAN_TEST_v5 as fortran2
start = time.time()

for loop in loops:
    g6 = make_array()
    g_norm6 = fortran2.norm5_2(g6)

print('Time taken fortran v2: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

import MBx_FORTRAN_TEST_v6 as fortran3
start = time.time()

for loop in loops:
    g7 = make_array()
    g_norm7 = fortran3.norm5_3(g7)

print('Time taken fortran v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))
