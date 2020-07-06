# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 17:09:23 2020

@author: s150127
"""

from functools import lru_cache
import time

num_loop = 100000

loops = list(range(0, num_loop))

#%%

@lru_cache(maxsize=32)
def fib_cache(n):
    if n < 2:
        return n
    return fib_cache(n-1) + fib_cache(n-2)

start = time.time()
    
for loop in loops:
    [fib_cache(n) for n in range(10)]
    
print('Time taken v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

# Output: [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]

#%%
def fib(n):
    if n < 2:
        return n
    return fib(n-1) + fib(n-2)

start = time.time()
    
for loop in loops:
    [fib(n) for n in range(10)]
    
print('Time taken v3: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))


