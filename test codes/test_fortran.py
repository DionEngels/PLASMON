# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 21:15:33 2020

@author: s150127
"""

import numpy as np

import fib2

print(fib2.fib.__doc__)

a = np.zeros(11)
print(fib2.fib(8))

import fib3

test = fib3.fib(8)

import gauss4 as gauss

data = np.zeros((9,9))

data = gauss.gaussian(1, 4, 4, 1, 1)

print(data)