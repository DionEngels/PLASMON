# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 10:39:50 2020

@author: s150127
"""

from collections import OrderedDict
import time

class LimitedSizeDict(OrderedDict):
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)
                
class LimitedSizeDict2(OrderedDict):
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        #self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)               
num = 1000000
limit = num

start = time.time()

dictionary = LimitedSizeDict(size_limit = limit)
for i in range(0,num):
    dictionary[tuple([i,i,i,i])]= [i,i,i,i,i,i,i,i,i,i,i,i]
    
print(str(round(time.time() - start, 3)))  
start = time.time()

dictionary_2 = LimitedSizeDict2(size_limit = limit)
for i in range(0,num):
    dictionary_2[tuple([i,i,i,i])]= [i,i,i,i,i,i,i,i,i,i,i,i]
    if i%1000 ==0:
        dictionary_2._check_size_limit()

print(str(round(time.time() - start, 3)))  
start = time.time()

norm_dict = {}
for i in range(0,num):
    norm_dict[tuple([i,i,i,i])]= [i,i,i,i,i,i,i,i,i,i,i,i]
    
print(str(round(time.time() - start, 3)))