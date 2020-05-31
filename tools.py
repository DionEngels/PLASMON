# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
MBx Python Data Analysis

Tools

----------------------------

v1: Save to CSV & Mat: 31/05/2020

"""

import csv # to save to csv
import scipy.io as sio #to export for MATLAB

def SaveToCsvMat(name,values, directory):
    with open(directory + "/"+ name + '.csv', mode='w') as csv_file:
        fieldnames = [k[0] for k in values.items()]
        writer = csv.DictWriter(csv_file, fieldnames = fieldnames)
        
        writer.writeheader()
        writer.writerow(values)
        
        sio.savemat(directory + "/"+ name + '.mat',values)
    
    