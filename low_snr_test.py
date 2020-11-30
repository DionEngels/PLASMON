# -*- coding: utf-8 -*-
"""
Created on 28/11/2020

----------------------------

Author: s150127
Project: PLASMON
File: low_snr_test

Description:
Low SNR test. Creates Gaussians with lowSNR. Tries to fit them. Check how bad their fits were.
----
To run, you have to move this to the main folder
"""
import src.tt as fitting
import src.hsm as hsm
import numpy as np
import matplotlib.pyplot as plt


def create_gaussian(x, roi_size, fitter):
    data = np.zeros((roi_size, roi_size))
    return fitter.fun_gaussian(x, data)


def create_roi(x, roi_size, fitter):
    params = [x[0], 4, 4, 1, 1]
    bg = x[-1]
    gaussian = np.reshape(create_gaussian(params, roi_size, fitter), (roi_size, roi_size))
    bg = np.random.poisson(bg, size=(roi_size, roi_size))
    return gaussian + bg


class TTPartBeun:
    def __init__(self):
        self.frame_start = 0
        self.offset_from_base = [0, 0]


parameter_sets = [[50, 500], [100, 500], [200, 500], [300, 500], [500, 500], [750, 500], [1000, 500],
                  [50, 1000], [100, 1000], [200, 1000], [300, 1000], [500, 1000], [750, 1000], [1000, 1000],
                  [50, 2500], [100, 2500], [200, 2500], [300, 2500], [500, 2500], [750, 2500], [1000, 2500],
                  [50, 5000], [100, 5000], [200, 5000], [300, 5000], [500, 5000], [750, 5000], [1000, 5000]]
loops = 20
roi_size = 7

settings = {'roi_size': roi_size, 'method': 'Gaussian', 'rejection': True}

fitter_create = fitting.Gaussian(settings, 600, 5, [0, 0])
fitter = hsm.HSMFit((7 - 1) // 2)

# tt_part = TTPartBeun()
shape = np.ones(1)
energy_width = np.ones(1)

total_error = np.zeros(len(parameter_sets))
for parameter_index, parameters in enumerate(parameter_sets):
    sigma_error = np.zeros(loops)
    for loop_index, loop in enumerate(range(loops)):
        roi = np.reshape(create_roi(parameters, roi_size, fitter_create), (1, roi_size, roi_size))
        _, _, raw_fits = fitter.fitter(roi, shape, energy_width)
        sigma_error[loop_index] = np.abs((raw_fits[0, 3] + raw_fits[0, 4]) / 2.0 - 1.0)
    total_error[parameter_index] = np.nanmean(sigma_error)
    print(f'Done with parameter {parameter_index} of {len(parameter_sets)}')

fig = plt.figure()
ax1 = fig.add_subplot(111)

step = 7
ax1.scatter([50, 100, 200, 300, 500, 750, 1000], total_error[step * 0:step * 1], s=10, c='r', marker="o", label='Background = 500')
ax1.scatter([50, 100, 200, 300, 500, 750, 1000], total_error[step * 1:step * 2], s=10, c='g', marker="s", label='Background = 1000')
ax1.scatter([50, 100, 200, 300, 500, 750, 1000], total_error[step * 2:step * 3], s=10, c='k', marker="o", label='Background = 2500')
ax1.scatter([50, 100, 200, 300, 500, 750, 1000], total_error[step * 3:step * 4], s=10, c='y', marker="s", label='Background = 5000')
plt.legend(loc='lower left')
plt.title("Settings: Base. Averaged over 20")
plt.xscale('log')
plt.yscale('log')
plt.grid(which='major')
plt.grid(which='minor', linestyle='--', linewidth=0.5)
plt.minorticks_on()
plt.xlabel('Peak intensity Gaussian')
plt.ylabel('Difference between sigma set and measured')
plt.show()
