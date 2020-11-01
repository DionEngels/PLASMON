import numpy as np
from numpy import swapaxes, flip, moveaxis
from scipy.io import savemat


def raw_to_matlab(raw):

    if raw.ndim < 3:
        # raw = swapaxes(raw, 0, 1)
        # raw = flip(raw, axis=0)
        pass
    else:
        raw = moveaxis(raw, 0, -1)
        #  raw = swapaxes(raw, 0, 1)
        # raw = swapaxes(raw, 1, 2)
        # raw = swapaxes(raw, 0, 1)
        # raw = flip(raw, axis=0)
    return raw


# %% Not functions but actual gaussians
def make_gaussian(size, fwhm=3, center=None):
    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    x0 = center[0]
    y0 = center[1]

    return np.exp(-4 * np.log(2) * ((x - x0) ** 2 + (y - y0) ** 2) / fwhm ** 2)


num_stack = 5
for i in range(num_stack):
    if i == 0:
        result = make_gaussian(9, center=[np.random.rand()*7 + 1, np.random.rand()*7 + 1])
    else:
        roi = make_gaussian(9, center=[np.random.rand()*7 + 1, np.random.rand()*7 + 1])
        result = np.dstack((result, roi))

result = moveaxis(result, -1, 0)

result2 = raw_to_matlab(result)

result2_dict = {'test': result2}

savemat('C:/Users/s150127/OneDrive - TU Eindhoven/Werk/MBx Python analysis/Code/test codes/raw_output_test.mat',
        result2_dict, long_field_names=True)
