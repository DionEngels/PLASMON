def test(x):
    y = 0
    for i in range(x):
        y += i
    return y

import numpy as np
from math import exp
def gaussian(h, c_x, c_y, w_x, w_y, s, d):
    out = np.zeros(s*s)

    for i in range(s):
        for j in range(s):
            out[i*s+j] = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2) -d[i, j]
    return out