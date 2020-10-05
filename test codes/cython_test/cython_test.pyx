# %%cython -a --cplus
cimport numpy as np
import numpy as np
import cython

cdef extern from "math.h":
    cpdef double exp(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef gaussian_c(double h, double c_x, double c_y, double w_x, double w_y, int s, np.ndarray [np.uint16_t, ndim=2] d):
    cdef int i, j
    cdef np.ndarray [np.float64_t, ndim=2] out = np.empty_like(d, dtype=np.float)
    for i in range(s):
        for j in range(s):
            out[i, j] = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2) -d[i, j]
    return out

def gaussian(double h, double c_x, double c_y, double w_x, double w_y, int s, np.ndarray [np.uint16_t, ndim=2] d):
    return gaussian_c(h, c_x, c_y, w_x, w_y, s, d)