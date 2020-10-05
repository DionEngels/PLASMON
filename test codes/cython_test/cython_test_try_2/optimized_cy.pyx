import numpy as np
cimport numpy as np
import cython

cdef extern from "math.h":
    cpdef double exp(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def gaussian(double h, double c_x, double c_y, double w_x, double w_y, int s, np.ndarray [np.uint16_t, ndim=2] d):
    cdef int i, j
    cdef np.ndarray [double, ndim=1] out = np.empty(s*s, dtype=np.float)

    for i in range(s):
        for j in range(s):
            out[i*s+j] = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2) -d[i, j]
    return out