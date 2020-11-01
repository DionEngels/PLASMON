import numpy as np
cimport numpy as np
import cython
from statistics import median

cdef extern from "stdlib.h":
    cpdef int qsort()

cdef quick_select_median(int arr[], int n):
    cdef int low, high
    cdef int median;
    cdef int middle, ll, hh;
    low = 0
    high = n-1
    median = (low + high) / 2
    middle = (low + high) / 2
    if arr[middle] > arr[high]:
        ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])
    ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])
    ELEM_SWAP(arr[middle], arr[low]) ;
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
    do ll++; while (arr[low] > arr[ll]) ;
    do hh--; while (arr[hh] > arr[low]) ;
    if (hh < ll)
    break;
    ELEM_SWAP(arr[ll], arr[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;
    /* Re-set active partition */
    if (hh <= median)
    low = ll;
    if (hh >= median)
    high = hh - 1;
    }
    return arr[median] ;
}

@cython.boundscheck(False)
@cython.wraparound(False)
def median_filter(np.ndarray [np.uint16_t, ndim=2] d, int size, int size_y, int size_x):
    cdef int i, j, y_loop, x_loop, size_1d, y_select, x_select
    cdef np.ndarray [np.uint16_t, ndim=2] out = np.zeros_like(d, dtype=np.uint16)
    cdef np.ndarray [np.uint16_t, ndim=1] array = np.zeros(size*size, dtype=np.uint16)

    size_1d = (size - 1) / 2

    for i in range(size_y):
        for j in range(size_x):
            for y_loop in range(size):
                for x_loop in range(size):
                    if i + y_loop - size_1d < 0:
                        y_select = -(i + y_loop - size_1d)
                    elif i + y_loop - size_1d >= size_y:
                        y_select = (size_y - 1) - (i + y_loop - size_1d - (size_y - 1))
                    else:
                        y_select = i + y_loop - size_1d
                    if j + x_loop - size_1d < 0:
                        x_select = -(j + x_loop - size_1d)
                    elif j + x_loop - size_1d >= size_x:
                        x_select = (size_x - 1) - (j + x_loop - size_1d - (size_x - 1))
                    else:
                        x_select = j + x_loop - size_1d
                    array[y_loop*size+x_loop] = d[y_select, x_select]
            array = qsort(array, size*size, sizeof(int))
            out[i, j] = median(array)
    return out