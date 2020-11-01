import numpy as np
import optimized_cy
import time
from scipy.ndimage import median_filter

shape = (1024, 1024)
to_do = np.random.random(shape) * 10000
to_do = np.asarray(to_do, dtype=np.uint16)

#for i in range(shape[0]):
#    for j in range(shape[1]):
#        to_do[i, j] = i+j

num_loop = 1
size = 9

loops = list(range(0, num_loop))

start = time.time()

for loop in loops:
    res = median_filter(to_do, size=size)

print('Time taken base: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

from statistics import median


def median_filter_self(d, size, size_y, size_x):
    out = np.zeros_like(d, dtype=np.uint16)

    size_1d = int((size - 1) / 2)
    array = [0]*size*size

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
            out[i, j] = median(array)

    return out

#start = time.time()

#for loop in loops:
#    res_self = median_filter_self(to_do, size, shape[0], shape[1])

#print('Time taken self: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    res_c = optimized_cy.median_filter(to_do, size, shape[0], shape[1])

print('Time taken c: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))