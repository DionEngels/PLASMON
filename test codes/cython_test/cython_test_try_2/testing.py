import example_cy
import numpy as np
import mbx_fortran
import optimized_cy
import time

size = 9

height = 1
center_x = 4
center_y = 4
width_x = 2
width_y = 2
data = np.zeros((size, size), dtype=np.uint16)

num_loop = 10000

loops = list(range(0, num_loop))

start = time.time()

for loop in loops:
    res_c = example_cy.gaussian(height, center_x, center_y, width_x, width_y, size, data)

res_c2 = np.reshape(res_c, (size, size))
print('Time taken base c: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    res_f = mbx_fortran.gaussian(height, center_x, center_y, width_x, width_y, size, data)

res_f2 = np.reshape(res_f, (size, size))
print('Time taken base c: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

start = time.time()

for loop in loops:
    res_c_opt = optimized_cy.gaussian(height, center_x, center_y, width_x, width_y, size, data)


res_c_opt2 = np.reshape(res_c_opt, (size, size))
print('Time taken base c: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))
