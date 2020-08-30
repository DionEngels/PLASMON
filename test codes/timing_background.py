from pims import ND2Reader_SDK
from scipy.signal import medfilt
from scipy.ndimage import median_filter
import time

filenames = ("C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2",)

num_loop = 1

loops = list(range(0, num_loop))

if __name__ == "__main__":
    for name in filenames:
        path = name.split(".")[0]

        nd2 = ND2Reader_SDK(name)

        frame_zero = nd2[0]

        start = time.time()

        for loop in loops:
            background = medfilt(frame_zero, kernel_size=9)

        print('Time taken medfilt: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

        start = time.time()

        for loop in loops:
            background2 = median_filter(frame_zero, size=9, mode='constant')

        print('Time taken median_filter: ' + str(round(time.time() - start, 3)) + ' s. Loops: ' + str(len(loops)))

