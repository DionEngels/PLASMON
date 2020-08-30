from pims import ND2Reader_SDK
import pywt
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
import numpy as np

filenames = ("C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2",)

if __name__ == "__main__":
    for name in filenames:
        path = name.split(".")[0]

        nd2 = ND2Reader_SDK(name)

        frame_zero = np.asarray(nd2[0], dtype=np.int32)
        background = np.zeros(frame_zero.shape, dtype=np.int64)

        Wavelets = ['bior3.3', 'bior3.5', 'bior3.7', 'bior3.9', 'bior4.4', 'bior5.5', 'bior6.8']

        for wavelet in Wavelets:
            background_pywt = pywt.wavedec2(frame_zero, wavelet)
            if wavelet == 'bior3.3':
                background = pywt.waverec2(background_pywt, wavelet)
                background = background
            background += pywt.waverec2(background_pywt, wavelet)

        background /= len(Wavelets)
        background2 = median_filter(frame_zero, size=9, mode='constant')

        fig_base, ax_base = plt.subplots(1)
        ax_base.imshow(frame_zero, extent=[0, frame_zero.shape[1], frame_zero.shape[0], 0], aspect='auto')
        plt.show()

        result = frame_zero - background
        fig, ax = plt.subplots(1)
        ax.imshow(result, extent=[0, background.shape[1], background.shape[0], 0], aspect='auto')
        plt.show()

        result2 = frame_zero - background2
        fig2, ax2 = plt.subplots(1)
        ax2.imshow(result2, extent=[0, background2.shape[1], background2.shape[0], 0], aspect='auto')
        plt.show()
