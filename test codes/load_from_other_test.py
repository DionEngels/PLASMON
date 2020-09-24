import tkinter as tk
import numpy as np
from tkinter.filedialog import askopenfilenames
from os import getcwd
import sys
from scipy.io import loadmat
from scipy.ndimage.interpolation import shift
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.fft as fft
from scipy.signal import fftconvolve
from scipy.ndimage import median_filter

FILETYPES_LOAD_FROM_OTHER = [(".npy and .mat", ".npy"), (".npy and .mat", ".mat")]


def normxcorr2(b, a):
    """
    Correlation of similar size frames
    """
    def conv2(a, b):
        ma, na = a.shape
        mb, nb = b.shape
        return fft.ifft2(fft.fft2(a, [2 * ma - 1, 2 * na - 1]) * fft.fft2(b, [2 * mb - 1, 2 * nb - 1]))

    c = conv2(a, np.flipud(np.fliplr(b)))
    a = conv2(a ** 2, np.ones(b.shape))
    b = sum(b.flatten() ** 2)
    c = c / np.sqrt(a * b)
    return c


def normxcorr2_large(template, image, mode="full"):
    # If this happens, it is probably a mistake
    if np.ndim(template) > np.ndim(image) or \
            len([i for i in range(np.ndim(template)) if template.shape[i] > image.shape[i]]) > 0:
        print("normxcorr2: TEMPLATE larger than IMG. Arguments may be swapped.")

    template = template - np.mean(template)
    image = image - np.mean(image)

    a1 = np.ones(template.shape)
    # Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
    ar = np.flipud(np.fliplr(template))
    out = fftconvolve(image, ar.conj(), mode=mode)

    image = fftconvolve(np.square(image), a1, mode=mode) - \
            np.square(fftconvolve(image, a1, mode=mode)) / (np.prod(template.shape))

    # Remove small machine precision errors after subtraction
    image[np.where(image < 0)] = 0

    template = np.sum(np.square(template))
    out = out / np.sqrt(image * template)

    # Remove any divisions by 0 or very close to 0
    out[np.where(np.logical_not(np.isfinite(out)))] = 0

    return out


def plot_rois(frame, roi_locations, roi_size):
    """

    Parameters
    ----------
    frame : frame to plot
    roi_locations : locations to draw box around
    roi_size : Size of boxes to draw

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(1)
    ax.imshow(frame, extent=[0, frame.shape[1], frame.shape[0], 0], aspect='auto')
    roi_size_1d = int((roi_size - 1) / 2)

    roi_locations = roi_locations - roi_size_1d

    for roi in roi_locations:
        rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                 linewidth=0.5, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    plt.title("ROI locations")
    plt.show()


def switch_axis(array):
    new = np.zeros(array.shape)
    new[:, 1] = array[:, 0]
    new[:, 0] = array[:, 1]
    return new


def roi_to_python_coordinates(roi_locs, height):
    roi_locs[:, 1] = height - roi_locs[:, 1]
    roi_locs = switch_axis(roi_locs)

    return roi_locs

"""
tk.Tk().withdraw()
while True:
    filenames = askopenfilenames(filetypes=FILETYPES_LOAD_FROM_OTHER,
                                 title="Select both frame_zero and ROI locations that you want to use",
                                 initialdir=getcwd())
    if len(filenames) == 2:
        break
    else:
        check = tk.messagebox.askokcancel("ERROR",
                                          "Select one ROI locations file and one frame_zero. Try again?")
        if not check:
            sys.exit(0)
"""


def load_in_files(filenames):
    frame_zero_old = np.load([file for file in filenames if file.endswith('.npy')][0])
    roi_locations = loadmat([file for file in filenames if file.endswith('.mat')][0])
    roi_locations = [roi_locations[key] for key in roi_locations.keys() if not key.startswith('_')][0]
    roi_locations = roi_to_python_coordinates(roi_locations, frame_zero_old.shape[0])
    return roi_locations, frame_zero_old


def correlate_frames(frame_old, frame_new):
    if frame_old.shape == frame_new.shape:
        corr = normxcorr2(frame_old, frame_new)
        maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
        offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
    else:
        corr = normxcorr2_large(frame_old, frame_new)
        maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
        offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
    return offset


filenames = ('C:/Users/s150127/Downloads/___MBx/datasets/AuNR_11bpAtto647N_1nM_637nm_100mW_dock/frame_zero.npy',
             'C:/Users/s150127/Downloads/___MBx/datasets/AuNR_11bpAtto647N_1nM_637nm_100mW_dock/ROI_locations.mat')

filename_new = ('C:/Users/s150127/Downloads/___MBx/datasets/AuNR_11bpAtto647N_1nM_637nm_10mW_dock/frame_zero.npy',)

roi_locations, frame_zero_old = load_in_files(filenames)

frame_zero = np.load(filename_new[0])

plot_rois(frame_zero_old, roi_locations, 9)

background = median_filter(frame_zero_old, size=9)
frame_zero_old = frame_zero_old.astype(np.int16) - background

plot_rois(frame_zero_old, roi_locations, 9)

background = median_filter(frame_zero, size=9)
frame_zero = frame_zero.astype(np.int16) - background

offset = correlate_frames(frame_zero_old, frame_zero)

if offset[0] < 50 and offset[1] < 50:
    try:
        print(offset)
        roi_locations += offset
    except:
        print("ERROR, You did not select a proper ROI locations file. Try again.")
else:
    load_in_files(filenames)
    frame_zero = np.load(filename_new[0])
    background = median_filter(frame_zero_old, size=9, mode='constant', cval=np.mean(frame_zero_old))
    frame_zero_old = frame_zero_old.astype(np.int16) - background
    background = median_filter(frame_zero, size=9, mode='constant', cval=np.mean(frame_zero))
    frame_zero = frame_zero.astype(np.int16) - background
    offset = correlate_frames(frame_zero_old, frame_zero)
    try:
        print(offset)
        roi_locations += offset
    except:
        print("ERROR, You did not select a proper ROI locations file. Try again.")

plot_rois(frame_zero, roi_locations, 9)
