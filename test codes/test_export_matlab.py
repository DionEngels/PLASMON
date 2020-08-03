from scipy.io import savemat  # to export for MATLAB
from os import getcwd
import numpy as np

def save_to_csv_mat_roi(name, rois):
    """
    Saves the ROIs to a .mat and .csv

    Parameters
    ----------
    name : name to save to
    rois : rois to save
    height : height of video to switch y-axis
    path : path to save

    Returns
    -------
    None.

    """

    path = getcwd()

    rois_dict = {"rois": rois}
    savemat(path + "/" + name + '.mat', rois_dict, do_compression=True)

test = np.zeros((7,2))


save_to_csv_mat_roi("_test", test)