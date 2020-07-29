from pims import ND2Reader_SDK, ND2_Reader
from nd2reader import ND2Reader
from pims_nd2 import ND2_Reader

from csv import DictWriter  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from os import mkdir

filenames = ("C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2",)


class ND2ReaderSelf(ND2_Reader):
    """
    Small class to read in ND2 using a prebuild ND2 Reader. Slightly edited to prevent it giving a warning
    """
    def __init__(self, filename, series=0, channel=0):
        self._clear_axes()
        self._get_frame_dict = dict()
        super().__init__(filename, series=series, channel=channel)


def save_to_csv_mat(name, values, path):
    """
    Basic saver to .csv and .mat, only used by metadata

    Parameters
    ----------
    name : name to save to
    values : values to save
    path : path to save

    Returns
    -------
    None.

    """
    with open(path + "/" + name + '.csv', mode='w') as csv_file:
        fieldnames = [k[0] for k in values.items()]
        writer = DictWriter(csv_file, fieldnames=fieldnames)

        #  writer.writeheader()
        writer.writerow(values)

        savemat(path + "/" + name + '.mat', values)


if __name__ == "__main__":
    for name in filenames:
        path = name.split(".")[0]

        directory_try = 0
        directory_success = False
        while not directory_success:
            try:
                mkdir(path)
                directory_success = True
            except:
                directory_try += 1
                if directory_try == 1:
                    path += "_%03d" % directory_try
                else:
                    path = path[:-4]
                    path += "_%03d" % directory_try

        #  nd2_new = ND2Reader(name)
        #  nd2_old = ND2Reader_SDK(name)
        #  nd2_alt = ND2_Reader(name)
        nd2_self = ND2ReaderSelf(name)

        #  metadata_new = nd2_new.metadata
        #  metadata_old = nd2_old.metadata
        #  metadata_alt = nd2_alt.metadata
        metadata_self = nd2_self.metadata

        #  metadata_old_filtered = {k: v for k, v in metadata_old.items() if v is not None}

        #  del metadata_old_filtered['time_start']
        #  del metadata_old_filtered['time_start_utc']

        #  metadata_new_filtered = {k: v for k, v in metadata_new.items() if v is not None}

       #   metadata_new_filtered.pop('rois')
       #   metadata_new_filtered.pop('z_levels')
    #  metadata_new_filtered.pop('frames')

        #  metadata_alt_filtered = {k: v for k, v in metadata_alt.items() if v is not None}

        #  del metadata_alt_filtered['time_start']
        #  del metadata_alt_filtered['time_start_utc']

        metadata_self_filtered = {k: v for k, v in metadata_self.items() if v is not None}

        del metadata_self_filtered['time_start']
        del metadata_self_filtered['time_start_utc']

      #    save_to_csv_mat('metadata_new', metadata_new_filtered, path)
        #  save_to_csv_mat('metadata_old', metadata_old_filtered, path)
        #  save_to_csv_mat('metadata_alt', metadata_alt_filtered, path)
        save_to_csv_mat('metadata_self', metadata_self_filtered, path)

      #    nd2_new.close()
        #  nd2_old.close()

