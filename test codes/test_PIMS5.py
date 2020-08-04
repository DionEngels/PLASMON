from pims import ND2Reader_SDK, ND2_Reader
from nd2reader import ND2Reader
from nd2reader.parser import Parser
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


class ND2ReaderForMetadata(ND2Reader):

    def __init__(self, filename):
        super(ND2Reader, self).__init__()
        self.filename = filename

        # first use the parser to parse the file
        self._fh = open(filename, "rb")
        self._parser = Parser(self._fh)

        # Setup metadata
        self.metadata = self._parser.metadata

        # Set data type
        self._dtype = self._parser.get_dtype_from_metadata()

        # Setup the axes
        self._setup_axes()

        # Other properties
        self._timesteps = None

    def get_metadata(self):
        metadata_dict = self.metadata

        metadata_dict.pop('rois')
        metadata_dict.pop('z_levels')
        metadata_dict.pop('frames')
        metadata_dict.pop('date')

        metadata_dict['pfs_status'] = self._parser._raw_metadata.pfs_status
        metadata_dict['pfs_offset'] = self._parser._raw_metadata.pfs_offset

        metadata_dict['timesteps'] = self.timesteps
        metadata_dict['frame_rate'] = self.frame_rate

        return metadata_dict


class ND2ReaderSelfV2(ND2_Reader):
    """
    Class to read in ND2 using a prebuild ND2 Reader.
    """
    def __init__(self, filename, series=0, channel=0):
        self._clear_axes()
        self._get_frame_dict = dict()
        super().__init__(filename, series=series, channel=channel)

    def get_metadata(self):
        metadata_dict = self.metadata
        metadata_dict_filtered = {k: v for k, v in metadata_dict.items() if v is not None}
        del metadata_dict_filtered['time_start']
        del metadata_dict_filtered['time_start_utc']

        nd2_part_2 = ND2ReaderForMetadata(self.filename)
        metadata_dict_part2 = nd2_part_2.get_metadata()
        total_metadata = {**metadata_dict_filtered, **metadata_dict_part2}

        nd2_part_2.close()

        return total_metadata


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

        values_dict = {'metadata': values}

        savemat(path + "/" + name + '.mat', values_dict)


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

        nd2_new = ND2Reader(name)
        nd2_old = ND2Reader_SDK(name)
        nd2_alt = ND2_Reader(name)
        nd2_self = ND2ReaderSelf(name)
        nd2_self_v2 = ND2ReaderSelfV2(name)

        metadata_new = nd2_new.metadata
        metadata_old = nd2_old.metadata
        metadata_alt = nd2_alt.metadata
        metadata_self = nd2_self.metadata
        metadata_self_v2 = nd2_self_v2.get_metadata()

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
        save_to_csv_mat('metadata_self_v2', metadata_self_v2, path)

      #    nd2_new.close()
        #  nd2_old.close()

