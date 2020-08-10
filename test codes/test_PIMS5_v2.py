from pims import ND2Reader_SDK
from nd2reader import ND2Reader
from nd2reader.parser import Parser
from pims_nd2 import ND2_Reader

from csv import writer  # to save to csv
from scipy.io import savemat  # to export for MATLAB
from os import mkdir
import datetime

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


class ND2ReaderSelfV3(ND2Reader):

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
        try:
            metadata_dict['z_levels'] = list(metadata_dict.pop('z_levels'))
            metadata_dict['z_coordinates'] = metadata_dict.pop('z_coordinates')
        except Exception:
            pass
        metadata_dict.pop('frames')
        metadata_dict.pop('date')

        metadata_dict['pfs_status'] = self._parser._raw_metadata.pfs_status
        metadata_dict['pfs_offset'] = self._parser._raw_metadata.pfs_offset

        metadata_dict['timesteps'] = self.timesteps
        metadata_dict['frame_rate'] = self.frame_rate

        info_to_parse = self.parser._raw_metadata.image_text_info
        metadata_text_dict = self.parse_text_info(info_to_parse)

        metadata_dict = {**metadata_dict, **metadata_text_dict}

        info_to_parse = self.parser._raw_metadata.image_metadata_sequence
        metadata_dict_sequence = self.parse_sequence_info(info_to_parse)
        try:
            metadata_dict['EnableGainMultiplier'] = metadata_dict_sequence.pop('EnableGainMultiplier')
            metadata_dict['GainMultiplier'] = metadata_dict.pop('Multiplier')
            metadata_dict['Conversion_Gain'] = metadata_dict.pop('Conversion_Gain')
        except Exception:
            pass
        metadata_dict['Others'] = metadata_dict_sequence

        return metadata_dict

    def recursive_add_to_dict(self, dictionary):
        metadata_text_dict = {}
        for key_decoded, value_decoded in dictionary.items():
            if type(key_decoded) is bytes:
                key_decoded = key_decoded.decode("utf-8")
            if type(value_decoded) is bytes:
                value_decoded = value_decoded.decode("utf-8")

            if type(value_decoded) == dict:
                return_dict = self.recursive_add_to_dict(value_decoded)
                metadata_text_dict = {**metadata_text_dict, **return_dict}
            elif type(value_decoded) != str:
                metadata_text_dict[key_decoded] = value_decoded
            else:
                pass

        return metadata_text_dict

    def parse_sequence_info(self, info_to_parse):
        main_part = info_to_parse[b'SLxPictureMetadata']

        metadata_text_dict = self.recursive_add_to_dict(main_part)

        return metadata_text_dict

    @staticmethod
    def parse_text_info(info_to_parse):
        main_part = info_to_parse[b'SLxImageTextInfo']
        metadata_text_dict = {}

        for key, value in main_part.items():
            value_string = value.decode("utf-8")
            if value_string != '':
                split_string = value_string.split('\r\n')
                for line_number, line in enumerate(split_string):
                    if line == '':
                        continue
                    if line_number == 0:  # these are the headers, they do not have a value, only key
                        key = key.decode("utf-8")
                        if '5' in key or '6' in key:
                            continue
                        elif '9' in key:
                            value = value_string
                            key = 'date'
                        elif '13' in key:
                            value = value_string
                            key = 'Objective'
                        metadata_text_dict[key] = value
                        continue
                    try:
                        key, value = line.split(':')  # try to split, works for most
                    except Exception as e:
                        if "too many" in str(e):
                            split_line = line.split(':')  # microscope name has a : in it
                            key = split_line[0]
                            value = ":".join(split_line[1:])
                        elif "not enough" in str(e):
                            continue
                    if key == 'Metadata:' or key == '':  # remove emtpy stuff and the metadata header key
                        continue
                    key = key.lstrip()  # remove the spaces at the start of some keys
                    if type(value) is str:
                        value = value.lstrip()  # remove the spaces at the start of some value that are strings
                    key = key.replace(", ", "_")
                    key = key.replace(" ", "_")  # prevent spaces in key name, matlab does not like that
                    if ',' in value:
                        value = float(value.replace(',', '.'))
                        metadata_text_dict[key] = value
                    else:
                        try:
                            metadata_text_dict[key] = int(value)
                        except Exception:
                            metadata_text_dict[key] = value

        return metadata_text_dict


def save_to_csv_mat_metadata(name, values, path):
    """
    Saver to .csv and .mat for metadata

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
        writer_item = writer(csv_file)
        writer_item.writerows(values.items())

        values_dict = {name: values}

        savemat(path + "/" + name + '.mat', values_dict, long_field_names=True)


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

        nd2_old = ND2Reader_SDK(name)
        nd2_self_v2 = ND2ReaderSelfV2(name)
        nd2_self_v3 = ND2ReaderSelfV3(name)

        metadata_old = nd2_old.metadata
        metadata_self_v2 = nd2_self_v2.get_metadata()

        metadata_old_filtered = {k: v for k, v in metadata_old.items() if v is not None}
        del metadata_old_filtered['time_start']
        del metadata_old_filtered['time_start_utc']

        metadata_self_v3 = nd2_self_v3.get_metadata()

        save_to_csv_mat_metadata('metadata_old', metadata_old_filtered, path)
        save_to_csv_mat_metadata('metadata_self_v2', metadata_self_v2, path)
        save_to_csv_mat_metadata('metadata_self_v3', metadata_self_v3, path)
