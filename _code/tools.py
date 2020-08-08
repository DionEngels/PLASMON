# -*- coding: utf-8 -*-
"""
Created on Sun May 31 23:27:01 2020

@author: Dion Engels
MBx Python Data Analysis

Tools

Some additional tools used by MBx Python

----------------------------

v0.1: Save to CSV & Mat: 31/05/2020
v0.2: also switch array: 04/06/2020
v0.3: cleaned up: 24/07/2020
v0.4: settings and results text output: 25/07/2020
v0.5: own ND2Reader class to prevent warnings: 29/07/2020
v0.6: save drift and save figures: 03/08/2020
v0.6.1: better MATLAB ROI and Drift output
v0.7: metadata v2
v0.8: metadata v3
v1.0: bug fixes
v1.0.1: metadata goldmine found, to be implemented

"""
from numpy import zeros
from datetime import datetime

from pims_nd2 import ND2_Reader
from nd2reader import ND2Reader
from nd2reader.parser import Parser


# %% Class to load in ND2s


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

        metadata_dict.pop('rois', None)
        metadata_dict.pop('z_levels', None)
        metadata_dict.pop('frames', None)
        metadata_dict.pop('date', None)

        metadata_dict['pfs_status'] = self.parser._raw_metadata.pfs_status
        metadata_dict['pfs_offset'] = self.parser._raw_metadata.pfs_offset

        info_to_parse = self.parser._raw_metadata.image_text_info
        metadata_text_dict = self.parse_text_info(info_to_parse)

        info_to_parse2 = self.parser._raw_metadata.image_metadata_sequence
        metadata_text_dict2 = self.parse_text_info2(info_to_parse2)

        metadata_dict['timesteps'] = self.timesteps.tolist()
        metadata_dict['frame_rate'] = self.frame_rate

        return metadata_dict, metadata_text_dict

    def recursive_add_to_dict(self, dictionary):
        metadata_text_dict = {}
        for key_decoded, value_decoded in dictionary.items():
            if type(key_decoded) is bytes:
                key_decoded = key_decoded.decode("utf-8")
            if type(value_decoded) is bytes:
                value_decoded = value_decoded.decode("utf-8")

            if type(value_decoded) == dict:
                return_dict = self.recursive_add_to_dict(value_decoded)
                metadata_text_dict[key_decoded] = return_dict
            elif type(value_decoded) != str:
                metadata_text_dict[key_decoded] = value_decoded
            else:
                pass

        return metadata_text_dict

    def parse_text_info2(self, info_to_parse):
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
                for line in split_string:
                    if line == '':
                        continue
                    try:
                        key, value = line.split(':')  # try to split, works for most
                    except Exception as e:
                        if "too many" in str(e):
                            try:
                                _ = datetime.strptime(line, "%d-%m-%Y  %H:%M:%S")  # date
                                key = 'date'
                                value = line
                            except:
                                split_line = line.split(':')  # microscope name has a : in it
                                key = split_line[0]
                                value = ":".join(split_line[1:])
                        elif "not enough" in str(e):
                            continue
                    if key == 'Metadata' or key == '':  # remove emtpy stuff and the metadata header key
                        continue
                    key = key.lstrip()  # remove the spaces at the start of some keys
                    if type(value) is str:
                        value = value.lstrip()  # remove the spaces at the start of some value that are strings
                    key = key.replace(" ", "_")  # prevent spaces in key name, matlab does not like that
                    metadata_text_dict[key] = value

        return metadata_text_dict


class ND2ReaderSelf(ND2_Reader):
    """
    Small class to read in ND2 using a prebuild ND2 Reader. Slightly edited to prevent it giving a warning
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
        metadata_dict_part2, metadata_text_dict = nd2_part_2.get_metadata()
        total_metadata = {**metadata_dict_filtered, **metadata_text_dict, **metadata_dict_part2}

        nd2_part_2.close()

        return total_metadata


def switch(array):
    """
    Switches a single arrays values

    Parameters
    ----------
    array : array to switch

    Returns
    -------
    new : switched array

    """
    new = zeros(array.shape)
    new[:, 1] = array[:, 0]
    new[:, 0] = array[:, 1]
    return new
