# -*- coding: utf-8 -*-
"""
Created on Mon August 10 2020

@author: Dion Engels
MBx Python Data Analysis

nd2 reading

Everything related to reading in the .nd2 files

----------------------------

v1.0: split from tools and new version with only nd2reader and no longer nd2_reader (that is a different package)
v1.1: None prevention: 10/08/2020

"""

from pims_nd2 import ND2_Reader
from nd2reader import ND2Reader
from nd2reader.parser import Parser
__self_made__ = True


class ND2ReaderForMetadata(ND2Reader):
    """
    The ND2Reader based reader purely for metadata
    """
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
        """
        Get metadata. Reads out nd2 and returns the metadata
        """
        metadata_dict = self.metadata

        metadata_dict.pop('rois', None)
        try:
            metadata_dict['z_levels'] = list(metadata_dict.pop('z_levels'))
            metadata_dict['z_coordinates'] = metadata_dict.pop('z_coordinates')
        except Exception:
            pass
        metadata_dict.pop('frames', None)
        metadata_dict.pop('date', None)

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

        for key, value in metadata_dict.items(): # prevent None values by making None string
            if value is None:
                metadata_dict[key] = str(value)

        metadata_dict['Others'] = metadata_dict_sequence

        return metadata_dict

    def recursive_add_to_dict(self, dictionary):
        """
        Function to parse dictionaries in metadata since many include dictionaries in dictionaries

        :param dictionary: the dictionary to be parsed
        :return: the result of the dictionary
        """
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
        """
        Parses the metadata info of the sequence

        :param info_to_parse: the info to parse from nd2 file
        :return: the parsed metadata
        """
        main_part = info_to_parse[b'SLxPictureMetadata']

        metadata_text_dict = self.recursive_add_to_dict(main_part)

        return metadata_text_dict

    @staticmethod
    def parse_text_info(info_to_parse):
        """
        Parses the metadata info of the image

        :param info_to_parse: the info to parse from nd2 file
        :return: the parsed metadata
        """
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


class ND2ReaderSelf(ND2_Reader):
    """
    Self-made ND2 reader with improved metadata and less warnings
    """
    def __init__(self, filename, series=0, channel=0):
        self._clear_axes()
        self._get_frame_dict = dict()
        super().__init__(filename, series=series, channel=channel)

    def get_metadata(self):
        metadata_nd2 = ND2ReaderForMetadata(self.filename)
        metadata_dict = metadata_nd2.get_metadata()
        metadata_nd2.close()

        return metadata_dict