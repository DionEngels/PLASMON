from pims_nd2 import ND2_Reader
from nd2reader import ND2Reader
from nd2reader.parser import Parser
import numpy as np


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

    def get_metadata(self, verbose=True):
        """
        Get metadata. Reads out nd2 and returns the metadata
        """
        # set warnings to show
        simplefilter('always', DataWarning)
        # get base metadata
        metadata_dict = self.metadata

        # remove ROIs
        metadata_dict.pop('rois', None)

        # try to find z levels
        try:
            metadata_dict['z_levels'] = list(metadata_dict.pop('z_levels'))
            metadata_dict['z_coordinates'] = metadata_dict.pop('z_coordinates')
        except Exception:
            if verbose:
                warn("Z-levels missing from metadata", DataWarning)

        # remove frames and date (for now)
        metadata_dict.pop('frames', None)
        metadata_dict.pop('date', None)

        # check pfs status
        try:
            metadata_dict['pfs_status'] = self._parser._raw_metadata.pfs_status
            metadata_dict['pfs_offset'] = self._parser._raw_metadata.pfs_offset
        except Exception:
            if verbose:
                warn("PFS data missing from metadata", DataWarning)

        # check timesteps and frame rate
        try:
            metadata_dict['timesteps'] = self.timesteps
            metadata_dict['frame_rate'] = self.frame_rate
        except Exception:
            if verbose:
                warn("Timestep data missing from metadata", DataWarning)

        # add more info
        try:
            info_to_parse = self.parser._raw_metadata.image_text_info
            metadata_text_dict = self.parse_text_info(info_to_parse)
            metadata_dict = {**metadata_dict, **metadata_text_dict}
        except Exception:
            if verbose:
                warn("Detailed metadata missing", DataWarning)

        # add raw metadata
        try:
            info_to_parse = self.parser._raw_metadata.image_metadata_sequence
            metadata_dict_sequence = self.parse_sequence_info(info_to_parse)
            try:
                # move some important stuff to main metadata
                metadata_dict['EnableGainMultiplier'] = metadata_dict_sequence.pop('EnableGainMultiplier')
                metadata_dict['GainMultiplier'] = metadata_dict.pop('Multiplier')
                metadata_dict['Conversion_Gain'] = metadata_dict.pop('Conversion_Gain')
            except Exception:
                pass
            # move rest to others
            metadata_dict['Others'] = metadata_dict_sequence
        except Exception:
            if verbose:
                warn("Raw metadata missing", DataWarning)

        # prevent None values by making None string
        for key, value in metadata_dict.items():
            if value is None:
                metadata_dict[key] = str(value)

        return metadata_dict

    def recursive_add_to_dict(self, dictionary):
        """
        Function to parse dictionaries in metadata since many include dictionaries in dictionaries
        ----------------------------------------------
        :param dictionary: the dictionary to be parsed
        :return: the result of the dictionary
        """
        metadata_text_dict = {}
        for key_decoded, value_decoded in dictionary.items():
            # decode keys and values
            if type(key_decoded) is bytes:
                key_decoded = key_decoded.decode("utf-8")
            if type(value_decoded) is bytes:
                value_decoded = value_decoded.decode("utf-8")

            # if dict, add to restart function. Otherwise to add
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
        --------------------------------------------
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
        --------------------------------------------
        :param info_to_parse: the info to parse from nd2 file
        :return: the parsed metadata
        """
        main_part = info_to_parse[b'SLxImageTextInfo']
        metadata_text_dict = {}

        for key, value in main_part.items():
            # decode
            value_string = value.decode("utf-8")
            if value_string != '':
                # if not empty, split
                split_string = value_string.split('\r\n')
                for line_number, line in enumerate(split_string):
                    if line == '':
                        continue
                    if line_number == 0:  # these are the headers, they do not have a value, only key
                        key = key.decode("utf-8")
                        # we do not want those
                        if '5' in key or '6' in key:
                            continue
                        elif '9' in key:
                            # this is the date
                            value = value_string
                            key = 'date'
                        elif '13' in key:
                            # this is the objective
                            value = value_string
                            key = 'Objective'
                        # otherwise just add
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
    Self-made ND2 reader with improved metadata and fewer warnings
    """
    def __init__(self, filename, series=0, channel=0):
        self._clear_axes()
        self._get_frame_dict = dict()
        super().__init__(filename, series=series, channel=channel)

    def get_metadata(self, verbose=True):
        """
        Gets metadata of nd2 using self-made class.
        ---------------------------------------
        :param verbose: Verbose if you want all the outputs
        :return: metadata_dict: the dictionary of all metadata
        """
        metadata_nd2 = ND2ReaderForMetadata(self.filename)
        metadata_dict = metadata_nd2.get_metadata(verbose=verbose)
        metadata_nd2.close()

        return metadata_dict


class Roi:
    """
    ROI class. Used to determine region of interest
    """
    def __init__(self, x, y):
        """
        Initialization of ROI class.
        ---------------------------
        :param x: x position
        :param y: y position
        """
        self.x = x
        self.y = y
        self.index = None

        self.results = {}

    def set_index(self, index):
        """
        Sets index of ROI
        ----------------
        :param index: index to set
        """
        self.index = index

    def get_roi(self, frame, roi_size_1d, offset):
        """
        Gets ROI for a certain frame, offset, and ROI size
        ------------------------------------
        :param frame: frame to get ROI of
        :param roi_size_1d: ROI size
        :param offset: offset of current ROI in that frame
        :return: Something by Something ROI around the x/y position of this ROI
        """
        return frame[self.y + offset[0] - roi_size_1d:self.y + offset[0] + roi_size_1d + 1,
                     self.x + offset[1] - roi_size_1d:self.x + offset[1] + roi_size_1d + 1]

    def get_frame_stack(self, frames, roi_size_1d, offset):
        """
        Gets ROI for a certain frame stack, offset, and ROI size
        ------------------------------------
        :param frames: frames to get ROI of
        :param roi_size_1d: ROI size
        :param offset: offset of current ROI in that frame
        :return: Something by Something ROI around the x/y position of this ROI, in time
        """
        return frames[:, self.y + offset[0] - roi_size_1d:self.y + offset[0] + roi_size_1d + 1,
                      self.x + offset[1] - roi_size_1d:self.x + offset[1] + roi_size_1d + 1]

    def get_frame_stack_new(self, frames, roi_size_1d, offset):
        shape = np.zeros(2)
        shape[0] = frames.sizes['y']
        shape[1] = frames.sizes['x']

        return np.asarray(pims.process.crop(frames, ((self.y - roi_size_1d, shape[0] - self.y - roi_size_1d - 1),
                                                     (self.x - roi_size_1d, shape[1] - self.x - roi_size_1d - 1))))

    def get_frame_stack_new_tt(self, name, roi_size_1d, offset):
        frames = ND2ReaderSelf(name)
        shape = np.zeros(2)
        shape[0] = frames.sizes['y']
        shape[1] = frames.sizes['x']

        return np.asarray(pims.process.crop(frames, ((self.y - roi_size_1d, shape[0] - self.y - roi_size_1d - 1),
                                                     (self.x - roi_size_1d, shape[1] - self.x - roi_size_1d - 1))))

    def in_frame(self, shape, offset, margin):
        """
        Checks whether or not this ROI is in the frame
        --------------------------
        :param shape: Shape of frame
        :param offset: offset of frame
        :param margin: margin required to edge
        :return: in_frame_boolean: whether or not in frame
        """
        if self.x + offset[1] < margin or self.x + offset[1] > shape[1] - margin:
            in_frame_boolean = False
        elif self.y + offset[0] < margin or self.y + offset[0] > shape[0] - margin:
            in_frame_boolean = False
        else:
            in_frame_boolean = True

        return in_frame_boolean


import pims
import tracemalloc

ROI_SIZE = 7
ROI_SIZE_1D = (ROI_SIZE - 1) // 2

tt_name = "C:/Users/s150127/Downloads/___MBx/datasets/1nMimager_newGNRs_100mW.nd2"

loops = 10

tracemalloc.start()

nd2 = ND2ReaderSelf(tt_name)

test_rois = [Roi(i+20, i+40) for i in range(0, 61, 10)]

test_roi = Roi(50, 70)
import time
start_time = time.time()
for loop in range(loops):
    frame_stack_new = test_roi.get_frame_stack_new(nd2, ROI_SIZE_1D, [0, 0])
print("Time taken: {} s".format((time.time() - start_time) / 10))


from concurrent.futures import ThreadPoolExecutor
start_time = time.time()
for loop in range(loops):
    with ThreadPoolExecutor() as executor:
        frame_stacks_new = [executor.submit(roi.get_frame_stack_new_tt, tt_name, ROI_SIZE_1D, [0, 0]) for roi in test_rois]

print("Time taken: {} s".format((time.time() - start_time) / 10 / 7))

#nd2_array = np.asarray(nd2)

#frame_stack = test_roi.get_frame_stack(nd2_array, ROI_SIZE_1D, [0, 0])


current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage {current / 1e6}MB; Peak: {peak / 1e6}MB")
