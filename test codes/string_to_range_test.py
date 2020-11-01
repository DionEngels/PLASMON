import numpy as np


def check_correct_chars(string):
    chars = set(string)

    for i in range(0, 10):
        chars.discard(str(i))

    chars.discard('[')
    chars.discard(']')
    chars.discard('.')
    chars.discard(',')
    chars.discard(':')
    chars.discard(' ')

    if len(chars) > 0:
        return False
    else:
        return True


def parse_string_to_numpy_array(string):
    wavelength_list = []
    string = string.replace('[', '')
    string = string.replace(']', '')

    split_string = string.split(',')

    for split in split_string:
        split = split.replace(' ','')
        try:
            split_int = int(split)
            wavelength_list.append(split_int)
        except ValueError:
            range_split = split.split(':')
            range_split = list(map(int, range_split))
            wavelength_list.extend(np.arange(range_split[0], range_split[2]+range_split[1], range_split[1]))

    return np.asarray(wavelength_list)


string = "[500:10:520, 532, 540:10:800]"

if check_correct_chars(string):
    wavelength = parse_string_to_numpy_array(string)
else:
    print("you done fd up")

