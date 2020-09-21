import numpy as np


class Results:

    def __init__(self, created_by, roi_size, length):
        if created_by == "hsm":
            self.hsm_results = np.zeros(4)
            self.hsm_intensities = np.zeros(length)
            self.raw_data = np.zeros((length, roi_size, roi_size))
        elif created_by == "tt":
            self.tt: np.zeros(length)
            self.tt_drift: np.zeros(length)
            self.drift: np.zeros(length)
            self.raw_data = np.zeros((length, roi_size, roi_size))

    def to_dict(self):
        return self.__dict__


result_test = Results("hsm", 9, 10)

result_dict = result_test.to_dict()

test_result2 = dict(result_test)