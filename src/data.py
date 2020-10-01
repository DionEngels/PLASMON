class Roi:

    def __init__(self, x, y):

        self.x = x
        self.y = y
        self.index = None

        self.results = {}

    def set_index(self, index):
        self.index = index

    def get_roi(self, frame, roi_size_1d):
        return frame[self.y - roi_size_1d:self.y + roi_size_1d + 1,
                     self.x - roi_size_1d:self.x + roi_size_1d + 1]

    def in_frame(self, shape, offset):
        if self.x + offset[0] < 0 or self.x + offset[0] > shape[0]:
            in_frame_boolean = False
        elif self.y + offset[1] < 0 or self.y + offset[1] > shape[1]:
            in_frame_boolean = False
        else:
            in_frame_boolean = True

        return in_frame_boolean


