# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
MBx Python Data Analysis

v1.0: first version GUI: 13/07/2020
v1.1: Adding ROI size, wavelength to GUI: 13/07/2020
v1.2: Small additions, 7x7 ROI
v1.3: Disabling buttons, removing from grid
v1.4: Save and change wavelength per data set
v1.5: Allow for nm and pixel saving and MT
v1.6: conforming to standards of Python
v1.7: further cleanup
v1.8: complete reconstruction of labels, sliders, and buttons
v1.9: tweaked max sigma and max intensity
v1.10: histograms
v1.11: interactive histograms
v1.12: Status for MP
v1.13: removed wavelength
v2.0: new ROI finding method: 20/07/2020
v2.1: adding additional options for new ROI method
v2.2: new ROI method finished, new options finished, and new Fitter options
v2.3: Rejection options
v2.4: new grid
v2.5: styling
v2.6: styling ttk
v2.7: styling and warnings
v2.8: New method of destroy and updating

"""

# GENERAL IMPORTS
import os  # to get standard usage
from sys import exit
import time  # for timekeeping
from win32api import GetSystemMetrics  # Get sys info

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.patches as patches

# GUI
import tkinter as tk  # for GUI
from tkinter import ttk  # GUI styling
from tkinter.filedialog import askopenfilenames  # for popup that asks to select .nd2's

# ND2 related
from pims import ND2_Reader  # reader of ND2 files

# Own code
import _code.roi_finding as roi_finding
import _code.fitters as fitting
import _code.tools as tools

# Multiprocessing
import multiprocessing as mp

mpl.use("TkAgg")  # set back end to TK

# %% Initializations

FILETYPES = [("ND2", ".nd2")]

FONT_HEADER = "Verdana 14 bold"
FONT_SUBHEADER = "Verdana 11 bold"
FONT_STATUS = "Verdana 12"
FONT_BUTTON = "Verdana 9"
FONT_LABEL = "Verdana 10"
FONT_DROP = "Verdana 10"
FONT_BUTTON_BIG = "Verdana 20"
PAD_BIG = 30
PAD_SMALL = 10
INPUT_BIG = 25
INPUT_SMALL = 5

width = GetSystemMetrics(0)
height = GetSystemMetrics(1)
GUI_WIDTH = int(width * 0.70)
GUI_HEIGHT = int(height * 0.70)
GUI_WIDTH_START = int(width * 0.15)
GUI_HEIGHT_START = int(height * 0.15)
DPI = 100

# %% Options

fit_options = ["Gaussian - Fit bg", "Gaussian - Estimate bg",
               "Phasor + Intensity", "Phasor + Sum", "Phasor"]

rejection_options = ["Strict", "Loose", "None"]

roi_size_options = ["7x7", "9x9"]


# %% Multiprocessing main


def mt_main(name, fitter, frames_split, roi_locations, shared, q):
    nd2 = ND2_Reader(name)
    frames = nd2
    metadata = nd2.metadata
    metadata['sequence_count'] = len(frames_split)
    frames = frames[frames_split]

    local_result = fitter.main(frames, metadata, roi_locations,
                               q=q, start_frame=frames_split[0])

    for result_index, result in enumerate(local_result):
        shared[9 * result_index:9 * (result_index + 1)] = result[:]


# %% Own buttons / fields


class BigButton(ttk.Frame):
    def __init__(self, parent, height=None, width=None, text="", command=None, state='enabled'):
        ttk.Frame.__init__(self, parent, height=height, width=width)

        self.pack_propagate(0)
        self._btn = ttk.Button(self, text=text, command=command, state=state)
        self._btn.pack(fill=tk.BOTH, expand=1)
        self._btn["style"] = "Big.TButton"

    def updater(self, state='enabled'):
        self._btn['state'] = state


class FigureFrame(ttk.Frame):
    def __init__(self, parent, height=None, width=None, dpi=DPI):
        ttk.Frame.__init__(self, parent, height=height, width=width)

        self.pack_propagate(0)
        self.fig = Figure(figsize=(height / dpi, width / dpi), dpi=dpi)

        self.dpi = dpi
        self.width = width
        self.height = height

class EntryPlaceholder(ttk.Entry):
    def __init__(self, master=None, placeholder="PLACEHOLDER", *args, **kwargs):
        super().__init__(master, *args, style="Placeholder.TEntry", **kwargs)
        self.placeholder = placeholder

        self.insert("0", self.placeholder)
        self.bind("<FocusIn>", self._clear_placeholder)
        self.bind("<FocusOut>", self._add_placeholder)

    def _clear_placeholder(self, e):
        if self["style"] == "Placeholder.TEntry":
            self.delete("0", "end")
            self["style"] = "TEntry"

    def _add_placeholder(self, e):
        if not self.get():
            self.insert("0", self.placeholder)
            self["style"] = "Placeholder.TEntry"

    def updater(self, text=None):
        self.delete("0", "end")

        if text is None:
            self.insert("0", self.placeholder)
            self["style"] = "Placeholder.TEntry"
        else:
            self.insert("0", text)
            self["style"] = "TEntry"


class NormalButton:
    def __init__(self, parent, text=None, row=None, column=None,
                 rowspan=1, columnspan=1, command=None, state='enabled', sticky=None, padx=0, pady=0):
        self._btn = ttk.Button(parent, text=text, command=command, state=state)
        self.parent = parent
        self.text = text
        self.row = row
        self.column = column
        self.rowspan = rowspan
        self.columnspan = columnspan
        self.sticky = sticky
        self.padx = padx
        self.pady = pady
        self._btn.grid(row=row, column=column, rowspan=rowspan, columnspan=columnspan,
                       sticky=sticky, padx=padx, pady=pady)

    def updater(self, command=None, state='enabled', text=None):
        if text is None:
            text = self.text
        self._btn['text'] = text
        self._btn['state'] = state
        self._btn['command'] = command


class NormalSlider:
    def __init__(self, parent, from_=0, to=np.inf, resolution=1, start=0,
                 row=None, column=None, rowspan=1, columnspan=1, sticky=None, padx=0, pady=0):
        self._scale = tk.Scale(parent, from_=from_, to=to, orient='horizontal',
                               resolution=resolution, bg='white', borderwidth=0, highlightthickness=0)
        self._scale.set(start)
        self.parent = parent
        self.from_ = from_
        self.to = to
        self.resolution = resolution
        self.start = start
        self.row = row
        self.column = column
        self.rowspan = rowspan
        self.columnspan = columnspan
        self.sticky = sticky
        self.padx = padx
        self.pady = pady
        self._scale.grid(row=self.row, column=self.column, rowspan=self.rowspan, columnspan=self.columnspan,
                         sticky=sticky, padx=padx, pady=pady)

    def updater(self, from_=None, to=None, start=None):
        if from_ is None:
            from_ = self.from_
        if to is None:
            to = self.to
        if start is None:
            start = self.start
        self._scale.configure(from_=from_, to=to)
        self._scale.set(start)

        self.from_ = from_
        self.to = to
        self.start = start

    def get(self):
        return self._scale.get()


class NormalLabel:
    def __init__(self, parent, text=None, font=None, bd=None, relief=None,
                 row=None, column=None, rowspan=1, columnspan=1, sticky=None, padx=0, pady=0):
        self._label = tk.Label(parent, text=text, font=font, bd=bd, relief=relief, bg='white')
        self._label.grid(row=row, column=column, rowspan=rowspan, columnspan=columnspan,
                         sticky=sticky, padx=padx, pady=pady)
        self.parent = parent
        self.text = text
        self.font = font
        self.bd = bd
        self.relief = relief
        self.row = row
        self.column = column
        self.rowspan = rowspan
        self.columnspan = columnspan
        self.sticky = sticky
        self.padx = padx
        self.pady = pady

    def updater(self, text=None):
        if text is None:
            text = self.text
        self._label['text'] = text

# %% Container


class MbxPython(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # tk.Tk.iconbitmap(self, default = "")
        tk.Tk.wm_title(self, "MBx Python")
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        frame_tuple = (LoadPage, FittingPage)

        for to_load_frame in frame_tuple:
            frame = to_load_frame(container, self)
            self.frames[to_load_frame] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_load_frame(LoadPage)

    def show_load_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

    def show_frame(self, cont, filenames):
        frame = self.frames[cont]
        frame.update_init(self, cont, filenames)
        frame.tkraise()


# %% Quit

def quit_program():
    global gui
    gui.destroy()
    exit(0)


# %% Loading page

class LoadPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.configure(bg='white')

        button1 = BigButton(self, text="LOAD", height=int(GUI_HEIGHT / 4),
                            width=int(GUI_WIDTH / 4),  # style= 'Big.TButton',
                            command=lambda: self.load_nd2(controller))
        button1.grid(row=0, column=0)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

    def load_nd2(self, controller):
        global filenames
        tk.Tk().withdraw()
        filenames = askopenfilenames(filetypes=FILETYPES,
                                     title="Select file",
                                     initialdir=os.getcwd())

        if len(filenames) == 0:
            return

        controller.show_frame(FittingPage, filenames)


# %% Fitting page, initial update (after loading)

# noinspection PyBroadException
class FittingPage(tk.Frame):

    def update_init(self, parent, controller, filenames):

        self.__init__(parent, controller, reset=True)

        self.filenames = filenames
        self.dataset_index = 0

        self.dataset_roi_status.updater(text="Dataset " + str(self.dataset_index + 1) + " of " + str(len(filenames)))
        self.roi_status.updater(text="0 of " + str(len(filenames)) + " have settings")

        self.nd2 = ND2_Reader(filenames[self.dataset_index])
        self.frames = self.nd2
        self.metadata = self.nd2.metadata

        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(self.frames[0], extent=[0, self.frames[0].shape[1], self.frames[0].shape[0], 0],
                       aspect='auto')

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=40, columnspan=10, rowspan=12, sticky='EW')

        self.restore_default()

        if self.dataset_index == len(filenames) - 1:
            self.button_right.updater(state='disabled')
        else:
            self.button_right.updater(command=lambda: self.next_dataset())

        if self.dataset_index == 0:
            self.button_left.updater(state='disabled')
        else:
            self.button_left.updater(command=lambda: self.previous_dataset())

    # %% Fitting page, fit ROIs

    def fit_rois(self):

        int_min = self.min_int_slider.get()
        int_max = self.max_int_slider.get()
        sigma_min = self.min_sigma_slider.get()
        sigma_max = self.max_sigma_slider.get()
        corr_min = self.min_corr_slider.get()
        pixel_min = self.min_pixel_slider.get()
        roi_size = int(self.roi_var.get()[0])

        filter_size = int(self.filter_input.get())
        roi_side = int(self.roi_side_input.get())
        inter_roi = int(self.inter_roi_input.get())

        if roi_side < int((roi_size - 1) / 2):
            tk.messagebox.showerror("ERROR", "Distance to size cannot be smaller than 1D ROI size")
            return False
        if filter_size % 2 != 1:
            tk.messagebox.showerror("ERROR", "Filter size should be odd")
            return False

        self.roi_fitter = fitting.Gaussian(roi_size, {}, "None", "Gaussian", 5)

        self.roi_finder.change_settings(int_min, int_max, corr_min, pixel_min,
                                        sigma_min, sigma_max, roi_size,
                                        filter_size, roi_side, inter_roi)

        self.temp_roi_locations = self.roi_finder.main(self.roi_fitter)

        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)

        fig_sub.imshow(self.frames[0], extent=[0, self.frames[0].shape[1], self.frames[0].shape[0], 0],
                       aspect='auto')

        roi_locations_temp = self.temp_roi_locations - self.roi_finder.roi_size_1d

        for roi in roi_locations_temp:
            rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                     linewidth=0.5, edgecolor='r', facecolor='none')
            fig_sub.add_patch(rect)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=40, columnspan=10, rowspan=12, sticky='EW')

        return True

    # %% Fitting page, save ROI settings

    def save_roi_settings(self):

        success = self.fit_rois()

        if not success:
            return

        self.roi_locations[self.dataset_index] = self.temp_roi_locations

        int_min = self.min_int_slider.get()
        int_max = self.max_int_slider.get()
        sigma_min = self.min_sigma_slider.get()
        sigma_max = self.max_sigma_slider.get()
        corr_min = self.min_corr_slider.get()
        pixel_min = self.min_pixel_slider.get()
        roi_size = int(self.roi_var.get()[0])

        filter_size = int(self.filter_input.get())
        roi_side = int(self.roi_side_input.get())
        inter_roi = int(self.inter_roi_input.get())

        settings = {'int_max': int_max, 'int_min': int_min,
                    'sigma_min': sigma_min, 'sigma_max': sigma_max,
                    'corr_min': corr_min, 'pixel_min': pixel_min,
                    'roi_size': roi_size, 'filter_size': filter_size,
                    'roi_side': roi_side, 'inter_roi': inter_roi}

        self.saved_settings[self.dataset_index] = settings

        self.roi_status.updater(text=str(len(self.roi_locations)) + " of " + str(len(
            self.filenames)) + " have settings")
        self.button_restore_saved.updater(command=lambda: self.restore_saved())

        if len(self.roi_locations) == len(self.filenames):
            self.button_fit.updater(state='enabled')

    # %% Fitting page, restore default settings

    def restore_default(self):

        self.roi_fitter = fitting.Gaussian(7, {}, "None", "Gaussian", 5)

        self.roi_finder = roi_finding.RoiFinder(9, self.frames[0], self.roi_fitter)

        self.min_int_slider.updater(from_=0, to=self.roi_finder.int_max / 4,
                                    start=self.roi_finder.int_min)
        self.max_int_slider.updater(from_=0, to=self.roi_finder.int_max,
                                    start=self.roi_finder.int_max)
        self.min_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=self.roi_finder.sigma_min)
        self.max_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=self.roi_finder.sigma_max)
        self.min_corr_slider.updater(from_=0, to=1, start=self.roi_finder.corr_min)
        self.min_pixel_slider.updater(from_=0, to=np.max(self.frames[0]),
                                      start=self.roi_finder.pixel_min)
        self.roi_var.set(roi_size_options[0])

        self.filter_input.updater()
        self.roi_side_input.updater()
        self.inter_roi_input.updater()

    # %% Fitting page, switch between datasets

    def next_dataset(self):

        self.temp_roi_locations = None
        self.histogram = None
        self.to_hist = None

        self.dataset_index += 1
        self.dataset_roi_status.updater(text="Dataset " + str(self.dataset_index + 1) + " of "
                                             + str(len(self.filenames)))

        self.nd2.close()
        self.nd2 = ND2_Reader(self.filenames[self.dataset_index])
        self.frames = self.nd2
        self.metadata = self.nd2.metadata

        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(self.frames[0], extent=[0, self.frames[0].shape[1], self.frames[0].shape[0], 0],
                       aspect='auto')

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=40, columnspan=10, rowspan=12, sticky='EW')

        self.restore_default()

        if self.dataset_index == len(self.filenames) - 1:
            self.button_right.updater(state='disabled')
        else:
            self.button_right.updater(command=lambda: self.next_dataset())

        if self.dataset_index == 0:
            self.button_left.updater(state='disabled')
        else:
            self.button_left.updater(command=lambda: self.previous_dataset())

        if self.dataset_index in self.saved_settings:
            self.button_restore_saved.updater(command=lambda: self.restore_saved())
        else:
            self.button_restore_saved.updater(state='disabled')

    def previous_dataset(self):

        self.temp_roi_locations = None
        self.histogram = None
        self.to_hist = None

        self.dataset_index -= 1
        self.dataset_roi_status.updater(text="Dataset " +
                                             str(self.dataset_index + 1) + " of "
                                             + str(len(self.filenames)))

        self.nd2.close()
        self.nd2 = ND2_Reader(filenames[self.dataset_index])
        self.frames = self.nd2
        self.metadata = self.nd2.metadata

        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(self.frames[0], extent=[0, self.frames[0].shape[1], self.frames[0].shape[0], 0],
                       aspect='auto')

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=40, columnspan=10, rowspan=12, sticky='EW')

        self.restore_default()

        if self.dataset_index == len(self.filenames) - 1:
            self.button_right.updater(state='disabled')
        else:
            self.button_right.updater(command=lambda: self.next_dataset())

        if self.dataset_index == 0:
            self.button_left.updater(state='disabled')
        else:
            self.button_left.updater(command=lambda: self.previous_dataset())

        if self.dataset_index in self.saved_settings:
            self.button_restore_saved.updater(command=lambda: self.restore_saved())
        else:
            self.button_restore_saved.updater(state='disabled')

    # %% Fitting page, start fitting

    def start_fitting(self):

        if len(self.roi_locations) != len(self.filenames):
            tk.messagebox.showerror("ERROR", "Not all datasets have settings yet, cannot start")
            return

        check = tk.messagebox.askokcancel("Are you sure?",
                                          "Fitting may take a while. Are you sure everything is set up correctly?")
        if not check:
            return

        n_processes = self.cores_var.get()
        basedir = os.getcwd()
        start_frame = self.frame_begin_input.get()
        end_frame = self.frame_end_input.get()
        method = self.method_var.get()
        rejection_type = self.rejection_var.get()

        if rejection_type != "None" and (method == "Phasor + Sum" or method == "Phasor"):
            rejection_check = tk.messagebox.askokcancel("Just a heads up",
                                                        """These Phasor methods have no rejection options,
            so "None" rejection will be used""")
            if not rejection_check:
                return

        if rejection_type == "Loose" and method == "Phasor + Intensity":
            rejection_check = tk.messagebox.askokcancel("Just a heads up",
                                                        """This Phasor method can only do "Strict" or
            "None" rejection. "Loose" will be changed to "None".""")
            rejection_type = "None"
            if not rejection_check:
                return

        if n_processes > 1 and (method == "Phasor + Intensity" or method == "Phasor + Sum" or method == "Phasor"):
            cores_check = tk.messagebox.askokcancel("Just a heads up",
                                                    """Phasor will be used with one core since the
            overhead only slows it down""")
            if not cores_check:
                return
            n_processes = 1

        directory = ""
        metadata_filtered = {}
        roi_locations_dict = {}
        results_dict = {}
        self.start_time = time.time()
        results_counter = 0

        for dataset_index, filename in enumerate(self.filenames):
            self.nd2.close()
            self.nd2 = ND2_Reader(filenames[self.dataset_index])
            self.frames = self.nd2
            self.metadata = self.nd2.metadata
            self.dataset_index = dataset_index

            roi_size = self.saved_settings[dataset_index]['roi_size']

            roi_locations = self.roi_locations[dataset_index]

            if method == "Phasor + Intensity":
                fitter = fitting.Phasor(roi_size, self.saved_settings[dataset_index], rejection_type, method)
            elif method == "Phasor":
                fitter = fitting.PhasorDumb(roi_size, self.saved_settings[dataset_index], rejection_type, method)
            elif method == "Gaussian - Fit bg":
                fitter = fitting.GaussianBackground(roi_size, self.saved_settings[dataset_index], rejection_type,
                                                    method, 6)
            elif method == "Gaussian - Estimate bg":
                fitter = fitting.Gaussian(roi_size, self.saved_settings[dataset_index], rejection_type, method, 5)
            else:
                fitter = fitting.PhasorSum(roi_size, self.saved_settings[dataset_index], rejection_type, method)

            if start_frame == "Leave empty for start" and end_frame == "Leave empty for end":
                end_frame = self.metadata['sequence_count']
                start_frame = 0
                to_fit = self.frames
            elif start_frame == "Leave empty for start" and end_frame != "Leave empty for end":
                start_frame = 0
                end_frame = int(end_frame)
                to_fit = self.frames[:end_frame]
            elif start_frame != "Leave empty for start" and end_frame == "Leave empty for end":
                start_frame = int(start_frame)
                end_frame = self.metadata['sequence_count']
                to_fit = self.frames[start_frame:]
            else:  # start_frame != "Leave empty for start" and end_frame != "Leave empty for end":
                start_frame = int(start_frame)
                end_frame = int(end_frame)
                to_fit = self.frames[start_frame:end_frame]

            num_frames = end_frame - start_frame

            if n_processes > 1:
                shared = [None] * n_processes
                results = np.zeros((len(roi_locations) * len(self.frames), 9))

                frames_split = np.array_split(list(range(start_frame, end_frame)), n_processes)
                processes = [None] * n_processes
                q = mp.Queue()

                self.update_status(0, num_frames)

                for i in range(0, n_processes):
                    shared[i] = mp.Array('d', int(9 * len(roi_locations) * len(frames_split[i])))
                    processes[i] = (mp.Process(target=mt_main,
                                               args=(filename, fitter,
                                                     frames_split[i], roi_locations,
                                                     shared[i], q)))

                for p in processes:
                    p.start()

                queue_dict = {}
                frames_fitted = 0
                while frames_fitted < num_frames * 0.95:
                    for p in processes:
                        queue = q.get()
                        queue_dict[queue[0]] = queue[1]
                        frames_fitted = sum(queue_dict.values())
                        self.update_status(frames_fitted, num_frames)

                for p in processes:
                    p.join()

                self.update_status(num_frames, num_frames)

                counter = 0
                for i, share in enumerate(shared):
                    arr = np.frombuffer(share.get_obj())  # mp_arr and arr share the same memory
                    result = arr.reshape((len(roi_locations) * len(frames_split[i]), 9))
                    results[counter:counter + len(roi_locations) * len(frames_split[i])] = result
                    counter += len(roi_locations) * len(frames_split[i])

                results = results[results[:, 3] != 0]
            else:
                self.update_status(0, num_frames)
                results = fitter.main(to_fit, self.metadata, roi_locations,
                                      gui=self, n_frames=num_frames)

            nm_or_pixels = self.dimension.get()
            if nm_or_pixels == "nm":
                pixelsize_nm = self.metadata['calibration_um'] * 1000
                results[:, 2] = results[:, 2] * pixelsize_nm
                results[:, 3] = results[:, 3] * pixelsize_nm

            # create folder for output

            directory = filename.split(".")[0].split("/")[-1]
            path = os.path.join(basedir, directory)
            try:
                os.mkdir(path)
            except:
                pass

            # Save

            metadata_filtered = {k: v for k, v in self.metadata.items() if v is not None}
            del metadata_filtered['time_start']
            del metadata_filtered['time_start_utc']

            # ROI_locations dict
            roi_locations_dict = dict(zip(['x', 'y'], roi_locations.T))

            # Localization dict
            results_dict = {'Localizations': results}
            results_counter += results.shape[0]

        end_message = 'Time taken: ' + str(round(time.time() - self.start_time, 3)) \
                      + ' s. Fits done: ' + str(results_counter)

        tk.messagebox.showinfo("Done!", end_message)

        tools.save_to_csv_mat('metadata', metadata_filtered, directory)
        tools.save_to_csv_mat('ROI_locations', roi_locations_dict, directory)
        tools.save_to_csv_mat('Localizations', results_dict, directory)

        self.nd2.close()

    # %% Fitting page, update the status

    def update_status(self, progress, comparator):

        num_files = len(self.filenames)
        base_progress = self.dataset_index / num_files * 100

        file_progress = progress / comparator * 100 / num_files
        progress = base_progress + file_progress

        progress_text = str(round(progress, 1)) + "% done"

        if progress == 0:
            self.time_status_label.updater(text="TBD")
        else:
            current_time = time.time()
            time_taken = current_time - self.start_time
            time_done_estimate = time_taken * 100 / progress + self.start_time
            tr = time.localtime(time_done_estimate)
            time_text = str(tr[3]) + ":" + str(tr[4]) + ":" + str(tr[5]) + " " + str(tr[2]) + "/" + str(tr[1])

            self.time_status_label.updater(text=time_text)
        self.progress_status_label.updater(text=progress_text)

        self.update()

    # %% Fitting page, return to load page

    def load_new(self, controller):

        controller.show_load_frame(LoadPage)

    # %% Fitting page, restore saved settings

    def restore_saved(self):

        settings = self.saved_settings[self.dataset_index]

        self.min_int_slider.updater(from_=0, to=self.roi_finder.int_max / 4,
                                    start=settings['int_min'])
        self.max_int_slider.updater(from_=0, to=self.roi_finder.int_max,
                                    start=settings['int_max'])
        self.min_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=settings['sigma_min'])
        self.max_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=settings['sigma_max'])
        self.min_corr_slider.updater(from_=0, to=1, start=settings['corr_min'])
        self.min_pixel_slider.updater(from_=0, to=np.max(self.frames[0]),
                                      start=settings['pixel_min'])
        if settings['roi_size'] == 7:
            self.roi_var.set(roi_size_options[0])
        else:
            self.roi_var.set(roi_size_options[1])

        self.filter_input.updater(settings['filter_size'])
        self.roi_side_input.updater(settings['roi_side'])
        self.inter_roi_input.updater(settings['inter_roi'])

    # %% Histogram of sliders

    def make_histogram(self, variable):

        fig_sub = self.histogram.add_subplot(111)
        hist, bins, _ = fig_sub.hist(self.to_hist, bins='auto')

        min_int = self.min_int_slider.get()
        max_int = self.max_int_slider.get()
        min_sigma = self.min_sigma_slider.get()
        max_sigma = self.max_sigma_slider.get()
        min_peak = self.min_pixel_slider.get()
        min_corr = self.min_corr_slider.get()

        if variable == "min_int" or variable == "max_int":
            self.histogram.clear()
            logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
            plt.hist(self.to_hist, bins=logbins)
            plt.title("Intensity. Use graph select to change threshold")
            plt.axvline(x=min_int, color='red', linestyle='--')
            plt.axvline(x=max_int, color='red', linestyle='--')
            plt.xscale('log')
        elif variable == "min_sigma" or variable == "max_sigma":
            plt.title("Sigma. Use graph select to change threshold")
            plt.axvline(x=min_sigma, color='red', linestyle='--')
            plt.axvline(x=max_sigma, color='red', linestyle='--')
        elif variable == "peak_min":
            self.histogram.clear()
            logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
            plt.hist(self.to_hist, bins=logbins)
            plt.title("Minimum pixel value for peak. Use graph select to change threshold")
            plt.axvline(x=min_peak, color='red', linestyle='--')
            plt.xscale('log')
        else:
            self.histogram.clear()
            logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
            plt.hist(self.to_hist, bins=logbins)
            plt.title("Correlation values for ROIs. Use graph select to change threshold")
            plt.axvline(x=min_corr, color='red', linestyle='--')
            plt.xscale('log')

        plt.show()

    def fun_histogram(self, variable):

        if variable == "min_int" or variable == "max_int":
            self.to_hist = self.roi_finder.main(self.roi_fitter, return_int=True)
        elif variable == "peak_min":
            self.to_hist = np.ravel(self.frames[0])
        elif variable == "corr_min":
            self.to_hist = self.roi_finder.main(self.roi_fitter, return_corr=True)
        else:
            self.to_hist = self.roi_finder.main(self.roi_fitter, return_sigmas=True)

        self.histogram = plt.figure(figsize=(6.4 * 1.2, 4.8 * 1.2))

        self.make_histogram(variable)

    def histogram_select(self, variable):
        click = self.histogram.ginput(1)

        if variable == "min_int":
            self.min_int_slider.updater(start=int(click[0][0]))
        elif variable == 'max_int':
            self.max_int_slider.updater(start=int(click[0][0]))
        elif variable == "min_sigma":
            self.min_sigma_slider.updater(start=click[0][0])
        elif variable == "max_sigma":
            self.max_sigma_slider.updater(start=click[0][0])
        elif variable == "peak_min":
            self.min_pixel_slider.updater(start=click[0][0])
        else:
            self.min_corr_slider.updater(start=click[0][0])

        self.histogram.clear()
        self.make_histogram(variable)

    # %% Figure popout

    def figure_popout(self):

        popout = plt.figure(figsize=(GUI_WIDTH / 100, GUI_HEIGHT / 100), dpi=100)
        fig_sub = popout.add_subplot(111)

        fig_sub.imshow(self.frames[0], extent=[0, self.frames[0].shape[1], self.frames[0].shape[0], 0],
                       aspect='auto')

        try:
            roi_locations_temp = self.temp_roi_locations - self.roi_finder.roi_size_1d
            roi_size = self.roi_finder.roi_size

            for roi in roi_locations_temp:
                rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                         linewidth=0.5, edgecolor='r', facecolor='none')
                fig_sub.add_patch(rect)
        except:
            pass  # only when temp_roi_locations not yet defined, before first fit

        plt.show()

    # %% Fitting page, initial declaration

    def __init__(self, parent, controller, reset=False):

        if not reset:
            tk.Frame.__init__(self, parent)
            self.configure(bg='white')
        self.nd2 = None
        self.frames = None
        self.roi_finder = None
        self.metadata = None
        self.temp_roi_locations = None
        self.filenames = None
        self.start_time = None
        self.roi_locations = {}
        self.dataset_index = 0
        self.saved_settings = {}
        self.histogram = None
        self.to_hist = None
        self.roi_fitter = None

        roi_finding_label = tk.Label(self, text="ROI finding", font=FONT_HEADER, bg='white')
        roi_finding_label.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_corr_label = tk.Label(self, text="Minimum Correlation", font=FONT_LABEL, bg='white')
        min_corr_label.grid(row=0, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_corr_slider = NormalSlider(self, from_=0, to=1, resolution=0.005,
                                            row=1, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_corr_histogram = ttk.Button(self, text="Graph",
                                        command=lambda: self.fun_histogram("corr_min"))
        min_corr_histogram.grid(row=1, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_corr_histogram_select = ttk.Button(self, text="Graph select",
                                               command=lambda: self.histogram_select("corr_min"))
        min_corr_histogram_select.grid(row=2, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_pixel_label = tk.Label(self, text="Minimum pixel intensity", font=FONT_LABEL, bg='white')
        min_pixel_label.grid(row=0, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_pixel_slider = NormalSlider(self, from_=0, to=5000,
                                             row=1, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_pixel_histogram = ttk.Button(self, text="Graph",
                                         command=lambda: self.fun_histogram("peak_min"))
        min_pixel_histogram.grid(row=1, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_pixel_histogram_select = ttk.Button(self, text="Graph select",
                                                command=lambda: self.histogram_select("peak_min"))
        min_pixel_histogram_select.grid(row=2, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_int_label = tk.Label(self, text="Minimum Intensity", font=FONT_LABEL, bg='white')
        min_int_label.grid(row=3, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_int_slider = NormalSlider(self, from_=0, to=1000,
                                           row=4, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_int_histogram = ttk.Button(self, text="Graph",
                                       command=lambda: self.fun_histogram("min_int"))
        min_int_histogram.grid(row=4, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_int_histogram_select = ttk.Button(self, text="Select min",
                                              command=lambda: self.histogram_select("min_int"))
        min_int_histogram_select.grid(row=4, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_int_label = tk.Label(self, text="Maximum Intensity", font=FONT_LABEL, bg='white')
        max_int_label.grid(row=3, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.max_int_slider = NormalSlider(self, from_=0, to=5000,
                                           row=4, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_int_histogram_select = ttk.Button(self, text="Select max",
                                              command=lambda: self.histogram_select("max_int"))
        max_int_histogram_select.grid(row=4, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_sigma_label = tk.Label(self, text="Minimum Sigma", font=FONT_LABEL, bg='white')
        min_sigma_label.grid(row=5, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_sigma_slider = NormalSlider(self, from_=0, to=5, resolution=0.01,
                                             row=6, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_sigma_histogram = ttk.Button(self, text="Graph",
                                         command=lambda: self.fun_histogram("min_sigma"))
        min_sigma_histogram.grid(row=6, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_sigma_histogram_select = ttk.Button(self, text="Select min",
                                                command=lambda: self.histogram_select("min_sigma"))
        min_sigma_histogram_select.grid(row=6, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_sigma_label = tk.Label(self, text="Maximum Sigma", font=FONT_LABEL, bg='white')
        max_sigma_label.grid(row=5, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.max_sigma_slider = NormalSlider(self, from_=0, to=10, resolution=0.01,
                                             row=6, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_sigma_histogram_select = ttk.Button(self, text="Select max",
                                                command=lambda: self.histogram_select("max_sigma"))
        max_sigma_histogram_select.grid(row=6, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=7, column=0, rowspan=1, columnspan=40, sticky='we')

        advanced_label = tk.Label(self, text="Advanced settings", font=FONT_SUBHEADER, bg='white')
        advanced_label.grid(row=8, column=0, columnspan=40, sticky='EW', padx=PAD_SMALL)

        roi_size_label = tk.Label(self, text="ROI size", bg='white', font=FONT_LABEL)
        roi_size_label.grid(row=9, column=0, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.roi_var = tk.StringVar(self)
        # self.roi_var.set(roi_size_options[0])

        roi_drop = ttk.OptionMenu(self, self.roi_var, roi_size_options[0], *roi_size_options)
        roi_drop.grid(row=9, column=5, columnspan=5, sticky='EW', padx=PAD_SMALL)

        filter_label = tk.Label(self, text="Filter size", bg='white', font=FONT_LABEL)
        filter_label.grid(row=9, column=10, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.filter_input = EntryPlaceholder(self, "9", width=INPUT_SMALL)
        self.filter_input.grid(row=9, column=15, columnspan=5)

        roi_side_label = tk.Label(self, text="Spacing side", bg='white', font=FONT_LABEL)
        roi_side_label.grid(row=9, column=20, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.roi_side_input = EntryPlaceholder(self, "11", width=INPUT_SMALL)
        self.roi_side_input.grid(row=9, column=25, columnspan=5)

        inter_roi_label = tk.Label(self, text="Inter-ROI spacing", bg='white', font=FONT_LABEL)
        inter_roi_label.grid(row=9, column=30, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.inter_roi_input = EntryPlaceholder(self, "6", width=INPUT_SMALL)
        self.inter_roi_input.grid(row=9, column=35, columnspan=5)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=10, column=0, rowspan=1, columnspan=40, sticky='we')

        self.button_left = NormalButton(self, text="<<", row=11, column=0, columnspan=5, sticky='EW',
                                        padx=PAD_SMALL)
        self.dataset_roi_status = NormalLabel(self, text="TBD",
                                              row=11, column=5, columnspan=30, font=FONT_STATUS)
        self.button_right = NormalButton(self, ">>", row=11, column=35, columnspan=5, sticky='EW',
                                         padx=PAD_SMALL)

        button_fit = ttk.Button(self, text="Fit", command=lambda: self.fit_rois())
        button_fit.grid(row=12, column=0, columnspan=10, sticky='EW', padx=PAD_BIG)

        button_restore = ttk.Button(self, text="Restore default",
                                    command=lambda: self.restore_default())
        button_restore.grid(row=12, column=10, columnspan=10, sticky='EW', padx=PAD_BIG)

        self.button_restore_saved = NormalButton(self, text="Restore saved",
                                                 state='disabled',
                                                 command=lambda: self.restore_saved(),
                                                 row=12, column=20,
                                                 columnspan=10, sticky='EW', padx=PAD_BIG)

        button_save = ttk.Button(self, text="Save",
                                 command=lambda: self.save_roi_settings())
        button_save.grid(row=12, column=30, columnspan=10, sticky='EW', padx=PAD_BIG)

        self.fig = Figure(figsize=(GUI_WIDTH / DPI * 0.4, GUI_WIDTH / DPI * 0.4), dpi=DPI)
        #  self.fig = FigureFrame(self, height=GUI_WIDTH*0.4, width=GUI_WIDTH*0.4, dpi=DPI)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=40, columnspan=10, rowspan=12, sticky='EW',
                                         padx=PAD_SMALL)

        button_popout = ttk.Button(self, text="Pop-out to zoom",
                                   command=lambda: self.figure_popout())
        button_popout.grid(row=12, column=40, columnspan=10)  # , sticky='EW', padx=PAD_BIG)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=18, column=0, rowspan=2, columnspan=50, sticky='we')

        fit_area_label = tk.Label(self, text="Fitting", font=FONT_HEADER, bg='white')
        fit_area_label.grid(row=21, column=0, columnspan=10, sticky='EW', padx=PAD_SMALL)

        self.roi_status = NormalLabel(self, text="TBD",
                                      row=21, column=10, columnspan=30, font=FONT_STATUS)

        method_label = tk.Label(self, text="Method", font=FONT_LABEL, bg='white')
        method_label.grid(row=23, column=0, columnspan=10, sticky='EW', padx=PAD_SMALL)

        self.method_var = tk.StringVar(self)
        # self.method_var.set(fit_options[1])

        method_drop = ttk.OptionMenu(self, self.method_var, fit_options[1], *fit_options)
        method_drop.grid(row=24, column=0, columnspan=10, sticky="ew")

        rejection_label = tk.Label(self, text="ROI size", bg='white', font=FONT_LABEL)
        rejection_label.grid(row=23, column=10, columnspan=10, sticky='EW', padx=PAD_SMALL)

        self.rejection_var = tk.StringVar(self)
        # self.rejection_var.set(rejection_options[1])

        rejection_drop = ttk.OptionMenu(self, self.rejection_var, rejection_options[1], *rejection_options)
        rejection_drop.grid(row=24, column=10, columnspan=10, sticky='EW', padx=PAD_SMALL)

        cores_label = tk.Label(self, text="#cores", font=FONT_LABEL, bg='white')
        cores_label.grid(row=23, column=20, columnspan=10, sticky='EW', padx=PAD_BIG)

        self.total_cores = mp.cpu_count()
        cores_options = [1, int(self.total_cores / 2), int(self.total_cores * 3 / 4), int(self.total_cores)]
        self.cores_var = tk.IntVar(self)
        # self.cores_var.set(cores_options[0])
        self.cores_drop = ttk.OptionMenu(self, self.cores_var, cores_options[0], *cores_options)
        self.cores_drop.grid(row=24, column=20, columnspan=10, sticky='EW', padx=PAD_BIG)

        dimensions_label = tk.Label(self, text="pixels or nm", font=FONT_LABEL, bg='white')
        dimensions_label.grid(row=23, column=30, columnspan=10, sticky='EW', padx=PAD_BIG)

        dimension_options = ["nm", "pixels"]
        self.dimension = tk.StringVar(self)
        # self.dimension.set(dimension_options[0])
        self.dimension_drop = ttk.OptionMenu(self, self.dimension, dimension_options[0], *dimension_options)
        self.dimension_drop.grid(row=24, column=30, columnspan=10, sticky='EW', padx=PAD_BIG)

        frame_begin_label = tk.Label(self, text="Begin frame", font=FONT_LABEL, bg='white')
        frame_begin_label.grid(row=27, column=0, columnspan=20, sticky='EW', padx=PAD_BIG)

        self.frame_begin_input = EntryPlaceholder(self, "Leave empty for start", width=INPUT_BIG)
        self.frame_begin_input.grid(row=28, column=0, columnspan=20)

        frame_end_label = tk.Label(self, text="End frame", font=FONT_LABEL, bg='white')
        frame_end_label.grid(row=27, column=20, columnspan=20, sticky='EW', padx=PAD_BIG)

        self.frame_end_input = EntryPlaceholder(self, "Leave empty for end", width=INPUT_BIG)
        self.frame_end_input.grid(row=28, column=20, columnspan=20)

        self.button_fit = BigButton(self, text="FIT", height=int(GUI_HEIGHT / 8),
                                    width=int(GUI_WIDTH / 8), state='disabled',
                                    command=lambda: self.start_fitting())  # , style= 'my.TButton')
        self.button_fit.grid(row=23, column=40, columnspan=5, rowspan=5)

        progress_label = tk.Label(self, text="Progress", font=FONT_LABEL, bg='white')
        progress_label.grid(row=23, column=45, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.progress_status_label = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                 row=24, column=45, columnspan=5, sticky="ew", font=FONT_LABEL)

        time_label = tk.Label(self, text="Estimated time done", font=FONT_LABEL, bg='white')
        time_label.grid(row=26, column=45, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.time_status_label = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                             row=27, column=45, columnspan=5, sticky="ew", font=FONT_LABEL)

        button_load_new = ttk.Button(self, text="Load new", command=lambda: self.load_new(parent))
        button_load_new.grid(row=50, column=45, columnspan=2, sticky='EW', padx=PAD_SMALL)

        button_quit = ttk.Button(self, text="Quit", command=quit_program)
        button_quit.grid(row=50, column=48, columnspan=2, sticky='EW', padx=PAD_SMALL)

        for i in range(0, 49):
            self.columnconfigure(i, weight=1, minsize=15)

        for i in range(0, 51):
            self.rowconfigure(i, weight=1)


# %% START GUI
if __name__ == '__main__':
    gui = MbxPython()
    gui.geometry(str(GUI_WIDTH) + "x" + str(GUI_HEIGHT) + "+" + str(GUI_WIDTH_START) + "+" + str(GUI_HEIGHT_START))

    ttk_style = ttk.Style(gui)
    ttk_style.configure("Big.TButton", font=FONT_BUTTON_BIG)
    ttk_style.configure("Placeholder.TEntry", foreground="Grey")
    ttk_style.configure("TButton", font=FONT_BUTTON, background="Grey")
    ttk_style.configure("TSeparator", background="black")
    ttk_style.configure("TMenubutton", font=FONT_DROP, background="White")

    #  print(ttk_style.layout("TMenubutton"))
    #  print(ttk_style.element_options('Menubutton.dropdown'))

    gui.mainloop()
