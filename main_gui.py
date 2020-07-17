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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.patches as patches

# GUI
import tkinter as tk  # for GUI
from tkinter import ttk  # GUI styling
from tkinter.filedialog import askopenfilenames  # for popup that asks to select .nd2's

# ND2 related
from pims import ND2_Reader  # reader of ND2 files
# from pims import Frame # for converting ND2 generator to one numpy array

# Own code
import _code.roi_finding as roi_finding
import _code.fitters as fitting
import _code.tools as tools

# Multiprocessing
import multiprocessing as mp

mpl.use("TkAgg")  # set back end to TK

# %% Inputs
ROI_SIZE = 9
THRESHOLD = 5  # X*sigma

# %% Initializations

FILETYPES = [("ND2", ".nd2")]

LARGE_FONT = ("Verdana", 12)
# style.use("ggplot")

fit_options = ["Phasor with intensity", "Phasor without intensity",
               "Gaussian with background", "Gaussian without background"]

roi_size_options = ["7x7", "9x9"]

width = GetSystemMetrics(0)
height = GetSystemMetrics(1)
GUI_WIDTH = int(width / 2)
GUI_HEIGHT = int(height / 2)
GUI_WIDTH_START = int(width / 4)
GUI_HEIGHT_START = int(height / 4)
DPI = 100

# %% Multiprocessing main


def mt_main(name, fitter, frames_split, roi_locations, shared):
    nd2 = ND2_Reader(name)
    frames = nd2
    metadata = nd2.metadata
    metadata['sequence_count'] = len(frames_split)
    frames = frames[frames_split]

    local_result = fitter.main(frames, metadata, roi_locations)
    local_result[:, 0] += frames_split[0]

    for result_index, result in enumerate(local_result):
        shared[9 * result_index:9 * (result_index + 1)] = result[:]

# %% Own buttons / fields


class BigButton(ttk.Frame):
    def __init__(self, parent, height=None, width=None, text="", command=None, state='enabled'):
        ttk.Frame.__init__(self, parent, height=height, width=width)

        self.pack_propagate(0)
        self._btn = ttk.Button(self, text=text, command=command, state=state)
        self._btn.pack(fill=tk.BOTH, expand=1)
            
            
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
        self["style"] = "TEntry"
        self.insert("0", text)

class NormalButton(ttk.Button):
    def __init__(self, parent, text=None, row=None, column=None, 
                 rowspan=1, columnspan=1, command=None, state='enabled'):
        super().__init__(parent, text=text, command=command, state=state)
        self.parent = parent
        self.text = text
        self.row = row
        self.column = column
        self.rowspan = rowspan
        self.columnspan = columnspan
        super().grid(row=row, column=column, 
                  rowspan=rowspan, columnspan=columnspan)
        
    def updater(self, command=None, state='enabled', text=None):
        if text == None:
            text = self.text
        super().grid_forget()
        super().__init__(self.parent, text=text, command=command, state=state)
        super().grid(row=self.row, column=self.column, rowspan=self.rowspan, 
                  columnspan=self.columnspan)
        
class NormalSlider(tk.Scale):
    def __init__(self, parent, from_=0, to=np.inf, resolution=1, start=0,
                 row=None, column=None, rowspan=1, columnspan=1):
        super().__init__(parent, from_=from_, to=to, orient='horizontal',
                       resolution=resolution)
        super().set(start)
        self.parent = parent
        self.from_ = from_
        self.to = to
        self.resolution = resolution
        self.start = start
        self.row = row
        self.column = column
        self.rowspan = rowspan
        self.columnspan = columnspan
        super().grid(row=self.row, column=self.column, rowspan=self.rowspan, 
                  columnspan=self.columnspan)
        
    def updater(self, from_=None, to=None, start=None):
        if from_ == None:
            from_ = self.from_
        if to == None:
            to = self.to
        if start == None:
            start = self.start
        super().grid_forget()
        super().__init__(self.parent, from_=from_, to=to, orient='horizontal',
                       resolution=self.resolution)
        super().set(start)
        super().grid(row=self.row, column=self.column, rowspan=self.rowspan, 
                  columnspan=self.columnspan)
        
        
class NormalLabel(tk.Label):
    def __init__(self, parent, text=None, font=None, bd=None, relief=None,
                 row=None, column=None, rowspan=1, columnspan=1):
        super().__init__(parent, text=text, font=font, bd=bd, relief=relief)
        super().grid(row=row, column=column, rowspan=rowspan, 
                  columnspan=columnspan)
        self.parent = parent
        self.text = text
        self.font = font
        self.bd = bd
        self.relief = relief
        self.row = row
        self.column = column
        self.rowspan = rowspan
        self.columnspan = columnspan
        
    def updater(self, text=None):
        super().grid_forget()
        super().__init__(self.parent, text=text, font=self.font, 
                         bd=self.bd, relief=self.relief)    
        super().grid(row=self.row, column=self.column, rowspan=self.rowspan, 
                  columnspan=self.columnspan)
        
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
        frame.tkraise()
        frame.update_init(filenames)


# %% Quit

def quit_program():
    global gui
    gui.destroy()
    exit(0)

# %% Loading page

class LoadPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        button1 = BigButton(self, text="LOAD", height=int(GUI_HEIGHT / 4),
                           width=int(GUI_WIDTH / 4),  # style= 'my.TButton',
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

class FittingPage(tk.Frame):

    def update_init(self, filenames):

        self.filenames = filenames
        self.roi_locations = {}
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
        self.canvas.get_tk_widget().grid(row=0, column=10, columnspan=3, rowspan=14, sticky='E')
        
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

        min_int = self.min_int_slider.get()
        max_int = self.max_int_slider.get()
        min_sigma = self.min_sigma_slider.get()
        max_sigma = self.max_sigma_slider.get()
        wavelength = self.wavelength_input.get()
        if wavelength == "wavelength in nm":
            tk.messagebox.showerror("ERROR", "Please give laser wavelength")
            return
        else:
            wavelength = int(wavelength)

        self.roi_finder.change_settings(intensity_min=min_int,
                                        intensity_max=max_int,
                                        sigma_min=min_sigma, sigma_max=max_sigma)

        fitter = fitting.scipy_last_fit_guess(self.metadata, ROI_SIZE,
                                              wavelength, THRESHOLD,
                                              "ScipyLastFitGuess", 5)

        self.temp_roi_locations = self.roi_finder.main(self.frames[0], fitter)

        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)

        fig_sub.imshow(self.frames[0], extent=[0, self.frames[0].shape[1], self.frames[0].shape[0], 0],
                       aspect='auto')

        roi_locations_temp = self.temp_roi_locations - self.roi_finder.roi_size_1d
        roi_size = self.roi_finder.roi_size

        for roi in roi_locations_temp:
            rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                     linewidth=0.5, edgecolor='r', facecolor='none')
            fig_sub.add_patch(rect)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=10, columnspan=3, rowspan=14, sticky='E')

    # %% Fitting page, save ROI setings

    def save_roi_settings(self):

        wavelength = self.wavelength_input.get()
        if wavelength == "wavelength in nm":
            tk.messagebox.showerror("ERROR", "Please give laser wavelength")
            return
        else:
            wavelength = int(wavelength)

        self.fit_rois()

        self.roi_locations[self.dataset_index] = self.temp_roi_locations

        min_int = self.min_int_slider.get()
        max_int = self.max_int_slider.get()
        min_sigma = self.min_sigma_slider.get()
        max_sigma = self.max_sigma_slider.get()

        settings = {'max_int': max_int, 'min_int': min_int, 
                    'min_sigma': min_sigma, 'max_sigma': max_sigma,
                    'wavelength': wavelength}

        self.saved_settings[self.dataset_index] = settings

        self.roi_status.updater(text=str(len(self.roi_locations)) + 
                                " of " + str(len(self.filenames)) + 
                                " have settings")
        self.button_restore_saved.updater(command=lambda: self.restore_saved())

        if len(self.roi_locations) == len(self.filenames):
            self.button_fit.grid_forget()
            self.button_fit = BigButton(self, text="FIT", height=int(GUI_HEIGHT / 8),
                                       width=int(GUI_WIDTH / 8),
                                       command=lambda: self.start_fitting())
            self.button_fit.grid(row=24, column=10, columnspan=2, rowspan=5)

    # %% Fitting page, restore default settings

    def restore_default(self):

        self.roi_finder = roi_finding.roi_finder(9, self.frames[0])
        
        self.min_int_slider.updater(from_=0, to=self.roi_finder.intensity_max / 4, 
                                    start=self.roi_finder.intensity_min)
        self.max_int_slider.updater(from_=0, to=self.roi_finder.intensity_max,
                                    start=self.roi_finder.intensity_max)
        self.min_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=self.roi_finder.sigma_min)
        self.max_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=self.roi_finder.sigma_max)

    # %% Fitting page, switch between datasets

    def next_dataset(self):

        self.dataset_index += 1
        self.dataset_roi_status.updater(text="Dataset " + 
                                        str(self.dataset_index + 1) + " of " 
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
        self.canvas.get_tk_widget().grid(row=0, column=10, columnspan=3, rowspan=14, sticky='E')
        
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

            wavelength = self.saved_settings[self.dataset_index]['wavelength']
            self.wavelength_input.updater(text=str(wavelength))
        else:
            self.button_restore_saved.updater(state='disabled')
            
            self.wavelength_input = EntryPlaceholder(self, "wavelength in nm")
            self.wavelength_input.grid(row=12, column=6, columnspan=3)

    def previous_dataset(self):

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
        self.canvas.get_tk_widget().grid(row=0, column=10, columnspan=3, rowspan=14, sticky='E')
        
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

            wavelength = self.saved_settings[self.dataset_index]['wavelength']
            self.wavelength_input.updater(text=str(wavelength))
        else:
            self.button_restore_saved.updater(state='disabled')
            
            self.wavelength_input = EntryPlaceholder(self, "wavelength in nm")
            self.wavelength_input.grid(row=12, column=6, columnspan=3)

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
        roi_size = int(self.roi_var.get()[0])

        if n_processes > 1 and (method == "Phasor with intensity" or method == "Phasor without intensity"):
            cores_check = tk.messagebox.askokcancel("Just a heads up",
                                                    """Phasor will be used with one core since the 
                                                    overhead only slowes it down""")
            if not cores_check:
                return
            n_processes = 1

        self.start_time = time.time()

        for dataset_index, filename in enumerate(self.filenames):
            self.nd2.close()
            self.nd2 = ND2_Reader(filenames[self.dataset_index])
            self.frames = self.nd2
            self.metadata = self.nd2.metadata
            self.dataset_index = dataset_index

            roi_locations = self.roi_locations[dataset_index]

            wavelength = self.saved_settings[dataset_index]['wavelength']
            min_intensity = self.saved_settings[dataset_index]['min_int']

            if method == "Phasor with intensity":
                fitter = fitting.phasor_only_ROI_loop(self.metadata, roi_size, wavelength, min_intensity, method)
            elif method == "Phasor without intensity":
                fitter = fitting.phasor_only_ROI_loop_dumb(self.metadata, roi_size, wavelength, min_intensity, method)
            elif method == "Gaussian without background":
                fitter = fitting.scipy_last_fit_guess(self.metadata, roi_size, wavelength, min_intensity, method, 5)
            elif method == "Gaussian with background":
                fitter = fitting.scipy_last_fit_guess_background(self.metadata, roi_size, wavelength, min_intensity,
                                                                 method, 6)

            if start_frame == "Leave empty for start" and end_frame == "Leave empty for end":
                to_fit = self.frames
                end_frame = self.metadata['sequence_count']
                start_frame = 0
            elif start_frame == "Leave empty for start" and end_frame != "Leave empty for end":
                to_fit = self.frames[:int(end_frame)]
                start_frame = 0
            elif start_frame != "Leave empty for start" and end_frame == "Leave empty for end":
                to_fit = self.frames[int(start_frame):]
                end_frame = self.metadata['sequence_count']
            else:  # start_frame != "Leave empty for start" and end_frame != "Leave empty for end":
                to_fit = self.frames[int(start_frame):int(end_frame)]

            if n_processes > 1:
                shared = [None] * n_processes
                results = np.zeros((len(roi_locations) * len(self.frames), 9))

                frames_split = np.array_split(list(range(start_frame, end_frame)), n_processes)
                processes = [None] * n_processes

                for i in range(0, n_processes):
                    shared[i] = mp.Array('d', int(9 * len(roi_locations) * len(frames_split[i])))
                    processes[i] = (
                        mp.Process(target=mt_main, args=(filename, fitter, frames_split[i], roi_locations, shared[i])))

                for p in processes:
                    p.start()

                for p in processes:
                    p.join()

                counter = 0
                for i, share in enumerate(shared):
                    arr = np.frombuffer(share.get_obj())  # mp_arr and arr share the same memory
                    result = arr.reshape((len(roi_locations) * len(frames_split[i]), 9))
                    results[counter:counter + len(roi_locations) * len(frames_split[i])] = result
                    counter += len(roi_locations) * len(frames_split[i])

                results = results[results[:, 3] != 0]
            else:
                results = fitter.main(to_fit, self.metadata, roi_locations,
                                      gui=self)

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

        tk.messagebox.showinfo("Done!", "Done!")

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

        current_time = time.time()
        time_taken = current_time - self.start_time
        time_done_estimate = time_taken * 100 / progress + self.start_time
        tr = time.localtime(time_done_estimate)

        progress_text = str(round(progress, 1)) + "% done"
        time_text = str(tr[3]) + ":" + str(tr[4]) + ":" + str(tr[5]) + " " 
        + str(tr[2]) + "/" + str(tr[1])
        
        self.progress_status_label.updater(text=progress_text)
        self.time_status_label.udpater(text=time_text)

        self.update()

    # %% Fitting page, return to load page

    def load_new(self, controller):

        controller.show_load_frame(LoadPage)

    # %% Fitting page, restore saved settings

    def restore_saved(self):

        settings = self.saved_settings[self.dataset_index]
        
        self.min_int_slider.updater(from_=0, to=self.roi_finder.intensity_max / 4, 
                                    start=settings['min_int'])
        self.max_int_slider.updater(from_=0, to=self.roi_finder.intensity_max,
                                    start=settings['max_int'])
        self.min_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=settings['min_sigma'])
        self.max_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                      start=settings['max_sigma'])

        wavelength = settings['wavelength']
        self.wavelength_input.updater(text=str(wavelength))

    # %% Fitting page, initial declaration

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)
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
        
        min_int_label = tk.Label(self, text="Minimum Intensity", font=LARGE_FONT)
        min_int_label.grid(row=0, column=0, columnspan=3)
        
        self.min_int_slider = NormalSlider(self, from_=0, to=1000, 
                                           row=1, column=0, columnspan=3)

        max_int_label = tk.Label(self, text="Maximum Intensity", font=LARGE_FONT)
        max_int_label.grid(row=2, column=0, columnspan=3)
        
        self.max_int_slider = NormalSlider(self, from_=0, to=5000, 
                                           row=3, column=0, columnspan=3)

        min_sigma_label = tk.Label(self, text="Minimum Sigma", font=LARGE_FONT)
        min_sigma_label.grid(row=0, column=3, columnspan=3)
        
        self.min_sigma_slider = NormalSlider(self, from_=0, to=5, resolution=0.01,
                                             row=1, column=3, columnspan=3)

        max_sigma_label = tk.Label(self, text="Maximum Sigma", font=LARGE_FONT)
        max_sigma_label.grid(row=2, column=3, columnspan=3)
        
        self.max_sigma_slider = NormalSlider(self, from_=0, to=10, resolution=0.01,
                                             row=3, column=3, columnspan=3)
        
        self.button_left = NormalButton(self, text="<<", row=9, column=0)
        self.dataset_roi_status = NormalLabel(self, text= "TBD",
                                              row=9, column=1, columnspan=4)
        self.button_right = NormalButton(self, ">>", row=9, column=5)

        button_fit = ttk.Button(self, text="Fit", command=lambda: self.fit_rois())
        button_fit.grid(row=12, column=0, columnspan=1)

        button_restore = ttk.Button(self, text="Restore default",
                                    command=lambda: self.restore_default())
        button_restore.grid(row=12, column=1, columnspan=2)
        
        self.button_restore_saved = NormalButton(self, text="Restore saved",
                                                 state='disabled',
                                                 command=lambda: self.restore_saved(),
                                                 row=12, column=3, 
                                                  columnspan=2)

        button_save = ttk.Button(self, text="Save",
                                 command=lambda: self.save_roi_settings())
        button_save.grid(row=12, column=5, columnspan=1)

        self.fig = Figure(figsize=(GUI_HEIGHT / DPI * 0.7, GUI_HEIGHT / DPI * 0.7), dpi=DPI)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=10, columnspan=3, rowspan=14, sticky='E')

        # toolbar = NavigationToolbar2Tk(canvas, self)
        # toolbar.update()
        # canvas._tkcanvas.grid(row = 15, column = 10)

        wavelength_label = tk.Label(self, text="Laser Wavelength", font=LARGE_FONT)
        wavelength_label.grid(row=10, column=6, columnspan=3)

        self.wavelength_input = EntryPlaceholder(self, "wavelength in nm")
        self.wavelength_input.grid(row=12, column=6, columnspan=3)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(column=0, row=18, rowspan=2, columnspan=10, sticky='we')
        
        self.roi_status = NormalLabel(self, text="TBD", 
                                      row=21, column=0, columnspan=6)

        method_label = tk.Label(self, text="Method", font=LARGE_FONT)
        method_label.grid(column=0, row=23, columnspan=3)

        self.method_var = tk.StringVar(self)
        self.method_var.set(fit_options[2])

        method_drop = tk.OptionMenu(self, self.method_var, *fit_options)
        method_drop.grid(column=0, row=24, columnspan=3)

        roi_size_label = tk.Label(self, text="ROI size", font=LARGE_FONT)
        roi_size_label.grid(column=3, row=23, columnspan=3)

        self.roi_var = tk.StringVar(self)
        self.roi_var.set(roi_size_options[0])

        roi_drop = tk.OptionMenu(self, self.roi_var, *roi_size_options)
        roi_drop.grid(column=3, row=24, columnspan=3)

        frame_begin_label = tk.Label(self, text="Begin frame", font=LARGE_FONT)
        frame_begin_label.grid(row=27, column=0, columnspan=3)

        self.frame_begin_input = EntryPlaceholder(self, "Leave empty for start")
        self.frame_begin_input.grid(row=28, column=0, columnspan=3)

        frame_end_label = tk.Label(self, text="End frame", font=LARGE_FONT)
        frame_end_label.grid(row=27, column=3, columnspan=3)

        self.frame_end_input = EntryPlaceholder(self, "Leave empty for end")
        self.frame_end_input.grid(row=28, column=3, columnspan=3)

        cores_label = tk.Label(self, text="#cores", font=LARGE_FONT)
        cores_label.grid(row=23, column=6)

        self.total_cores = mp.cpu_count()
        cores_options = [1, int(self.total_cores / 2), int(self.total_cores * 3 / 4), int(self.total_cores)]
        self.cores_var = tk.IntVar(self)
        self.cores_var.set(cores_options[0])
        self.cores_drop = tk.OptionMenu(self, self.cores_var, *cores_options)
        self.cores_drop.grid(row=24, column=6)

        dimensions_label = tk.Label(self, text="pixels or nm", font=LARGE_FONT)
        dimensions_label.grid(row=27, column=6)

        dimension_options = ["nm", "pixels"]
        self.dimension = tk.StringVar(self)
        self.dimension.set(dimension_options[0])
        self.dimension_drop = tk.OptionMenu(self, self.dimension, *dimension_options)
        self.dimension_drop.grid(row=28, column=6)

        self.button_fit = BigButton(self, text="FIT", height=int(GUI_HEIGHT / 8),
                                   width=int(GUI_WIDTH / 8), state='disabled',
                                   command=lambda: self.start_fitting())  # , style= 'my.TButton')
        self.button_fit.grid(row=24, column=10, columnspan=2, rowspan=5)

        progress_label = tk.Label(self, text="Progress", font=LARGE_FONT)
        progress_label.grid(row=24, column=12)
        
        self.progress_status_label = NormalLabel(self, text="Not yet started",
                                              bd=1, relief='sunken',
                                              row=25, column=12)

        time_label = tk.Label(self, text="Estimated time done", font=LARGE_FONT)
        time_label.grid(row=27, column=12)
        
        self.time_status_label = NormalLabel(self, text="Not yet started",
                                          bd=1, relief='sunken',
                                          row=28, column=12)

        button_quit = ttk.Button(self, text="Quit", command=quit_program)
        button_quit.grid(row=50, column=8)

        button_load_new = ttk.Button(self, text="Load new", command=lambda: self.load_new(controller))
        button_load_new.grid(row=50, column=6)

        for i in range(0, 13):
            self.columnconfigure(i, weight=1)

        for i in range(0, 51):
            self.rowconfigure(i, weight=1)


# %% START GUI
if __name__ == '__main__':
    gui = MbxPython()
    gui.geometry(str(GUI_WIDTH) + "x" + str(GUI_HEIGHT) + "+" + str(GUI_WIDTH_START) + "+" + str(GUI_HEIGHT_START))
    gui.mainloop()
    
    ttk_style = ttk.Style(gui)
    ttk_style.configure('my.TButton', font=('Verdana', 1000))
    ttk_style.configure("Placeholder.TEntry", foreground="#d5d5d5")
