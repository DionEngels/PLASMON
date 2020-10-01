# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
MBx Python Data Analysis

main_GUI

This package is for the GUI of Mbx Python.

----------------------------

v0.1: first version GUI: 13/07/2020
v0.2: new ROI finding method: 20/07/2020
v0.3: Styling, ready for review on functional level
v0.4: Bug fix to batch loading, MATLAB-ready output, settings and results text file output.
v0.4.1: different directory output
v0.4.2: prevent overwriting output
v0.4.3: settings dict, one command to change dataset
v0.4.4: saved standard values
v0.5: removed pixel min, moved min_corr in GUI
v0.5.1: fixed size GUI
v0.6: matplotlib fix, PIMS 0.5, and own ND2Reader class to prevent warnings
v0.7: drift correction v1: 31/07/2020
v0.7.1: add button to restore to different saved settings
v0.7.2: save figures and save drift
v0.7.3: new metadata
v0.8: GUI part of HSM
v0.8.1: metadata v3
v1.0: roll-out version one
v1.0.1: instant bugfix open all .nd2s
v1.1: Bugfixes and improved figures (WIP)
v1.1.1: tiny bugfixes: 10/08/2020
v1.2: GUI and output improvement based on Sjoerd's feedback, HSM: 27/08/2020 - 13/09/2020
v1.2.0.1: small bugfix
v1.2.1: HSM checks, bugfixes
v1.2.2: load from other bugfixes
v1.3: HSM to eV: 24/09/2020
v1.4: HSM output back to nm, while fitting in eV: 29/09/2020
"""
__version__ = "1.4"
__self_made__ = True

# GENERAL IMPORTS
from os import getcwd, mkdir, environ, listdir  # to get standard usage
from tempfile import mkdtemp
import sys
import time  # for timekeeping
from win32api import GetSystemMetrics  # Get sys info

environ['MPLCONFIGDIR'] = mkdtemp()

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.patches as patches
from scipy.io import loadmat
from scipy.ndimage import median_filter

# GUI
import tkinter as tk  # for GUI
from tkinter import ttk  # GUI styling
from tkinter.filedialog import askopenfilenames, askdirectory  # for popup that asks to select .nd2's or folders

# Own code
import src.roi_finding as roi_finding
import src.tt as fitting
import src.tools as tools
import src.drift_correction as drift_correction
import src.figure_making as figuring
import src.output as outputting
import src.nd2_reading as nd2_reading
from src.hsm import normxcorr2, normxcorr2_large
import src.hsm as hsm

# Multiprocessing
import multiprocessing as mp

mpl.use("TkAgg")  # set back end to TK

# %% Initializations. Defining filetypes, fonts, paddings, input sizes, and GUI sizes.

FILETYPES = [("ND2", ".nd2")]
FILETYPES_LOAD_FROM_OTHER = [(".npy and .mat", ".npy"), (".npy and .mat", ".mat")]

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
GUI_WIDTH = 1344  # int(width * 0.70)
GUI_HEIGHT = 756  # int(height * 0.70)
GUI_WIDTH_START = int((width - GUI_WIDTH) / 2)
GUI_HEIGHT_START = int((height - GUI_HEIGHT) / 2)
DPI = 100

# %% Options for dropdown menus

fit_options = ["Gaussian - Fit bg", "Gaussian - Estimate bg",
               "Phasor + Intensity", "Phasor + Sum", "Phasor"]

rejection_options = ["Strict", "Loose", "None"]

roi_size_options = ["7x7", "9x9"]


# %% Multiprocessing main


def mt_main(name, fitter, frames_split, roi_locations, shared, q):
    """
    Main for multiprocessing. Loads .nd2, takes certain frames of this and sends this to fitter.
    Communicates using shared memory (for results) and queue for status updates

    Parameters
    ----------
    name : name of ND2 file
    fitter : Fitter to be used
    frames_split : What frames this instance has to fit of the total .nd2
    roi_locations : ROI locations
    shared : Shared memory. The place where the results will be placed
    q : Queue. Used to send status updates to GUI.

    Returns
    -------
    Fills shared memory

    """
    nd2 = nd2_reading.ND2ReaderSelf(name)
    frames = nd2
    metadata = nd2.get_metadata()
    metadata['num_frames'] = len(frames_split)
    frames = frames[frames_split]

    local_result = fitter.main(frames, metadata, roi_locations,
                               q=q, start_frame=frames_split[0])

    for result_index, result in enumerate(local_result):
        shared[9 * result_index:9 * (result_index + 1)] = result[:]


# %% Divert errors


def show_error_critical(self, exc, val, tb):
    show_error(True)


def show_error(critical):
    exc_type, exc_value, exc_traceback = sys.exc_info()  # most recent (if any) by default
    while True:
        try:
            self_made = exc_traceback.tb_next.tb_frame.f_globals['__self_made__']
        except:
            self_made = None
        if self_made is None:
            break
        elif self_made:
            exc_traceback = exc_traceback.tb_next
        else:
            break
    traceback_details = {
        'filename': exc_traceback.tb_frame.f_code.co_filename,
        'lineno': exc_traceback.tb_lineno,
        'name': exc_traceback.tb_frame.f_code.co_name,
        'type': exc_type.__name__,
        'message': exc_value
    }
    if critical:
        tk.messagebox.showerror("Critical error. Send screenshot to Dion. PROGRAM WILL STOP",
                                message=str(traceback_details))
    else:
        tk.messagebox.showerror("Error. Send screenshot to Dion. PROGRAM WILL CONTINUE", message=str(traceback_details))


# %% Close GUI

def quit_gui(gui):
    gui.quit()
    sys.exit(0)

# %% Own buttons / fields


class BigButton(ttk.Frame):
    """
    Big button, used for FIT and LOAD
    """

    def __init__(self, parent, height=None, width=None, text="", command=None, state='enabled'):
        ttk.Frame.__init__(self, parent, height=height, width=width)

        self.pack_propagate(0)
        self._btn = ttk.Button(self, text=text, command=command, state=state)
        self._btn.pack(fill=tk.BOTH, expand=1)
        self._btn["style"] = "Big.TButton"

    def updater(self, state='enabled'):
        self._btn['state'] = state


class FigureFrame(tk.Frame):
    """
    Frame in which Figure sits.
    """

    def __init__(self, parent, height=None, width=None, dpi=DPI):
        tk.Frame.__init__(self, parent, height=height + 40, width=width,
                          highlightbackground="black", highlightthickness=2)

        self.pack_propagate(0)
        self.fig = Figure(figsize=(height / dpi, width / dpi), dpi=dpi)

        self.parent = parent
        self.dpi = dpi
        self.width = width
        self.height = height

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.toolbar.configure(background="White")

    def updater(self, frame, roi_locations=None, roi_size=None):
        """
        Updater. Takes existing frame with figure and places new figure in it

        Parameters
        ----------
        frame : New frame to be shown
        roi_locations : optional, possible ROI locations to be highlighted. The default is None.
        roi_size : optional, ROI size in case ROIs are highlighted. The default is None.

        Returns
        -------
        Updated figure.

        """
        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(frame, extent=[0, frame.shape[1], frame.shape[0], 0], aspect='auto')

        if roi_locations is not None and roi_size is not None:
            roi_locations_temp = roi_locations - roi_size

            for roi_index, roi in enumerate(roi_locations_temp):
                rect = patches.Rectangle((roi[1], roi[0]), roi_size * 2 + 1, roi_size * 2 + 1,
                                         linewidth=0.5, edgecolor='r', facecolor='none')
                fig_sub.add_patch(rect)
                fig_sub.text(roi[1], roi[0], str(roi_index+1), color='red', fontsize='small')

        self.canvas.draw()
        self.toolbar.update()


class EntryPlaceholder(ttk.Entry):
    """
    Entry with a placeholder text in grey
    """

    def __init__(self, master=None, placeholder="PLACEHOLDER", *args, **kwargs):
        super().__init__(master, *args, style="Placeholder.TEntry", **kwargs)
        self.placeholder = placeholder

        self.insert("0", self.placeholder)
        self.bind("<FocusIn>", self._clear_placeholder)
        self.bind("<FocusOut>", self._add_placeholder)

    def _clear_placeholder(self, e):
        self.delete("0", "end")
        if self["style"] == "Placeholder.TEntry":
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
    """
    My normal button, again with an updater function to update the button.
    Only buttons that need updating use this class
    """

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
    """
    My normal slider, again with an updater function to update the slider.
    """

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
    """
    My normal label, again with an updater function to update the label.
    """

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


# %% Controller


class MbxPython(tk.Tk):
    """
    Controller of GUI. This container calls the page we need (load page or fitting page)
    """

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

    def show_fitting_frame(self, cont, filenames):
        frame = self.frames[cont]
        frame.update_init(self, cont, filenames)
        frame.tkraise()


# %% Loading page


class LoadPage(tk.Frame):
    """
    Loading page. On this page, there is only a big button to load ND2s
    """

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
                                     initialdir=getcwd())

        if len(filenames) == 0:
            return

        controller.show_fitting_frame(FittingPage, filenames)


# %% Fitting page, initial update (after loading)


class FittingPage(tk.Frame):
    """
    Fitting page. This is main page of the program
    """

    def __init__(self, parent, container, reset=False):
        """
        Initial declaration of fitting page. Buttons/sliders/etc. that are to be updated are linked to self using
        self-made classes.

        Parameters
        ----------
        parent : Parent page, the controller. MBxPython.
        container : The dict containing all the pages
        reset : optional, if reset it TRUE, this means that update_init is calling this and page will only be reset.
        The default is False.

        Returns
        -------
        None. GUI.

        """

        if not reset:
            tk.Frame.__init__(self, parent)
            self.configure(bg='white')
        self.nd2 = None
        self.frames = None
        self.roi_finder = None
        self.metadata = None
        self.temp_roi_locations = None
        self.filenames = None
        self.filename_short = None
        self.start_time = None
        self.roi_locations = {}
        self.dataset_index = 0
        self.saved_settings = {}
        self.default_settings = {}
        self.histogram = None
        self.to_hist = None
        self.roi_fitter = None
        self.saved_settings_list = []
        self.hsm_folder_full = None
        self.hsm_folder_show = None

        roi_finding_label = tk.Label(self, text="ROI finding", font=FONT_HEADER, bg='white')
        roi_finding_label.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_int_label = tk.Label(self, text="Minimum Intensity", font=FONT_LABEL, bg='white')
        min_int_label.grid(row=1, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_int_slider = NormalSlider(self, from_=0, to=1000,
                                           row=2, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_int_histogram = ttk.Button(self, text="Graph",
                                       command=lambda: self.fun_histogram("min_int"))
        min_int_histogram.grid(row=2, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_int_histogram_select = ttk.Button(self, text="Select min",
                                              command=lambda: self.histogram_select("min_int"))
        min_int_histogram_select.grid(row=2, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_int_label = tk.Label(self, text="Maximum Intensity", font=FONT_LABEL, bg='white')
        max_int_label.grid(row=1, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.max_int_slider = NormalSlider(self, from_=0, to=5000,
                                           row=2, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_int_histogram_select = ttk.Button(self, text="Select max",
                                              command=lambda: self.histogram_select("max_int"))
        max_int_histogram_select.grid(row=2, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_sigma_label = tk.Label(self, text="Minimum Sigma", font=FONT_LABEL, bg='white')
        min_sigma_label.grid(row=3, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_sigma_slider = NormalSlider(self, from_=0, to=5, resolution=0.01,
                                             row=4, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_sigma_histogram = ttk.Button(self, text="Graph",
                                         command=lambda: self.fun_histogram("min_sigma"))
        min_sigma_histogram.grid(row=4, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_sigma_histogram_select = ttk.Button(self, text="Select min",
                                                command=lambda: self.histogram_select("min_sigma"))
        min_sigma_histogram_select.grid(row=4, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_sigma_label = tk.Label(self, text="Maximum Sigma", font=FONT_LABEL, bg='white')
        max_sigma_label.grid(row=3, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.max_sigma_slider = NormalSlider(self, from_=0, to=10, resolution=0.01,
                                             row=4, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        max_sigma_histogram_select = ttk.Button(self, text="Select max",
                                                command=lambda: self.histogram_select("max_sigma"))
        max_sigma_histogram_select.grid(row=4, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=5, column=0, rowspan=1, columnspan=40, sticky='we')

        advanced_label = tk.Label(self, text="Advanced settings", font=FONT_SUBHEADER, bg='white')
        advanced_label.grid(row=6, column=0, columnspan=40, sticky='EW', padx=PAD_SMALL)

        min_corr_label = tk.Label(self, text="Minimum Correlation", font=FONT_LABEL, bg='white')
        min_corr_label.grid(row=6, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.min_corr_slider = NormalSlider(self, from_=0, to=1, resolution=0.005,
                                            row=7, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_corr_histogram = ttk.Button(self, text="Graph",
                                        command=lambda: self.fun_histogram("corr_min"))
        min_corr_histogram.grid(row=7, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        min_corr_histogram_select = ttk.Button(self, text="Graph select",
                                               command=lambda: self.histogram_select("corr_min"))
        min_corr_histogram_select.grid(row=7, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        roi_size_label = tk.Label(self, text="ROI size", bg='white', font=FONT_LABEL)
        roi_size_label.grid(row=9, column=0, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.roi_var = tk.StringVar(self)
        roi_drop = ttk.OptionMenu(self, self.roi_var, roi_size_options[0], *roi_size_options)
        roi_drop.grid(row=9, column=5, columnspan=5, sticky='EW', padx=PAD_SMALL)

        filter_label = tk.Label(self, text="Filter size", bg='white', font=FONT_LABEL)
        filter_label.grid(row=9, column=10, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.filter_input = EntryPlaceholder(self, "9", width=INPUT_SMALL)
        self.filter_input.grid(row=9, column=15, columnspan=5)

        roi_side_label = tk.Label(self, text="Side spacing", bg='white', font=FONT_LABEL)
        roi_side_label.grid(row=9, column=20, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.roi_side_input = EntryPlaceholder(self, "11", width=INPUT_SMALL)
        self.roi_side_input.grid(row=9, column=25, columnspan=5)

        inter_roi_label = tk.Label(self, text="ROI spacing", bg='white', font=FONT_LABEL)
        inter_roi_label.grid(row=9, column=30, columnspan=5, sticky='EW', padx=PAD_SMALL)

        self.inter_roi_input = EntryPlaceholder(self, "6", width=INPUT_SMALL)
        self.inter_roi_input.grid(row=9, column=35, columnspan=5)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=10, column=0, rowspan=1, columnspan=40, sticky='we')

        hsm_label = tk.Label(self, text="HSM", font=FONT_HEADER, bg='white')
        hsm_label.grid(row=11, column=0, columnspan=40, sticky='EW', padx=PAD_SMALL)

        self.hsm_load = NormalButton(self, text="Load HSM data", row=12, column=0, columnspan=8, sticky='EW',
                                     padx=PAD_SMALL, command=lambda: self.load_hsm())

        self.hsm_clear = NormalButton(self, text="Clear HSM", row=13, column=0, columnspan=8, sticky='EW',
                                      padx=PAD_SMALL, command=lambda: self.clear_hsm())

        hsm_disp_label = tk.Label(self, text="Loaded folder:", font=FONT_LABEL, bg='white', anchor='e')
        hsm_disp_label.grid(row=12, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.hsm_folder_disp = tk.Label(self, text="", font=FONT_LABEL, bg='white')
        self.hsm_folder_disp.grid(row=12, column=16, columnspan=24, sticky='EW', padx=PAD_SMALL)

        hsm_correct_label = tk.Label(self, text="Correction file:", font=FONT_LABEL, bg='white', anchor='e')
        hsm_correct_label.grid(row=13, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        path_hsm_corrections = getcwd() + "/spectral_corrections"

        self.hsm_correct_options = listdir(path_hsm_corrections)
        self.hsm_correct_options = [option[:-4] for option in self.hsm_correct_options]  # remove .mat in name

        self.hsm_correct_var = tk.StringVar(self)
        self.hsm_correct_drop = ttk.OptionMenu(self, self.hsm_correct_var, [], *self.hsm_correct_options)
        self.hsm_correct_drop.grid(row=13, column=16, columnspan=24, sticky="ew")

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=15, column=0, rowspan=1, columnspan=40, sticky='we')

        self.button_left = NormalButton(self, text="<<", row=17, column=8, columnspan=5, sticky='EW',
                                        padx=PAD_SMALL)
        self.dataset_roi_status = NormalLabel(self, text="TBD",
                                              row=17, column=13, columnspan=22, font=FONT_LABEL, sticky='EW')
        self.button_right = NormalButton(self, ">>", row=17, column=35, columnspan=5, sticky='EW',
                                         padx=PAD_SMALL)

        button_find_rois = ttk.Button(self, text="Find ROIs", command=lambda: self.fit_rois())
        button_find_rois.grid(row=16, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.number_of_rois = NormalLabel(self, text="TBD", row=17, column=0, columnspan=8, font=FONT_LABEL,
                                          padx=PAD_SMALL, sticky='EW')

        button_restore = ttk.Button(self, text="Restore default",
                                    command=lambda: self.restore_default())
        button_restore.grid(row=16, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.button_load_other_video = NormalButton(self, text="Load from other",
                                                    command=lambda: self.load_from_other(),
                                                    row=16, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.button_restore_saved = NormalButton(self, text="Restore saved",
                                                 state='disabled',
                                                 command=lambda: self.restore_saved(),
                                                 row=16, column=24,
                                                 columnspan=8, sticky='EW', padx=PAD_SMALL)

        button_save = ttk.Button(self, text="Save",
                                 command=lambda: self.save_roi_settings())
        button_save.grid(row=16, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.fig = FigureFrame(self, height=GUI_WIDTH * 0.4, width=GUI_WIDTH * 0.4, dpi=DPI)
        self.fig.grid(row=0, column=40, columnspan=10, rowspan=18, sticky='EW', padx=PAD_SMALL)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=18, column=0, rowspan=2, columnspan=50, sticky='we')

        fit_area_label = tk.Label(self, text="Fitting", font=FONT_HEADER, bg='white')
        fit_area_label.grid(row=21, column=0, columnspan=10, sticky='EW', padx=PAD_SMALL)

        self.roi_status = NormalLabel(self, text="TBD",
                                      row=21, column=10, columnspan=30, font=FONT_STATUS)

        method_label = tk.Label(self, text="Method", font=FONT_LABEL, bg='white')
        method_label.grid(row=23, column=0, columnspan=12, sticky='EW', padx=PAD_SMALL)

        self.method_var = tk.StringVar(self)
        method_drop = ttk.OptionMenu(self, self.method_var, fit_options[1], *fit_options)
        method_drop.grid(row=24, column=0, columnspan=12, sticky="ew")

        rejection_label = tk.Label(self, text="Rejection", bg='white', font=FONT_LABEL)
        rejection_label.grid(row=23, column=12, columnspan=7, sticky='EW', padx=PAD_SMALL)

        self.rejection_var = tk.StringVar(self)
        rejection_drop = ttk.OptionMenu(self, self.rejection_var, rejection_options[1], *rejection_options)
        rejection_drop.grid(row=24, column=12, columnspan=7, sticky='EW', padx=PAD_SMALL)

        cores_label = tk.Label(self, text="#cores", font=FONT_LABEL, bg='white')
        cores_label.grid(row=23, column=19, columnspan=7, sticky='EW', padx=PAD_BIG)

        self.total_cores = mp.cpu_count()
        cores_options = [1, int(self.total_cores / 2), int(self.total_cores * 3 / 4), int(self.total_cores)]
        self.cores_var = tk.IntVar(self)
        self.cores_drop = ttk.OptionMenu(self, self.cores_var, cores_options[0], *cores_options)
        self.cores_drop.grid(row=24, column=19, columnspan=7, sticky='EW', padx=PAD_BIG)

        dimensions_label = tk.Label(self, text="pixels or nm", font=FONT_LABEL, bg='white')
        dimensions_label.grid(row=23, column=26, columnspan=7, sticky='EW', padx=PAD_BIG)

        dimension_options = ["nm", "pixels"]
        self.dimension = tk.StringVar(self)
        self.dimension_drop = ttk.OptionMenu(self, self.dimension, dimension_options[0], *dimension_options)
        self.dimension_drop.grid(row=24, column=26, columnspan=7, sticky='EW', padx=PAD_BIG)

        figures_label = tk.Label(self, text="All Figures", font=FONT_LABEL, bg='white')
        figures_label.grid(row=23, column=33, columnspan=7, sticky='EW', padx=PAD_BIG)

        self.figures_var = tk.StringVar(self, value="Overview")
        self.figures_check = tk.Checkbutton(self, variable=self.figures_var, onvalue="All", offvalue="Overview",
                                            bg="white")
        self.figures_check.grid(row=24, column=33, columnspan=7, sticky='EW', padx=PAD_BIG)

        frame_begin_label = tk.Label(self, text="Begin frame", font=FONT_LABEL, bg='white')
        frame_begin_label.grid(row=27, column=0, columnspan=16, sticky='EW', padx=PAD_BIG)

        self.frame_begin_input = EntryPlaceholder(self, "Leave empty for start", width=INPUT_BIG)
        self.frame_begin_input.grid(row=28, column=0, columnspan=16)

        frame_end_label = tk.Label(self, text="End frame", font=FONT_LABEL, bg='white')
        frame_end_label.grid(row=27, column=16, columnspan=16, sticky='EW', padx=PAD_BIG)

        self.frame_end_input = EntryPlaceholder(self, "Leave empty for end", width=INPUT_BIG)
        self.frame_end_input.grid(row=28, column=16, columnspan=16)

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

        version_label = tk.Label(self, text="MBx Python, version: " + __version__, font=FONT_LABEL, bg='white',
                                 anchor='w')
        version_label.grid(row=50, column=0, columnspan=20, sticky='EW', padx=PAD_SMALL)

        button_load_new = ttk.Button(self, text="Load new", command=lambda: self.load_new(parent))
        button_load_new.grid(row=50, column=45, columnspan=2, sticky='EW', padx=PAD_SMALL)

        button_quit = ttk.Button(self, text="Quit", command=lambda: quit_gui(self))
        button_quit.grid(row=50, column=48, columnspan=2, sticky='EW', padx=PAD_SMALL)

        for i in range(0, 49):
            self.columnconfigure(i, weight=1, minsize=18)

        for i in range(0, 51):
            self.rowconfigure(i, weight=1)

    # %% Update after loading of ND2s

    def update_init(self, parent, container, filenames):
        """
        Initial update. The init function is already called by the controller, so this is called after the ND2 are
        loaded in. Updates all kinds of things such as the sliders.

        Parameters
        ----------
        parent : This is controller page. MBxPython
        container : A dict in which all the pages reside. Only used for __init__
        filenames : Filenames of the loaded ND2s

        Returns
        -------
        Updated page

        """
        self.__init__(parent, container, reset=True)

        self.filenames = filenames
        self.dataset_index = 0

        self.filename_short = filenames[self.dataset_index].split("/")[-1][:-4]

        self.dataset_roi_status.updater(text=self.filename_short + " (" + str(self.dataset_index + 1) + " of " +
                                        str(len(filenames)) + ")")
        self.roi_status.updater(text="0 of " + str(len(filenames)) + " have settings")

        self.nd2 = nd2_reading.ND2ReaderSelf(filenames[self.dataset_index])
        self.frames = self.nd2
        self.metadata = self.nd2.get_metadata()

        self.fig.updater(self.frames[0])

        self.restore_default()

        if self.dataset_index == len(filenames) - 1:
            self.button_right.updater(state='disabled')
        else:
            self.button_right.updater(command=lambda: self.change_dataset(1))

        if self.dataset_index == 0:
            self.button_left.updater(state='disabled')
        else:
            self.button_left.updater(command=lambda: self.change_dataset(-1))

    # %% Fitting page, return to load page

    def load_new(self, controller):
        """
        Return to load page

        Parameters
        ----------
        controller : Calls controller to return to load page

        Returns
        -------
        None.

        """
        controller.show_load_frame(LoadPage)

    # %% Fitting page, restore default settings

    def restore_default(self):
        """
        Restores defaults. Called when calling the page for the first time, and when switching dataset

        Returns
        -------
        None.

        """

        if self.dataset_index in self.default_settings:
            defaults = self.default_settings[self.dataset_index]
            self.roi_finder = roi_finding.RoiFinder(self.frames[0], self.roi_fitter, settings=defaults)

            self.min_int_slider.updater(from_=0, to=self.roi_finder.int_max / 4,
                                        start=defaults['int_min'])
            self.max_int_slider.updater(from_=0, to=self.roi_finder.int_max,
                                        start=defaults['int_max'])
            self.min_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                          start=defaults['sigma_min'])
            self.max_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                          start=defaults['sigma_max'])
            self.min_corr_slider.updater(from_=0, to=1, start=defaults['corr_min'])
        else:
            self.roi_fitter = fitting.Gaussian(7, {}, "None", "Gaussian", 5, 300)

            self.roi_finder = roi_finding.RoiFinder(self.frames[0], self.roi_fitter)

            self.min_int_slider.updater(from_=0, to=self.roi_finder.int_max / 4,
                                        start=self.roi_finder.int_min)
            self.max_int_slider.updater(from_=0, to=self.roi_finder.int_max,
                                        start=self.roi_finder.int_max)
            self.min_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                          start=self.roi_finder.sigma_min)
            self.max_sigma_slider.updater(from_=0, to=self.roi_finder.sigma_max,
                                          start=self.roi_finder.sigma_max)
            self.min_corr_slider.updater(from_=0, to=1, start=self.roi_finder.corr_min)

        self.roi_var.set(roi_size_options[0])

        self.filter_input.updater()
        self.roi_side_input.updater()
        self.inter_roi_input.updater()

        self.update()

        self.fit_rois()

        settings = self.read_out_settings()

        settings['processed_frame'] = self.roi_finder.frame

        self.default_settings[self.dataset_index] = settings

    # %% Read out sliders and other settings, return dictionary

    def read_out_settings(self):
        """
        Function that reads out all the settings, and saves it to a dict which it returns

        Returns
        -------
        settings: a dictionary of all read-out settings
        """
        int_min = self.min_int_slider.get()
        int_max = self.max_int_slider.get()
        sigma_min = self.min_sigma_slider.get()
        sigma_max = self.max_sigma_slider.get()
        corr_min = self.min_corr_slider.get()
        roi_size = int(self.roi_var.get()[0])

        filter_size = int(self.filter_input.get())
        roi_side = int(self.roi_side_input.get())
        inter_roi = int(self.inter_roi_input.get())

        hsm_directory = self.hsm_folder_full
        hsm_correction = self.hsm_correct_var.get()

        settings = {'int_max': int_max, 'int_min': int_min,
                    'sigma_min': sigma_min, 'sigma_max': sigma_max,
                    'corr_min': corr_min, 'roi_size': roi_size, 'filter_size': filter_size,
                    'roi_side': roi_side, 'inter_roi': inter_roi,
                    'hsm_directory': hsm_directory, 'hsm_correction': hsm_correction}

        return settings

    # %% Fitting page, fit ROIs

    def fit_rois(self):
        """
        Function that takes all the inputs and uses it to fit ROIs

        Returns
        -------
        Updates the figures.
        Also returns boolean whether or not it was a success.

        """
        settings = self.read_out_settings()

        if settings['roi_side'] < int((settings['roi_size'] - 1) / 2):
            tk.messagebox.showerror("ERROR", "Distance to size cannot be smaller than 1D ROI size")
            return False
        if settings['filter_size'] % 2 != 1:
            tk.messagebox.showerror("ERROR", "Filter size should be odd")
            return False

        self.roi_fitter = fitting.Gaussian(settings['roi_size'], {}, "None", "Gaussian", 5, 300)

        self.roi_finder.change_settings(settings)

        self.temp_roi_locations = self.roi_finder.main(self.roi_fitter)

        self.fig.updater(self.frames[0], roi_locations=self.temp_roi_locations, roi_size=self.roi_finder.roi_size_1d)
        self.number_of_rois.updater(text=str(len(self.temp_roi_locations)) + " ROIs found")

        return True

    # %% Fitting page, save ROI settings

    def save_roi_settings(self):
        """
        Saves current settings of sliders etc. to dict

        Returns
        -------
        None.

        """
        success = self.fit_rois()

        if not success:
            return

        self.roi_locations[self.dataset_index] = self.temp_roi_locations.copy()

        max_its = self.roi_finder.find_snr(self.roi_fitter)

        settings = self.read_out_settings()
        settings['max_its'] = max_its

        self.saved_settings[self.dataset_index] = settings

        self.roi_status.updater(text=str(len(self.roi_locations)) + " of " + str(len(
            self.filenames)) + " have settings")
        self.button_restore_saved.updater(command=lambda: self.restore_saved())

        if len(self.roi_locations) == len(self.filenames):
            self.button_fit.updater(state='enabled')

    # %% Fitting page, restore saved settings

    def restore_saved(self):
        """
        Restores saved settings to sliders etc.

        Returns
        -------
        None, updates GUI

        """
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

        if settings['roi_size'] == 7:
            self.roi_var.set(roi_size_options[0])
        else:
            self.roi_var.set(roi_size_options[1])

        self.filter_input.updater(settings['filter_size'])
        self.roi_side_input.updater(settings['roi_side'])
        self.inter_roi_input.updater(settings['inter_roi'])

        self.roi_finder.change_settings(settings)
        self.roi_finder.frame = self.default_settings[self.dataset_index]['processed_frame']

        self.temp_roi_locations = self.roi_locations[self.dataset_index]

        self.fig.updater(self.frames[0], roi_locations=self.temp_roi_locations, roi_size=self.roi_finder.roi_size_1d)
        self.number_of_rois.updater(text=str(len(self.temp_roi_locations)) + " ROIs found")

        self.update()

    # %% Load from previously processed dataset

    def load_from_other(self):
        def load_in_files(filenames):
            frame_zero_old = np.load([file for file in filenames if file.endswith('.npy')][0])
            roi_locations = loadmat([file for file in filenames if file.endswith('.mat')][0])
            roi_locations = [roi_locations[key] for key in roi_locations.keys() if not key.startswith('_')][0]
            roi_locations = tools.roi_to_python_coordinates(roi_locations, frame_zero_old.shape[0])
            return roi_locations, frame_zero_old

        def correlate_frames(frame_old, frame_new):
            if frame_old.shape == frame_new.shape:
                corr = normxcorr2(frame_old, frame_new)
                maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
                offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
            else:
                corr = normxcorr2_large(frame_old, frame_new)
                maxima = np.transpose(np.asarray(np.where(corr == np.amax(corr))))[0]
                offset = maxima - np.asarray(frame_old.shape) + np.asarray([1, 1])
            return offset

        tk.Tk().withdraw()
        while True:
            filenames = askopenfilenames(filetypes=FILETYPES_LOAD_FROM_OTHER,
                                         title="Select both frame_zero and ROI locations that you want to use",
                                         initialdir=getcwd())
            if len(filenames) == 2:
                break
            else:
                check = tk.messagebox.askokcancel("ERROR",
                                                  "Select one ROI locations file and one frame_zero. Try again?")
                if not check:
                    return

        roi_locations, frame_zero_old = load_in_files(filenames)
        frame_zero = np.asarray(self.frames[0])

        background = median_filter(frame_zero_old, size=9)
        frame_zero_old = frame_zero_old.astype(np.int16) - background
        background = median_filter(frame_zero, size=9)
        frame_zero = frame_zero.astype(np.int16) - background

        offset = correlate_frames(frame_zero_old, frame_zero)

        try:
            if offset[0] < 50 and offset[1] < 50:
                roi_locations += offset
            else:  # other background mode as backup
                load_in_files(filenames)
                frame_zero = np.asarray(self.frames[0])
                background = median_filter(frame_zero_old, size=9, mode='constant', cval=np.mean(frame_zero_old))
                frame_zero_old = frame_zero_old.astype(np.int16) - background
                background = median_filter(frame_zero, size=9, mode='constant', cval=np.mean(frame_zero))
                frame_zero = frame_zero.astype(np.int16) - background
                offset = correlate_frames(frame_zero_old, frame_zero)
                roi_locations += offset
        except:
            tk.messagebox.askokcancel("ERROR", "You did not select a proper ROI locations file. Try again.")
            return

        self.temp_roi_locations = roi_locations
        self.roi_locations[self.dataset_index] = self.temp_roi_locations.copy()

        max_its = self.roi_finder.find_snr_load_from_other(self.roi_fitter, roi_locations)

        settings = self.read_out_settings()
        settings['max_its'] = max_its

        self.saved_settings[self.dataset_index] = settings

        self.roi_status.updater(text=str(len(self.roi_locations)) + " of " + str(len(
            self.filenames)) + " have settings")
        self.button_restore_saved.updater(command=lambda: self.restore_saved())

        if len(self.roi_locations) == len(self.filenames):
            self.button_fit.updater(state='enabled')

        self.fig.updater(self.frames[0], roi_locations=self.temp_roi_locations, roi_size=self.roi_finder.roi_size_1d)
        self.number_of_rois.updater(text=str(len(self.temp_roi_locations)) + " ROIs found")

        message = "ROIs loaded in. Be careful, the set variables for ROI finding and the used ROIs now no longer " \
                  "match. Therefore, using 'Strict' rejection rules is not advised. " \
                  "The loaded ROIs are automatically saved."

        tk.messagebox.showinfo("Done!", message)

    # %% Fitting page, switch between datasets

    def change_dataset(self, change):
        """
        Select next or previous dataset

        Parameters
        ------
        change: in what direction you want to change +1 or -1

        Returns
        -------
        None.

        """
        self.temp_roi_locations = None
        self.histogram = None
        self.to_hist = None

        self.dataset_index += change

        self.filename_short = filenames[self.dataset_index].split("/")[-1][:-4]
        self.dataset_roi_status.updater(text=self.filename_short + " (" + str(self.dataset_index + 1) + " of " +
                                        str(len(filenames)) + ")")

        self.nd2.close()
        self.nd2 = nd2_reading.ND2ReaderSelf(self.filenames[self.dataset_index])
        self.frames = self.nd2
        self.metadata = self.nd2.get_metadata()

        if self.dataset_index in self.saved_settings:
            self.restore_saved()
            self.button_restore_saved.updater(command=lambda: self.restore_saved())
            self.hsm_correct_var.set(self.saved_settings[self.dataset_index]['hsm_correction'])
            if self.saved_settings[self.dataset_index]['hsm_directory'] is not None:
                hsm_folder_show = '/'.join(self.saved_settings[self.dataset_index]['hsm_directory'].split('/')[-2:])
            else:
                hsm_folder_show = ""
            self.hsm_folder_disp['text'] = hsm_folder_show
        else:
            self.button_restore_saved.updater(state='disabled')
            self.clear_hsm()
            self.restore_default()

        if self.dataset_index == len(self.filenames) - 1:
            self.button_right.updater(state='disabled')
        else:
            self.button_right.updater(command=lambda: self.change_dataset(1))

        if self.dataset_index == 0:
            self.button_left.updater(state='disabled')
        else:
            self.button_left.updater(command=lambda: self.change_dataset(-1))

    # %% Fitting page, start fitting

    def start_fitting(self):
        """
        Start fitting. Takes all the inputs and fits.

        Returns
        -------
        None officially. Outputs files.

        """
        if len(self.roi_locations) != len(self.filenames):
            tk.messagebox.showerror("ERROR", "Not all datasets have settings yet, cannot start")
            return

        check = tk.messagebox.askokcancel("Are you sure?",
                                          "Fitting may take a while. Are you sure everything is set up correctly?")
        if not check:
            return

        n_processes = self.cores_var.get()
        method = self.method_var.get()
        rejection_type = self.rejection_var.get()
        figures_option = self.figures_var.get()

        if rejection_type == "Strict" and (method == "Phasor + Sum" or method == "Phasor" or
                                           method == "Phasor + Intensity"):
            rejection_check = tk.messagebox.askokcancel("Just a heads up",
                                                        """Phasor methods have no strict rejection option,
            so "Loose" rejection will be used""")
            if not rejection_check:
                return
            else:
                rejection_type = "Loose"

        if n_processes > 1 and (method == "Phasor + Intensity" or method == "Phasor + Sum" or method == "Phasor"):
            cores_check = tk.messagebox.askokcancel("Just a heads up",
                                                    """Phasor will be used with one core since the
            overhead only slows it down""")
            if not cores_check:
                return
            n_processes = 1

        self.start_time = time.time()
        results_counter = 0

        dataset_index_viewing = self.dataset_index

        for self.dataset_index, filename in enumerate(self.filenames):
            dataset_time = time.time()
            self.nd2.close()
            self.nd2 = nd2_reading.ND2ReaderSelf(filenames[self.dataset_index])
            self.frames = self.nd2
            self.metadata = self.nd2.get_metadata()

            dataset_settings = self.saved_settings[self.dataset_index]

            roi_size = dataset_settings['roi_size']
            max_its = dataset_settings['max_its']

            roi_locations = self.roi_locations[self.dataset_index]

            if method == "Phasor + Intensity":
                fitter = fitting.Phasor(roi_size, dataset_settings, rejection_type, method)
            elif method == "Phasor":
                fitter = fitting.PhasorDumb(roi_size, dataset_settings, rejection_type, method)
            elif method == "Gaussian - Fit bg":
                fitter = fitting.GaussianBackground(roi_size, dataset_settings, rejection_type, method, 6, max_its)
            elif method == "Gaussian - Estimate bg":
                fitter = fitting.Gaussian(roi_size, dataset_settings, rejection_type, method, 5, max_its)
            else:
                fitter = fitting.PhasorSum(roi_size, dataset_settings, rejection_type, method)

            start_frame = self.frame_begin_input.get()
            end_frame = self.frame_end_input.get()

            if start_frame == "Leave empty for start" and end_frame == "Leave empty for end":
                end_frame = self.metadata['num_frames']
                start_frame = 0
                to_fit = self.frames
                time_axis = self.metadata['timesteps']
            elif start_frame == "Leave empty for start" and end_frame != "Leave empty for end":
                start_frame = 0
                end_frame = int(end_frame)
                to_fit = self.frames[:end_frame]
                time_axis = self.metadata['timesteps'][:end_frame]
            elif start_frame != "Leave empty for start" and end_frame == "Leave empty for end":
                start_frame = int(start_frame)
                end_frame = self.metadata['num_frames']
                to_fit = self.frames[start_frame:]
                time_axis = self.metadata['timesteps'][start_frame:]
            else:  # start_frame != "Leave empty for start" and end_frame != "Leave empty for end":
                start_frame = int(start_frame)
                end_frame = int(end_frame)
                to_fit = self.frames[start_frame:end_frame]
                time_axis = self.metadata['timesteps'][start_frame:end_frame]

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
            else:
                self.update_status(0, num_frames)
                results = fitter.main(to_fit, self.metadata, roi_locations,
                                      gui=self, n_frames=num_frames)

            nm_or_pixels = self.dimension.get()
            if nm_or_pixels == "nm":
                results = tools.change_to_nm(results, self.metadata, method)

            # drift correction

            self.progress_status_label.updater(text="Drift correction dataset " +
                                                    str(self.dataset_index + 1) + " of " + str(len(filenames)))
            self.update()

            drift_corrector = drift_correction.DriftCorrector(method)
            results_drift, drift, event_or_not = drift_corrector.main(results, roi_locations, num_frames)

            # HSM

            hsm_dir = (dataset_settings['hsm_directory'],)
            hsm_corr = dataset_settings['hsm_correction']
            hsm_result = None  # just to make Python shut up about potential reference before assignment
            hsm_intensity = None
            hsm_wavelengths = None
            hsm_raw = None

            if hsm_dir[0] is not None and hsm_corr != '':
                self.progress_status_label.updater(text="HSM dataset " +
                                                        str(self.dataset_index + 1) + " of " + str(len(filenames)))

                self.update()

                hsm_object = hsm.HSM(hsm_dir, np.asarray(self.frames[0], dtype=self.frames[0].dtype),
                                     roi_locations.copy(), self.metadata, hsm_corr)
                hsm_result, hsm_raw, hsm_intensity = hsm_object.main(verbose=False)

                hsm_wavelengths = hsm_object.wavelength

            # create folder for output

            path = filename.split(".")[0]
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

            # Settings output

            total_fits = results.shape[0]
            failed_fits = results[np.isnan(results[:, 3]), :].shape[0]
            time_taken = round(time.time() - dataset_time, 3)

            successful_fits = total_fits - failed_fits
            results_counter += successful_fits

            outputting.text_output(dataset_settings.copy(), method, rejection_type, nm_or_pixels,
                                   total_fits, failed_fits, time_taken, path)

            # Plotting

            self.progress_status_label.updater(text="Plotting dataset " +
                                                    str(self.dataset_index + 1) + " of " + str(len(filenames)))
            self.update()
            try:
                figuring.save_graphs(self.frames, results, results_drift, roi_locations, method, nm_or_pixels,
                                     figures_option, path, event_or_not, dataset_settings, time_axis.copy(),
                                     hsm_result, hsm_raw, hsm_intensity, hsm_wavelengths)
            except Exception as _:
                show_error(False)

            # Switch to MATLAB coordinates

            results = tools.switch_results_to_matlab_coordinates(results, self.frames[0].shape[0],
                                                                 method, nm_or_pixels, self.metadata)
            results_drift = tools.switch_results_to_matlab_coordinates(results_drift, self.frames[0].shape[0],
                                                                       method, nm_or_pixels, self.metadata)
            roi_locations = tools.roi_to_matlab_coordinates(roi_locations, self.frames[0].shape[0])
            drift = tools.switch_axis(drift)
            if hsm_result is not None:
                hsm_result, hsm_raw, hsm_intensity = tools.switch_to_matlab_hsm(hsm_result, hsm_raw, hsm_intensity)

            # Save

            self.progress_status_label.updater(text="Saving dataset " +
                                                    str(self.dataset_index + 1) + " of " + str(len(filenames)))
            self.update()

            outputting.save_first_frame(self.frames[0], path)
            outputting.save_to_csv_mat_metadata('metadata', self.metadata, path)
            outputting.save_to_csv_mat_roi('ROI_locations', roi_locations, path)
            outputting.save_to_csv_mat_drift('Drift_correction', drift, path)
            outputting.save_to_csv_mat_results('Localizations', results, method, path)
            outputting.save_to_csv_mat_results('Localizations_drift', results_drift, method, path)
            if hsm_result is not None:
                outputting.save_hsm(hsm_result, hsm_raw, hsm_intensity, path)

        end_message = 'Time taken: ' + str(round(time.time() - self.start_time, 3)) \
                      + ' s. Fits done: ' + str(results_counter)

        self.progress_status_label.updater(text="Done!")
        tk.messagebox.showinfo("Done!", end_message)

        self.nd2.close()

        self.dataset_index = dataset_index_viewing

        self.nd2 = nd2_reading.ND2ReaderSelf(self.filenames[self.dataset_index])
        self.frames = self.nd2
        self.metadata = self.nd2.get_metadata()

    # %% Fitting page, update the status

    def update_status(self, progress, comparator):
        """
        Updates the status visible to the user.

        Parameters
        ----------
        progress : How far the fitter has come
        comparator : How far it has to go

        Returns
        -------
        None, updates GUI.

        """
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
            time_text = "{:02d}:{:02d}:{:02d} {:02d}/{:02d}".format(tr[3], tr[4], tr[5], tr[2], tr[1])

            self.time_status_label.updater(text=time_text)
        self.progress_status_label.updater(text=progress_text)

        self.update()

    # %% Histogram of sliders

    def fun_histogram(self, variable):
        """
        Actually makes the histogram

        Parameters
        ----------
        variable : Variable to make histogram of

        Returns
        -------
        None, outputs figure

        """
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

    def make_histogram(self, variable):
        """
        Makes histograms of parameters.

        Parameters
        ----------
        variable : The parameter to make a histogram of

        Returns
        -------
        None, output figure

        """
        fig_sub = self.histogram.add_subplot(111)
        hist, bins, _ = fig_sub.hist(self.to_hist, bins='auto')

        min_int = self.min_int_slider.get()
        max_int = self.max_int_slider.get()
        min_sigma = self.min_sigma_slider.get()
        max_sigma = self.max_sigma_slider.get()
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
        else:
            self.histogram.clear()
            logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
            plt.hist(self.to_hist, bins=logbins)
            plt.title("Correlation values for ROIs. Use graph select to change threshold")
            plt.axvline(x=min_corr, color='red', linestyle='--')
            plt.xscale('log')

        self.histogram.show()

    def histogram_select(self, variable):
        """
        Allows user to select from histogram to change slider

        Parameters
        ----------
        variable : variable to change

        Returns
        -------
        None.

        """
        try:
            click = self.histogram.ginput(1)
        except AttributeError:
            tk.messagebox.showerror("You cannot do that", "You cannot click on a figure without having a figure open")
            return

        if variable == "min_int":
            self.min_int_slider.updater(start=int(click[0][0]))
        elif variable == 'max_int':
            self.max_int_slider.updater(start=int(click[0][0]))
        elif variable == "min_sigma":
            self.min_sigma_slider.updater(start=click[0][0])
        elif variable == "max_sigma":
            self.max_sigma_slider.updater(start=click[0][0])
        else:
            self.min_corr_slider.updater(start=click[0][0])

        self.histogram.clear()
        self.make_histogram(variable)

    def load_hsm(self):
        tk.Tk().withdraw()
        hsm_folder = askdirectory(title="Select folder in which HSM data is located", initialdir=getcwd())

        if len(hsm_folder) == 0:
            return
        else:
            self.hsm_folder_full = hsm_folder
            hsm_folder_show = '/'.join(hsm_folder.split('/')[-2:])
            self.hsm_folder_disp['text'] = hsm_folder_show

    def clear_hsm(self):
        self.hsm_folder_full = None
        self.hsm_folder_disp['text'] = ""
        self.hsm_correct_var.set("")


# %% START GUI and declare styles (how things look)


if __name__ == '__main__':
    mp.freeze_support()
    gui = MbxPython()
    gui.geometry(str(GUI_WIDTH) + "x" + str(GUI_HEIGHT) + "+" + str(GUI_WIDTH_START) + "+" + str(GUI_HEIGHT_START))
    gui.iconbitmap(getcwd() + "\ico.ico")
    gui.protocol("WM_DELETE_WINDOW", lambda: quit_gui(gui))

    ttk_style = ttk.Style(gui)
    ttk_style.configure("Big.TButton", font=FONT_BUTTON_BIG)
    ttk_style.configure("Placeholder.TEntry", foreground="Grey")
    ttk_style.configure("TButton", font=FONT_BUTTON, background="Grey")
    ttk_style.configure("TSeparator", background="black")
    ttk_style.configure("TMenubutton", font=FONT_DROP, background="White")

    tk.Tk.report_callback_exception = show_error_critical

    plt.ioff()
    gui.mainloop()
