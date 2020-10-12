# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
MBx Python Data Analysis

main_GUI

This package is for the GUI of Mbx Python.

----------------------------

v1.0: roll-out version one
v1.1: Bugfixes and improved figures (WIP)
v1.2: GUI and output improvement based on Sjoerd's feedback, HSM: 27/08/2020 - 13/09/2020
v1.3: HSM to eV: 24/09/2020
v1.4: HSM output back to nm, while fitting in eV: 29/09/2020
"""
__version__ = "2.0"
__self_made__ = True

# GENERAL IMPORTS
from os import getcwd, mkdir, environ, listdir  # to get standard usage
from tempfile import mkdtemp
import sys
import time  # for timekeeping
from win32api import GetSystemMetrics  # Get sys info
import warnings  # for warning diversion

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
from tkinter.filedialog import askopenfilename, askdirectory  # for popup that asks to select .nd2's or folders

# Own code v2
from main import DivertError
from src.class_experiment import Experiment
import src.figure_making as figuring
from src.warnings import InputWarning

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

rejection_options = ["Loose", "None"]

roi_size_options = ["7x7", "9x9"]


# %% Multiprocessing main

# TO DO

# %% Divert errors


class DivertorErrorsGUI(DivertError):
    @staticmethod
    def show(error, traceback_details):
        if error:
            tk.messagebox.showerror("Critical error. Send screenshot to Dion. PROGRAM WILL STOP",
                                    message=str(traceback_details))
        else:
            tk.messagebox.showerror("Warning. Take note. PROGRAM WILL CONTINUE",
                                    message=str(traceback_details))

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
    Controller of GUI. This container calls the page we need
    """

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "MBx Python")
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        frame_tuple = (MainPage, LoadPage, ROIPage, TTPage, HSMPage)

        for to_load_frame in frame_tuple:
            frame = to_load_frame(container, self)
            self.frames[to_load_frame] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_page(MainPage)

        self.experiments = []
        self.list_of_datasets = []
        self.experiment_to_link_name = None

    def show_page(self, page):
        frame = self.frames[page]
        frame.tkraise()

# %% Main page base


class MainPageBase(tk.Frame):
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        self.configure(bg='white')
        self.controller = controller

        label_version = tk.Label(self, text="MBx Python, version: " + __version__, font=FONT_LABEL, bg='white',
                                 anchor='w')
        label_version.grid(row=50, column=0, columnspan=20, sticky='EW', padx=PAD_SMALL)

        button_quit = ttk.Button(self, text="Quit", command=lambda: quit_gui(self.controller))
        button_quit.grid(row=50, column=43, columnspan=6, sticky='EW', padx=PAD_SMALL)

    def column_row_configure(self):
        for i in range(49):
            self.grid_columnconfigure(i, weight=1)
        for i in range(51):
            self.grid_rowconfigure(i, weight=1)


# %% Base page


class BasePage(MainPageBase):
    def __init__(self, container, controller):
        super().__init__(container, controller)

        button_cancel = ttk.Button(self, text="Cancel", command=lambda: self.cancel())
        button_cancel.grid(row=50, column=37, columnspan=6, sticky='EW', padx=PAD_SMALL)

    def cancel(self):
        self.controller.show_page(MainPage)

# %% Main page


class MainPage(MainPageBase):
    def __init__(self, container, controller):
        super().__init__(container, controller)

        label_new = tk.Label(self, text="New", font=FONT_HEADER, bg='white')
        label_new.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        button_new_experiment = BigButton(self, text="ADD EXPERIMENT", height=int(GUI_HEIGHT / 4),
                                          width=int(GUI_WIDTH / 4), command=lambda: self.add_experiment())
        button_new_experiment.grid(row=0, column=0)

        button_new_dataset = BigButton(self, text="ADD DATASET", height=int(GUI_HEIGHT / 4),
                                       width=int(GUI_WIDTH / 4), command=lambda: self.add_dataset())
        button_new_dataset.grid(row=1, column=0)

        label_loaded = tk.Label(self, text="Loaded", font=FONT_HEADER, bg='white')
        label_loaded.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.listbox_loaded = tk.Listbox(self)
        self.listbox_loaded.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        button_loaded_delete = ttk.Button(self, text="Delete", command=lambda: self.delete_experiment())
        button_loaded_delete.grid(row=2, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        button_loaded_deselect = ttk.Button(self, text="Delete", command=lambda: self.deselect_experiment())
        button_loaded_deselect.grid(row=2, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        label_queued = tk.Label(self, text="Queued", font=FONT_HEADER, bg='white')
        label_queued.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.listbox_queued = tk.Listbox(self)
        self.listbox_queued.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        button_run = ttk.Button(self, text="Run", command=lambda: self.run())
        button_run.grid(row=2, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        label_progress_task = tk.Label(self, text="Task Progress", font=FONT_HEADER, bg='white')
        label_progress_task.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        label_progress_overall = tk.Label(self, text="Overall Progress", font=FONT_HEADER, bg='white')
        label_progress_overall.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.label_progress_task_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                      row=24, column=45, columnspan=5, sticky="ew", font=FONT_LABEL)
        self.label_progress_overall_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                         row=24, column=45, columnspan=5, sticky="ew", font=FONT_LABEL)

        self.label_current_task = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                              row=24, column=45, columnspan=5, sticky="ew", font=FONT_LABEL)

        label_time_done = tk.Label(self, text="Time done", font=FONT_HEADER, bg='white')
        label_time_done.grid(row=0, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.label_time_done_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                  row=24, column=45, columnspan=5, sticky="ew", font=FONT_LABEL)

        # general setup
        self.column_row_configure()

    def add_experiment(self):
        self.controller.show_page(LoadPage)
        self.controller.experiment_to_link_name = None

    def add_dataset(self):
        self.controller.show_page(LoadPage)
        self.controller.experiment_to_link_name = "Test"  # TO DO

    def delete_experiment(self):
        pass

    def deselect_experiment(self):
        pass

    def run(self):
        pass

# %% Loading page


class LoadPage(BasePage):
    """
    Loading page. On this page, there are only two big buttons to select which type of dataset you want to load
    """
    def __init__(self, container, controller):
        super().__init__(container, controller)

        button1 = BigButton(self, text="TT", height=int(GUI_HEIGHT / 4),
                            width=int(GUI_WIDTH / 4),  # style= 'Big.TButton',
                            command=lambda: self.load_nd2("TT"))
        button1.grid(row=0, column=0)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        button2 = BigButton(self, text="HSM", height=int(GUI_HEIGHT / 4),
                            width=int(GUI_WIDTH / 4),  # style= 'Big.TButton',
                            command=lambda: self.load_nd2("HSM"))
        button2.grid(row=0, column=1)

        # general setup
        self.column_row_configure()

    def load_nd2(self, dataset_type):
        filename = askopenfilename(filetypes=FILETYPES,
                                   title="Select nd2",
                                   initialdir=getcwd())

        if len(filename) == 0:
            return

        if self.controller.experiment_to_link_name is None:
            experiment = Experiment(dataset_type, filename, self.controller.proceed_question,
                                    self.controller.progress_updater, self.controller.show_rois)
            self.controller.experiments.append(experiment)
            self.controller.show_page(ROIPage)
        else:
            experiment_to_link = [experiment for experiment in self.controller.experiments if
                                  self.controller.experiment_to_link_name in experiment.name][0]
            if dataset_type == "TT":
                experiment_to_link.init_new_tt(filename)
                self.controller.show_page(TTPage)
            else:
                experiment_to_link.init_new_hsm(filename)
                self.controller.show_page(HSMPage)

# %% ROIPage


class ROIPage(BasePage):
    def __init__(self, container, controller):
        super().__init__(container, controller)

# %% TTPage


class TTPage(BasePage):
    def __init__(self, container, controller):
        super().__init__(container, controller)

# %% HSMPage


class HSMPage(BasePage):
    def __init__(self, container, controller):
        super().__init__(container, controller)

# %% START GUI and declare styles (how things look)


if __name__ == '__main__':
    mp.freeze_support()
    divertor = DivertorErrorsGUI()
    warnings.showwarning = divertor.warning
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

    tk.Tk.report_callback_exception = divertor.error

    plt.ioff()
    gui.mainloop()
