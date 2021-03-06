# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
PLASMON Data Analysis

main_GUI

This package is for the GUI of PLASMON.

----------------------------

v1.0: roll-out version one
v1.1: Bugfixes and improved figures (WIP)
v1.2: GUI and output improvement based on Sjoerd's feedback, HSM: 27/08/2020 - 13/09/2020
v1.3: HSM to eV: 24/09/2020
v1.4: HSM output back to nm, while fitting in eV: 29/09/2020
v2.0 pre-1: First version of GUI v2.0: 15/10/2020
v2.0: GUI v2.0 ready for release: 30/10/2020
"""

__self_made__ = True

# GENERAL IMPORTS
from os import getcwd, environ, listdir, rmdir, mkdir, path  # to get standard usage
from tempfile import mkdtemp
import sys
import time  # for timekeeping
import ctypes  # Get sys info
import logging  # for logging warnings
from traceback import format_exc  # for error handling
import shelve

environ['MPLCONFIGDIR'] = mkdtemp()

# Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# GUI
import tkinter as tk  # for GUI
from tkinter import ttk  # GUI styling
from tkinter.filedialog import askopenfilename  # for popup that asks to select .nd2's or folders

# Own code
from main import ProgressUpdater, logging_setup
from src.class_experiment import Experiment
import src.figure_making as figuring
from src.nd2_reading import ND2ReaderSelf
from setup import __version__

# Multiprocessing
import multiprocessing as mp
import _thread

mpl.use("TkAgg")  # set back end to TK

# %% Initializations. Defining filetypes, fonts, paddings, input sizes, and GUI sizes.

FILETYPES = [("ND2", ".nd2")]
FILETYPES_LOAD_FROM_OTHER = [(".npy and .mat", ".npy"), (".npy and .mat", ".mat")]

user32 = ctypes.windll.user32
SCREEN_WIDTH = user32.GetSystemMetrics(0)
SCREEN_HEIGHT = user32.GetSystemMetrics(1)
DPI = 100

if SCREEN_WIDTH > 1280 and SCREEN_HEIGHT > 720:
    # 1080p version
    GUI_WIDTH = 1440  # int(width * 0.75)
    GUI_HEIGHT = 810  # int(height * 0.75)
    FONT_HEADER = "Verdana 14 bold"
    FONT_SUBHEADER = "Verdana 12 bold"
    FONT_STATUS = "Verdana 11"
    FONT_ENTRY = "Verdana 11"
    FONT_ENTRY_SMALL = "Verdana 9 "
    FONT_BUTTON = "Verdana 11"
    FONT_LABEL = "Verdana 11"
    FONT_DROP = "Verdana 11"
    FONT_LISTBOX = "Verdana 8"
    FONT_BUTTON_BIG = "Verdana 20 bold"
    FONT_TOOLTIP = "Verdana 8"
    PAD_BIG = 30
    PAD_SMALL = 10
    INPUT_BIG = 25
    INPUT_MEDIUM = 10
    INPUT_SMALL = 5
else:
    # 720p version
    GUI_WIDTH = 1152  # int(width * 0.90)
    GUI_HEIGHT = 648  # int(height * 0.90)
    FONT_HEADER = "Verdana 12 bold"
    FONT_SUBHEADER = "Verdana 10 bold"
    FONT_STATUS = "Verdana 9"
    FONT_ENTRY = "Verdana 9"
    FONT_ENTRY_SMALL = "Verdana 7"
    FONT_BUTTON = "Verdana 9"
    FONT_LABEL = "Verdana 9"
    FONT_DROP = "Verdana 9"
    FONT_LISTBOX = "Verdana 6"
    FONT_BUTTON_BIG = "Verdana 16 bold"
    FONT_TOOLTIP = "Verdana 6"
    PAD_BIG = 20
    PAD_SMALL = 5
    INPUT_BIG = 15
    INPUT_MEDIUM = 5
    INPUT_SMALL = 3

GUI_WIDTH_START = int((SCREEN_WIDTH - GUI_WIDTH) / 2)
GUI_HEIGHT_START = int((SCREEN_HEIGHT - GUI_HEIGHT) / 2)

# %% Options for dropdown menus

fit_options = ["Gaussian - Fit bg", "Gaussian - Estimate bg",
               "Phasor + Intensity", "Phasor + Sum", "Phasor"]
roi_size_options = ["7x7", "9x9"]
dimension_options = ["nm", "pixels"]

# %% Proceed Question


def proceed_question(title, text):
    """
    Asks user to proceed or not
    :param title: Title
    :param text: Text
    :return: True or False depending on proceed or not
    """
    check = tk.messagebox.askokcancel(title, text)
    return check

# %% Logging


class GUIStream(logging.Handler):

    def emit(self, record):
        msg = self.format(record)
        tk.messagebox.showwarning("Attention!", msg)


# %% Divert errors


class DivertorErrorsGUI:
    """
    GUI version of DivertorError
    """
    def __init__(self, logger):
        self.logger = logger

    def error(self, *params):
        """
        Called when error is given
        :param: All unused error info
        :return: Calls show
        """
        # try to log. If fails, simply show on screen
        try:
            traceback_details = self.extract_error()
            logger.error(f"Critical error occurred! PROGRAM WILL STOP.\nCheck logging file for more info.\n"
                         f"Short summary: {traceback_details}")
            logger.info("More info about critical error: " + format_exc(10, params[2]))
        except:
            self.show(format_exc(10, params[2]))

    @staticmethod
    def extract_error():
        """
        Gets error info in detail
        :return: traceback_details. Dict with info about error
        """
        exc_type, exc_value, exc_traceback = sys.exc_info()
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
        return traceback_details

    @staticmethod
    def show(traceback_details):
        """
        Shows the actual error or warning in Tkinter. Only used if logging fails
        :param traceback_details: Error details
        :return: Prints out
        """
        tk.messagebox.showerror("Error. Send screenshot to Dion. PROGRAM WILL STOP",
                                message=f"Logging failed. Full error:\n{traceback_details}")


# %% Close GUI


def quit_gui(gui):
    """
    Quits the GUI, clear empty directories created, and saved last opened directory
    :param gui: GUI to close
    :return: None. Closes program
    """
    # for all loaded experiments, try to delete directory if directory is made. This only succeeds if its empty
    for experiment in gui.experiments:
        if experiment.dir_made:
            try:
                rmdir(experiment.directory)
            except:
                pass
    # store last opened directory
    try:
        my_shelve = shelve.open(gui.shelve_file, 'n')
        my_shelve['dir_open'] = gui.dir_open
        my_shelve.close()
    except:
        pass
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

    def updater(self, frame, roi_locations=None, roi_size=None, roi_offset=None, overwrite=False):
        """
        Updater. Takes existing frame with figure and places new figure in it

        Parameters
        ----------
        frame : New frame to be shown
        roi_locations : optional, possible ROI locations to be highlighted. The default is None.
        roi_size : optional, ROI size in case ROIs are highlighted. The default is None.
        roi_offset: offset compared to ROI frame
        overwrite: Only called by ROI page, when the figure needs to be updated with new ROIs but same frame.

        Returns
        -------
        Updated figure.

        """
        if roi_offset is None:
            roi_offset = [0, 0]
        # if you want to overwrite, clear patches and place new ones
        if overwrite:
            fig_sub = self.fig.axes[0]
            fig_sub.patches = []
            figuring.plot_rois(fig_sub, frame, roi_locations, roi_size, roi_offset, overwrite=True)
        else:
            # else make new figure
            self.fig.clear()
            fig_sub = self.fig.add_subplot(111)
            figuring.plot_rois(fig_sub, frame, roi_locations, roi_size, roi_offset, overwrite=False)
        # always draw
        self.canvas.draw()
        # only update toolbar once
        if not overwrite:
            self.toolbar.update()


class EntryPlaceholder(ttk.Entry):
    """
    Entry with a placeholder text in grey
    """
    def __init__(self, master=None, placeholder="PLACEHOLDER", small=False, *args, **kwargs):
        if small:
            self.placeholder_style = "PlaceholderSmall.TEntry"
            self.normal_style = "Small.TEntry"
            self.font = FONT_ENTRY_SMALL
        else:
            self.placeholder_style = "Placeholder.TEntry"
            self.normal_style = "TEntry"
            self.font = FONT_ENTRY

        super().__init__(master, *args, style=self.placeholder_style, font=self.font, **kwargs)
        self.placeholder = placeholder

        self.insert("0", self.placeholder)
        self.bind("<FocusIn>", self._clear_placeholder)
        self.bind("<FocusOut>", self._add_placeholder)

    def _clear_placeholder(self, _):
        self.delete("0", "end")
        if self["style"] == self.placeholder_style:
            self["style"] = self.normal_style

    def _add_placeholder(self, _):
        if not self.get():
            self.insert("0", self.placeholder)
            self["style"] = self.placeholder_style

    def updater(self, text=None, placeholder=None):
        """
        Update Placeholder with new text or new placeholder
        :param text: new text
        :param placeholder: new placeholder
        :return: Updates Entry
        """
        self.delete("0", "end")
        if placeholder is not None:
            self.placeholder = placeholder

        if text is None:
            self.insert("0", self.placeholder)
            self["style"] = self.placeholder_style
        else:
            self.insert("0", text)
            self["style"] = self.normal_style


class EntrySlider(ttk.Entry):
    """
    Entry to be linked to a slider
    """
    def __init__(self, master=None, row=None, column=None, rowspan=1, columnspan=1, *args, **kwargs):
        super().__init__(master, *args, style="Small.TEntry", font=FONT_ENTRY_SMALL, **kwargs)
        self.grid(row=row, column=column, rowspan=rowspan, columnspan=columnspan)

    def binder(self, slider, variable):
        """
        Bind a slider to a variable
        :param slider: slider to couple
        :param variable: variable to couple
        """
        def change_slider(_):
            # if it is a number, it can become a int. If not, this will fail, and nothing will happen.
            try:
                int(variable.get())
                slider.updater(start=variable.get())
            except:
                tk.messagebox.showerror("Error", "Please enter a number.")

        self.bind("<Return>", change_slider)
        self.bind("<FocusOut>", change_slider)


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
        """
        Updates button with new command, state, or text
        :param command: new command
        :param state: new state
        :param text: new text
        :return: Updates button
        """
        if text is None:
            text = self.text
        self._btn['text'] = text
        self._btn['state'] = state
        self._btn['command'] = command


class NormalSlider:
    """
    My normal slider, again with an updater function to update the slider.
    """
    def __init__(self, parent, from_=0, to=np.inf, resolution=1, start=0, command=None,
                 row=None, column=None, rowspan=1, columnspan=1, sticky=None, padx=0, pady=0):
        self._scale = tk.Scale(parent, from_=from_, to=to, orient='horizontal', command=command, showvalue=0,
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
        """
        Updates scale with new from, to, and/or start
        :param from_: new from
        :param to: new to
        :param start: new start
        :return: Updated scale
        """
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
                 row=None, column=None, rowspan=1, columnspan=1, sticky=None, padx=0, pady=0, wrap=None):
        self._label = tk.Label(parent, text=text, font=font, bd=bd, relief=relief, bg='white', wrap=wrap)
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
        """
        Update label with new text
        :param text: new text
        :return: updated label
        """
        if text is None:
            text = self.text
        self._label['text'] = text
        # update page
        self.parent.update()

    def grid(self):
        self._label.grid()
        self.parent.update()

    def grid_remove(self):
        self._label.grid_remove()
        self.parent.update()

# %% Progress Updater


class ProgressUpdaterGUI(ProgressUpdater):
    """
    GUI version of ProgressUpdater
    """
    def __init__(self, gui, progress_task_status, progress_overall_status, current_task_status, time_done_status):
        """
        Initializer of ProgressUpdaterGUI. Also adds GUI and GUI elements to object
        :param gui: GUI called by
        :param progress_task_status: label with task status
        :param progress_overall_status: label with overall status
        :param current_task_status: label with current status
        :param time_done_status: label with eta
        """
        super().__init__()
        self.gui = gui
        self.style = ttk.Style(gui)
        self.progress_task_status = progress_task_status
        self.progress_overall_status = progress_overall_status
        self.current_task_status = current_task_status
        self.time_done_status = time_done_status
        self.start_time = time.time()

    def update(self, new_experiment, new_dataset, message_bool):
        """
        Update function. Does the actual communication.
        :param new_experiment: Boolean. Whether or not new experiment
        :param new_dataset: Boolean. Whether or not new dataset
        :param message_bool: Boolean. Whether or not message
        :return: prints out info in tkinter
        """
        progress_dataset = (self.current_dataset - 1) / self.total_datasets
        progress_per_dataset = 1 / self.total_datasets
        if new_experiment:
            pass
        elif new_dataset:
            overall_progress = progress_dataset

            self.progress_task_status.updater(text="0%")
            self.progress_overall_status.updater(text="{:.2f}%".format(overall_progress * 100))

            self.current_task_status.updater(text="Experiment #{}: Dataset #{}: {}".format(self.current_experiment + 1,
                                                                                           self.current_dataset,
                                                                                           self.current_type))
        elif message_bool:
            progress_task = 1
            progress_overall = progress_task * progress_per_dataset + progress_dataset

            self.progress_task_status.updater(text="{:.2f}%".format(progress_task * 100))
            self.progress_overall_status.updater(text="{:.2f}%".format(progress_overall * 100))

            self.current_task_status.updater(text="Experiment {}: ".format(self.current_experiment + 1) +
                                                  self.message_string)
        else:
            progress_task = self.progress / self.total / self.dataset_parts
            progress_overall = progress_task * progress_per_dataset + progress_dataset

            self.progress_task_status.updater(text="{:.2f}%".format(progress_task*100))
            self.progress_overall_status.updater(text="{:.2f}%".format(progress_overall*100))

            if progress_overall != 0:
                time_taken = time.time() - self.start_time
                time_done_estimate = time_taken * 1 / progress_overall + self.start_time
                time_done_estimate += 20  # 20 seconds for overview figure
                if self.experiment_ind_figures is True:
                    time_done_estimate += 1.5 * self.experiment_rois  # add 1.5 seconds per individual figure

                tr = time.localtime(time_done_estimate)
                time_text = "{:02d}:{:02d}:{:02d} {:02d}/{:02d}".format(tr[3], tr[4], tr[5], tr[2], tr[1])
                self.time_done_status.updater(text=time_text)
            else:
                self.time_done_status.updater(text="TBD")

        self.gui.update()

# %% Footer


class FooterBase(tk.Frame):
    """
    FooterBase class. Shows close and version number
    """
    def __init__(self, controller):
        tk.Frame.__init__(self, controller)
        self.configure(bg='white')
        self.controller = controller

        label_version = tk.Label(self, text="PLASMON v{}".format(__version__), font=FONT_LABEL, bg='white',
                                 anchor='w')
        label_version.grid(row=0, column=0, columnspan=20, sticky='EW', padx=PAD_SMALL)

        button_resize = ttk.Button(self, text="Reset window size", command=lambda: self.resize())
        button_resize.grid(row=0, column=40, columnspan=4, sticky='EW', padx=PAD_SMALL)

        button_quit = ttk.Button(self, text="Quit", command=lambda: quit_gui(self.controller))
        button_quit.grid(row=0, column=44, columnspan=4, sticky='EW', padx=PAD_SMALL)

        for i in range(48):
            self.grid_columnconfigure(i, weight=1)

    def resize(self):
        self.controller.geometry("{}x{}+{}+{}".format(GUI_WIDTH, GUI_HEIGHT, GUI_WIDTH_START, GUI_HEIGHT_START))


class Footer(FooterBase):
    """
    Footer class. Also shows cancel button
    """
    def __init__(self, controller):
        super().__init__(controller)

        button_cancel = ttk.Button(self, text="Cancel", command=lambda: self.cancel())
        button_cancel.grid(row=0, column=40, columnspan=4, sticky='EW', padx=PAD_SMALL)

        for i in range(48):
            self.grid_columnconfigure(i, weight=1)

    def cancel(self):
        """
        Cancel current page. Resets to MainPage and clears the page if needed
        :return: None. Changes page
        """
        # if not on load page
        if self.controller.current_page != LoadPage:
            # clear current page
            self.controller.pages[self.controller.current_page].clear_page()
            # if new experiment, clear directory and experiment
            if self.controller.experiment_to_link_name is None:
                if self.controller.experiments[-1].dir_made:
                    rmdir(self.controller.experiments[-1].directory)
                del self.controller.experiments[-1]
            else:
                # else clear last dataset
                experiment_to_link = [experiment for experiment in self.controller.experiments if
                                      self.controller.experiment_to_link_name in experiment.name][0]
                del experiment_to_link.datasets[-1]
        # always show MainPage
        self.controller.show_page(MainPage)


# %% Tooltips
TOOLTIP_MAIN_PROGRESS_TASK = "The progress of the current task"
TOOLTIP_MAIN_PROGRESS_OVERALL = "The overall progress of current analysis queue"
TOOLTIP_MAIN_CURRENT_TASK = "The task that is currently being analysed"
TOOLTIP_MAIN_TIME_DONE = "A rough estimate of when the analysis should be done.\nTake with a grain of salt."
TOOLTIP_ROI_NAME_EXPERIMENT = "The name that will be given to the experiment.\n" \
                              "Defaults to the name of the first dataset selected"
TOOLTIP_ROI_MIN_INTENSITY = "The minimum intensity of a particle to be valid to be selected.\n" \
                            "Found by fitting a 2D Gaussian to the particle."
TOOLTIP_ROI_MAX_INTENSITY = "The maximum intensity of a particle to be valid to be selected.\n" \
                            "Found by fitting a 2D Gaussian to the particle."
TOOLTIP_ROI_MIN_SIGMA = "The minimum sigma of a particle to be valid to be selected.\n" \
                            "Found by fitting a 2D Gaussian to the particle."
TOOLTIP_ROI_MAX_SIGMA = "The maximum sigma of a particle to be valid to be selected.\n" \
                            "Found by fitting a 2D Gaussian to the particle."
TOOLTIP_ROI_SETTINGS = "General settings for finding the correct ROIs around particles.\n" \
                       "A ROI has to pass all tests otherwise it is not selected."
TOOLTIP_ROI_ADVANCED_SETTINGS = "Advanced settings for finding the correct ROIS around particles.\n" \
                                "These might not be as intuitive as the settings above.\n" \
                                "A ROI has to pass all tests otherwise it is not selected."
TOOLTIP_ROI_CORRELATION = "The minimum correlation with a 2D Gaussian that a particle has to have.\n" \
                          "If set too low, ROIs will be found on noise.\n" \
                          "A ROI has to pass all tests otherwise it is not selected."
TOOLTIP_ROI_EDGE_DIST = "Minimum distance (pixels) between a ROI center and the edge of the FOV.\n"  \
                        "A ROI has to pass all tests otherwise it is not selected."
TOOLTIP_ROI_INTER_DIST = "Minimum distance (pixels) between different ROI centers.\n" \
                         "A ROI has to pass all tests otherwise it is not selected."
TOOLTIP_ROI_OTHER_SETTINGS = "Other settings, not directly related to finding ROIs\n" \
                             "but will be applied to the entire experiment"
TOOLTIP_ROI_INDIVIDUAL_FIGURES = "For each experiment, an overview figure will be created.\n" \
                                 "If this is checked, individual figures for each ROI will also be created.\n" \
                                 "Note, this can be slow."
TOOLTIP_ROI_BACKGROUND_FILTER = "The size used by the background median filter.\n9 works well for most datasets."
TOOLTIP_ANALYSIS_NAME_LOADED = "Name of the loaded ND2 file."
TOOLTIP_ANALYSIS_NAME_DATASET = "Name given to the dataset. Defaults to the name of the loaded ND2.\n"
TOOLTIP_ANALYSIS_PRECROP = "Here you have to find the ROIs found on the main dataset of the experiment for the " \
                           "new dataset.\nFor the first dataset, this is trivial, however, you still have to do it." \
                           "\nThe entry fields below allow you to crop the search area, helping the algorithm find " \
                           "where the particles are present in the new dataset."
TOOLTIP_ANALYSIS_MIN_X = "The minimum x-value to pre-crop with.\n" \
                         "This means that you know that the particles are not below this x-value."
TOOLTIP_ANALYSIS_MAX_X = "The maximum x-value to pre-crop with.\n" \
                         "This means that you know that the particles are not above this x-value."
TOOLTIP_ANALYSIS_MIN_Y = "The minimum y-value to pre-crop with.\n" \
                         "This means that you know that the particles are not below this y-value."
TOOLTIP_ANALYSIS_MAX_Y = "The maximum y-value to pre-crop with.\n" \
                         "This means that you know that the particles are not above this y-value."
TOOLTIP_TT_MAIN = "All the settings related to TT analysis."
TOOLTIP_TT_METHOD = "Method to fit the ROIs with.\nPhasor will only give location (and whatever is also shown)," \
                    "while Gaussian will fit a 2D Gaussian."
TOOLTIP_TT_REJECTION = "Enable rejection or not.\nIf rejection is enabled, " \
                       "impossible values (such as negative intensities, very large sigmas) are rejected."
TOOLTIP_TT_CORES = "The number of cores used for analysis. More is faster.\n" \
                   "The maximum makes your pc unusable while running. 3rd option is advised for fast analysis"
TOOLTIP_TT_PIXELS_OR_NM = "Do you want your output in pixels or in nm?\nPixelsize can be gathered from video metadata."
TOOLTIP_TT_CORRELATION_INTERVAL = "If you expect a lot of drift,\nyou can have the program adjust the ROI positions.\n" \
                                  "This will be done every X frame (with X being what you set below).\n" \
                                  "At this point, the program wil correlate the entire FOV\n" \
                                  " at that frame with X frames earlier,\nadjusting the ROI positions if need be.\n" \
                                  "This way, the particles will never leave the ROIs.\n" \
                                  "In short, a low value will prevent particles leaving ROIs,\n" \
                                  "but will slow the program down."
TOOLTIP_TT_USED_ROI_SPACING = "The ROI spacing you set for the experiment. For your information."
TOOLTIP_TT_ROI_SIZE = "The ROI size you want to use. 7x7 works usually, and is fastest.\n" \
                      "Larger ROI size is slower, but might work better if you expected large PSFs."
TOOLTIP_TT_FIRST_FRAME = "First frame to fit for the TT.\nYou can use this to crop the video."
TOOLTIP_TT_LAST_FRAME = "Last frame to fit for the TT.\nYou can use this to crop the video."
TOOLTIP_HSM_MAIN = "All the settings related to HSM analysis."
TOOLTIP_HSM_CORRECTION_FILE = "The correction file to use for HSM."
TOOLTIP_HSM_WAVELENGTHS = "The wavelengths that were used to created the HSM.\n" \
                          "Use MATLAB-like notation as shown in the placeholder"


class ToolTip:

    def __init__(self, widget):
        self.widget = widget
        self.tip_window = None
        self.text = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        """Display text in tooltip window"""
        self.text = text
        if self.tip_window or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() + 27
        self.tip_window = tk.Toplevel(self.widget)
        self.tip_window.wm_overrideredirect(1)
        self.tip_window.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tip_window, text=self.text, justify=tk.LEFT,
                         background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                         font=FONT_TOOLTIP)
        label.pack(ipadx=1)

    def hidetip(self):
        if self.tip_window:
            self.tip_window.destroy()
        self.tip_window = None


def create_tooltip(widget, text):
    tooltip = ToolTip(widget)

    def enter(_):
        tooltip.showtip(text)

    def leave(_):
        tooltip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

# %% Controller


class PLASMON(tk.Tk):
    """
    Controller of GUI. This container calls the pages we need and store overall data
    """
    def __init__(self, proceed_question=None, *args, **kwargs):
        """
        Initializer of total GUI
        :param proceed_question: Proceed Question to use
        :param args: other arguments
        :param kwargs: other arguments
        """
        tk.Tk.__init__(self, *args, **kwargs)
        # withdraw until done
        self.withdraw()

        # setup container
        tk.Tk.wm_title(self, "PLASMON")
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)

        # and setup footer
        self.footer = FooterBase(self)
        self.footer.pack(side="bottom", fill="both")

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # setup data
        self.proceed_question = proceed_question
        self.progress_updater = None
        self.experiments = []
        self.experiment_to_link_name = None
        self.thread_started = False

        # load in old working directory
        # check location in Documents
        userprofile = environ['USERPROFILE']
        dir_path = path.join(userprofile, "Documents", 'PLASMON')
        # check dir path
        if path.exists(dir_path):
            pass
        else:
            mkdir(dir_path)

        self.shelve_file = dir_path + "/Cache.out"
        try:
            my_shelve = shelve.open(self.shelve_file)
            self.dir_open = my_shelve['dir_open']
            my_shelve.close()
        except:
            self.dir_open = getcwd()

        # setup pages
        self.pages = {}
        page_tuple = (MainPage, LoadPage, ROIPage, TTPage, HSMPage)
        for to_load_page in page_tuple:
            page = to_load_page(container, self)
            self.pages[to_load_page] = page
            page.grid(row=0, column=0, sticky="nsew")
        self.current_page = MainPage
        self.show_page(MainPage)

        # additional settings and pop-up again
        self.additional_settings()
        self.deiconify()

    @staticmethod
    def show_rois(frame, figure=None, roi_locations=None, roi_size=None, roi_offset=None, overwrite=False):
        """
        Shows ROIs within python
        :param frame: frame to make figure of
        :param figure: figure (only used by GUI)
        :param roi_locations: ROI locations within frame
        :param roi_size: ROI size
        :param roi_offset: Offset of ROIs within dataset
        :param overwrite: overwrite the previous figure or not. Basically, if true, it keeps the zoom
        :return:
        """
        if figure is None:
            figure = plt.subplots(1)
        figure.updater(frame, roi_locations=roi_locations, roi_size=roi_size, roi_offset=roi_offset, overwrite=overwrite)

    def additional_settings(self):
        """
        Sets additional settings, such as styles and icon.
        :return: None. Edits GUI
        """
        self.geometry(str(GUI_WIDTH) + "x" + str(GUI_HEIGHT) + "+" + str(GUI_WIDTH_START) + "+" + str(GUI_HEIGHT_START))
        self.iconbitmap(getcwd() + "\ico.ico")
        self.protocol("WM_DELETE_WINDOW", lambda: quit_gui(gui))

        ttk_style = ttk.Style(self)
        ttk_style.configure("Big.TButton", font=FONT_BUTTON_BIG)
        ttk_style.configure("Placeholder.TEntry", foreground="Grey", font=FONT_ENTRY)
        ttk_style.configure("TEntry", font=FONT_ENTRY)
        ttk_style.configure("PlaceholderSmall.TEntry", foreground="Grey", font=FONT_ENTRY_SMALL)
        ttk_style.configure("Small.TEntry", font=FONT_ENTRY_SMALL)
        ttk_style.configure("TButton", font=FONT_BUTTON, background="White")
        ttk_style.configure("TSeparator", background="black")
        ttk_style.configure("TMenubutton", font=FONT_DROP, background="White")
        ttk_style.configure("TCheckbutton", background="White")

    def show_page(self, page, experiment=None):
        """
        Show other page
        :param page: Page to show
        :param experiment: experiment to sent to page if need be
        :return: None. Changes page
        """
        self.current_page = page
        if page == MainPage:
            self.footer.pack_forget()
            self.footer = FooterBase(self)
            self.footer.pack(side="bottom", fill="both")
        else:
            self.footer.pack_forget()
            self.footer = Footer(self)
            self.footer.pack(side="bottom", fill="both")
        page = self.pages[page]
        page.update_page(experiment=experiment)
        page.tkraise()

# %% Base page


class BasePage(tk.Frame):
    """
    BasePage with some standard initializations
    """
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        self.configure(bg='white')
        self.controller = controller

        self.column_row_configure()

    def column_row_configure(self):
        """
        Standard function for column and row scaling
        """
        for i in range(48):
            self.grid_columnconfigure(i, weight=1, minsize=18)
        for i in range(20):
            self.grid_rowconfigure(i, weight=1)

    def update_page(self, experiment=None):
        pass

# %% Main page


class MainPage(BasePage):
    def __init__(self, container, controller):
        """
        Initializes MainPage, including all buttons
        :param container: Container
        :param controller: Controller
        """
        super().__init__(container, controller)

        label_new = tk.Label(self, text="New", font=FONT_HEADER, bg='white')
        label_new.grid(row=0, column=0, columnspan=16, rowspan=1, sticky='EW', padx=PAD_SMALL)

        self.button_new_experiment = BigButton(self, text="ADD EXPERIMENT", height=int(GUI_HEIGHT / 4),
                                               width=int(GUI_WIDTH / 7), command=lambda: self.add_experiment())
        self.button_new_experiment.grid(row=1, column=0, columnspan=16, rowspan=4, sticky='EW', padx=PAD_SMALL)

        self.button_new_dataset = BigButton(self, text="ADD DATASET", height=int(GUI_HEIGHT / 4),
                                            width=int(GUI_WIDTH / 7), command=lambda: self.add_dataset())
        self.button_new_dataset.grid(row=5, column=0, columnspan=16, rowspan=4, sticky='EW', padx=PAD_SMALL)

        label_loaded = tk.Label(self, text="Loaded", font=FONT_HEADER, bg='white')
        label_loaded.grid(row=0, column=16, columnspan=16, sticky='EW', padx=PAD_SMALL)

        self.listbox_loaded = tk.Listbox(self, font=FONT_LISTBOX)
        self.listbox_loaded.grid(row=1, column=16, columnspan=16, rowspan=8, sticky='NSEW', padx=PAD_SMALL)
        self.listbox_loaded.configure(justify="center")

        self.button_loaded_delete = NormalButton(self, text="Delete", command=lambda: self.delete_experiment(),
                                                 row=9, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.button_loaded_deselect = NormalButton(self, text="Deselect", command=lambda: self.deselect_experiment(),
                                                   row=9, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        label_queued = tk.Label(self, text="Queued", font=FONT_HEADER, bg='white')
        label_queued.grid(row=0, column=32, columnspan=16, sticky='EW', padx=PAD_SMALL)

        self.listbox_queued = tk.Listbox(self, font=FONT_LISTBOX)
        self.listbox_queued.grid(row=1, column=32, columnspan=16, rowspan=8, sticky='NSEW', padx=PAD_SMALL)
        self.listbox_queued.configure(justify="center")
        self.listbox_queued.bindtags((self.listbox_queued, self, "all"))

        self.button_run = NormalButton(self, text="Run", command=lambda: self.run(),
                                       row=9, column=40, columnspan=8, sticky='EW', padx=PAD_SMALL)

        label_progress_task = tk.Label(self, text="Task Progress", font=FONT_HEADER, bg='white')
        label_progress_task.grid(row=13, column=0, columnspan=8, rowspan=2, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_progress_task, TOOLTIP_MAIN_PROGRESS_TASK)

        label_progress_overall = tk.Label(self, text="Overall Progress", font=FONT_HEADER, bg='white')
        label_progress_overall.grid(row=15, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_progress_overall, TOOLTIP_MAIN_PROGRESS_OVERALL)

        self.label_progress_task_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                      row=13, column=8, columnspan=8, rowspan=2,
                                                      sticky="ew", font=FONT_LABEL)
        self.label_progress_overall_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                         row=15, column=8, columnspan=8, rowspan=2,
                                                         sticky="ew", font=FONT_LABEL)

        label_current_task = tk.Label(self, text="Current Task", font=FONT_HEADER, bg='white')
        label_current_task.grid(row=13, column=24, columnspan=8, rowspan=2, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_current_task, TOOLTIP_MAIN_CURRENT_TASK)

        self.label_current_task_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                     row=13, column=32, columnspan=16, rowspan=2,
                                                     sticky="ew", font=FONT_LABEL)

        label_time_done = tk.Label(self, text="Time Done", font=FONT_HEADER, bg='white')
        label_time_done.grid(row=15, column=24, columnspan=8, rowspan=2, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_time_done, TOOLTIP_MAIN_TIME_DONE)

        self.label_time_done_status = NormalLabel(self, text="Not yet started", bd=1, relief='sunken',
                                                  row=15, column=32, columnspan=16, rowspan=2,
                                                  sticky="ew", font=FONT_LABEL)

        # set progress updater to control created labels
        self.controller.progress_updater = ProgressUpdaterGUI(self, self.label_progress_task_status,
                                                              self.label_progress_overall_status,
                                                              self.label_current_task_status,
                                                              self.label_time_done_status)

    def add_experiment(self):
        """
        Add experiment function. Simply shows LoadPage and tells it no experiment to link to
        """
        self.controller.experiment_to_link_name = None
        self.controller.show_page(LoadPage)

    def add_dataset(self):
        """
        Add dataset function. Shows LoadPage and also links to an experiment
        """
        try:
            # try to get experiment from listbox. If fail, none selected
            selected = self.listbox_loaded.get(self.listbox_loaded.curselection())
            self.controller.experiment_to_link_name = selected[7:]
            self.controller.show_page(LoadPage)
        except:
            tk.messagebox.showerror("ERROR", "No experiment selected, please select one to link dataset to")

    def delete_experiment(self):
        """
        Deletes an experiment
        """
        try:
            # try to get experiment from listbox
            selected = self.listbox_loaded.get(self.listbox_loaded.curselection())
            name = selected.split(" ")[-1]
            # delete the correct experiment
            for index, experiment in enumerate(self.controller.experiments):
                if name == experiment.name:
                    del self.controller.experiments[index]
                    break
            # update MainPage
            self.update_page()
        except:
            return

    def deselect_experiment(self):
        """
        Deselects listbox
        """
        self.listbox_loaded.selection_clear(0, "end")

    def run_thread(self):
        """
        The function that is run within the thread when run is called
        """
        # analyze experiments
        for exp_index, experiment in enumerate(self.controller.experiments):
            self.controller.progress_updater.new_experiment(exp_index, experiment.settings['All Figures'],
                                                            len(experiment.rois))
            experiment.run()

        # close down thread
        self.close_down()

    def run(self):
        """
        Wrapper run function. Sets some things and starts up thread
        """
        if self.controller.thread_started is False:
            self.controller.thread_started = True
            self.controller.progress_updater.start(self.controller.experiments)
            # disable all buttons and start analysis
            self.disable()
            _thread.start_new_thread(self.run_thread, ())

    def disable(self):
        """
        Disables all buttons
        """
        self.button_run.updater(state='disabled')
        self.button_new_experiment.updater(state='disabled')
        self.button_new_dataset.updater(state='disabled')
        self.button_loaded_delete.updater(state='disabled')
        self.button_loaded_deselect.updater(state='disabled')

    def close_down(self):
        """
        Closes down run thread by enabling all buttons and clearing experiments
        """
        self.button_run.updater(command=lambda: self.run(), state='enabled')
        self.button_new_experiment.updater(state='enabled')
        self.button_new_dataset.updater(state='enabled')
        self.button_loaded_delete.updater(state='enabled', command=lambda: self.delete_experiment())
        self.button_loaded_deselect.updater(state='enabled', command=lambda: self.deselect_experiment())

        self.controller.experiments = []
        self.update_page()
        self.controller.thread_started = False

        # reset backend for further use in GUI
        mpl.use("TkAgg", force=True)
        from matplotlib import pyplot as plt
        tk.messagebox.showinfo("Done", "Processing is done!")

    def update_page(self, experiment=None):
        """
        Update page by clearing and filling list boxes
        :param experiment: experiment to sent to page if need be
        """
        self.listbox_loaded.delete(0, 'end')
        self.listbox_queued.delete(0, 'end')

        for index, experiment in enumerate(self.controller.experiments, 1):
            self.listbox_loaded.insert('end', "Exp {}: {}".format(index, experiment.name))
            for dataset in experiment.datasets:
                self.listbox_queued.insert('end', "Exp {}: {} ({})".format(index, dataset.name, dataset.type))

# %% Loading page


class LoadPage(BasePage):
    """
    Loading page. On this page, there are only two big buttons to select which type of dataset you want to load
    """
    def __init__(self, container, controller):
        super().__init__(container, controller)

        button1 = BigButton(self, text="TT", height=int(GUI_HEIGHT / 6),
                            width=int(GUI_WIDTH / 8),
                            command=lambda: self.load_nd2("TT"))
        button1.grid(row=10, column=24, columnspan=1, rowspan=1, padx=PAD_SMALL)

        button2 = BigButton(self, text="HSM", height=int(GUI_HEIGHT / 6),
                            width=int(GUI_WIDTH / 8),
                            command=lambda: self.load_nd2("HSM"))
        button2.grid(row=10, column=25, columnspan=1, rowspan=1, padx=PAD_SMALL)

        self.label_wait = NormalLabel(self, text="HSM frames are being merged. Progress 0%", row=11, column=24,
                                      columnspan=2, font=FONT_LABEL, padx=PAD_SMALL, sticky='EW')
        self.label_wait.grid_remove()

    def load_nd2(self, dataset_type):
        """
        Function to load an nd2
        :param dataset_type: Type is given by which button you click
        """
        filename = askopenfilename(filetypes=FILETYPES,
                                   title="Select nd2",
                                   initialdir=self.controller.dir_open)

        if len(filename) == 0:
            return

        # save directory
        self.controller.dir_open = '/'.join(filename[:-4].split("/")[:-1])

        if dataset_type == "HSM":
            if self.bad_hsm_size(filename):
                return
            # if datatype is HSM, show wait label
            self.label_wait.grid()

        if self.controller.experiment_to_link_name is None:
            # if no experiment to link to, new experiment
            experiment = Experiment(dataset_type, filename, self.controller.proceed_question, tk.messagebox.showerror,
                                    self.controller.progress_updater, self.controller.show_rois, label=self.label_wait)
            self.controller.experiments.append(experiment)
            # show ROIPage
            self.controller.show_page(ROIPage, experiment=experiment)
        else:
            # otherwise, link to experiment
            experiment_to_link = [experiment for experiment in self.controller.experiments if
                                  self.controller.experiment_to_link_name == experiment.name][0]
            if dataset_type == "TT":
                experiment_to_link.init_new_tt(filename)
                self.controller.show_page(TTPage, experiment=experiment_to_link)
            else:
                experiment_to_link.init_new_hsm(filename, label=self.label_wait)
                self.controller.show_page(HSMPage, experiment=experiment_to_link)

        # remove wait label
        if dataset_type == "HSM":
            self.label_wait.grid_remove()
            self.label_wait.updater(text="HSM frames are being merged. Progress 0%")

    def bad_hsm_size(self, filename):
        nd2 = ND2ReaderSelf(filename)
        if len(nd2) >  50:
            nd2.close()
            return not self.controller.proceed_question("Are you sure?",
                                                        "This HSM is over 50 frames, which is unusually long.")
        else:
            nd2.close()
            return False

# %% ROIPage


class ROIPage(BasePage):
    """
    Page on which you can find the ROIs in your experiment
    """
    def __init__(self, container, controller):
        """
        Sets all data and GUI
        :param container: container
        :param controller: controller
        """
        super().__init__(container, controller)

        # set data
        self.experiment = None
        self.default_settings = None
        self.saved_settings = None
        self.histogram_fig = None
        self.to_hist = None

        label_name = tk.Label(self, text="Name Experiment", font=FONT_SUBHEADER, bg='white')
        label_name.grid(row=0, column=0, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_name, TOOLTIP_ROI_NAME_EXPERIMENT)

        self.entry_name = EntryPlaceholder(self, "TBD", width=INPUT_BIG, small=False)
        self.entry_name.grid(row=0, column=8, columnspan=40, sticky='EW')

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=1, column=0, rowspan=1, columnspan=40, sticky='we')

        label_settings = tk.Label(self, text="Settings", font=FONT_SUBHEADER, bg='white')
        label_settings.grid(row=2, column=16, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_settings, TOOLTIP_ROI_SETTINGS)

        label_min_int = tk.Label(self, text="Minimum Intensity", font=FONT_LABEL, bg='white')
        label_min_int.grid(row=2, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_min_int, TOOLTIP_ROI_MIN_INTENSITY)
        var_min_int = tk.IntVar()
        entry_min_int = EntrySlider(self, row=3, column=0, columnspan=8, textvariable=var_min_int, width=INPUT_MEDIUM)

        def slider_min_int(_):
            var_min_int.set(self.slider_min_int.get())
        self.slider_min_int = NormalSlider(self, from_=0, to=1000, command=slider_min_int,
                                           row=4, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)
        entry_min_int.binder(self.slider_min_int, var_min_int)

        self.button_min_int_histogram = NormalButton(self, text="Graph", command=lambda: self.fun_histogram("min_int"),
                                                     row=3, column=16, columnspan=8, rowspan=2,
                                                     sticky='EW', padx=PAD_SMALL)
        self.button_min_int_histogram_select = NormalButton(self, text="Select min", state='disabled', rowspan=2,
                                                            row=3, column=8, columnspan=8, sticky='EW', padx=PAD_SMALL)
        self.button_max_int_histogram_select = NormalButton(self, text="Select max", state='disabled', rowspan=2,
                                                            row=3, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)

        label_max_int = tk.Label(self, text="Maximum Intensity", font=FONT_LABEL, bg='white')
        label_max_int.grid(row=2, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_max_int, TOOLTIP_ROI_MAX_INTENSITY)
        var_max_int = tk.IntVar()
        entry_max_int = EntrySlider(self, row=3, column=32, columnspan=8, textvariable=var_max_int, width=INPUT_MEDIUM)

        def slider_max_int(_):
            var_max_int.set(self.slider_max_int.get())
        self.slider_max_int = NormalSlider(self, from_=0, to=5000, command=slider_max_int,
                                           row=4, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)
        entry_max_int.binder(self.slider_max_int, var_max_int)

        label_min_sigma = tk.Label(self, text="Minimum Sigma", font=FONT_LABEL, bg='white')
        label_min_sigma.grid(row=5, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_min_sigma, TOOLTIP_ROI_MIN_SIGMA)
        var_min_sigma = tk.DoubleVar()
        entry_min_sigma = EntrySlider(self, row=6, column=0, columnspan=8,
                                      textvariable=var_min_sigma, width=INPUT_MEDIUM)

        def slider_min_sigma(_):
            var_min_sigma.set(self.slider_min_sigma.get())
        self.slider_min_sigma = NormalSlider(self, from_=0, to=5, resolution=0.01, command=slider_min_sigma,
                                             row=7, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)
        entry_min_sigma.binder(self.slider_min_sigma, var_min_sigma)

        self.button_min_sigma_histogram = NormalButton(self, text="Graph",
                                                       command=lambda: self.fun_histogram("min_sigma"), rowspan=2,
                                                       row=6, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)
        self.button_min_sigma_histogram_select = NormalButton(self, text="Select min", state='disabled',
                                                              row=6, column=8, columnspan=8, sticky='EW', rowspan=2,
                                                              padx=PAD_SMALL)
        self.button_max_sigma_histogram_select = NormalButton(self, text="Select max", state='disabled', rowspan=2,
                                                              row=6, column=24, columnspan=8, sticky='EW',
                                                              padx=PAD_SMALL)

        label_max_sigma = tk.Label(self, text="Maximum Sigma", font=FONT_LABEL, bg='white')
        label_max_sigma.grid(row=5, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_max_sigma, TOOLTIP_ROI_MAX_SIGMA)
        var_max_sigma = tk.DoubleVar()
        entry_max_sigma = EntrySlider(self, row=6, column=32, columnspan=8,
                                      textvariable=var_max_sigma, width=INPUT_MEDIUM)

        def slider_max_sigma(_):
            var_max_sigma.set(self.slider_max_sigma.get())
        self.slider_max_sigma = NormalSlider(self, from_=0, to=10, resolution=0.01, command=slider_max_sigma,
                                             row=7, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)
        entry_max_sigma.binder(self.slider_max_sigma, var_max_sigma)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=8, column=0, rowspan=1, columnspan=40, sticky='we')

        label_advanced_settings = tk.Label(self, text="Advanced settings", font=FONT_SUBHEADER, bg='white')
        label_advanced_settings.grid(row=9, column=0, columnspan=40, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_advanced_settings, TOOLTIP_ROI_ADVANCED_SETTINGS)

        label_min_corr = tk.Label(self, text="Minimum Correlation", font=FONT_LABEL, bg='white')
        label_min_corr.grid(row=10, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_min_corr, TOOLTIP_ROI_CORRELATION)
        var_min_corr = tk.DoubleVar()
        entry_min_corr = EntrySlider(self, row=10, column=8, columnspan=8,
                                     textvariable=var_min_corr, width=INPUT_MEDIUM)

        def slider_min_corr(_):
            var_min_corr.set(self.slider_min_corr.get())
        self.slider_min_corr = NormalSlider(self, from_=0, to=1, resolution=0.005, command=slider_min_corr,
                                            row=10, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)
        entry_min_corr.binder(self.slider_min_corr, var_min_corr)

        self.button_min_corr_histogram = NormalButton(self, text="Graph",
                                                      command=lambda: self.fun_histogram("corr_min"), rowspan=1,
                                                      row=10, column=24, columnspan=8, sticky='EW', padx=PAD_SMALL)
        self.button_min_corr_histogram_select = NormalButton(self, text="Graph select", state='disabled', rowspan=1,
                                                             row=10, column=32, columnspan=8, sticky='EW',
                                                             padx=PAD_SMALL)

        label_roi_side = tk.Label(self, text="Minimum spacing ROIs to edge FOV", bg='white', font=FONT_LABEL)
        label_roi_side.grid(row=12, column=0, columnspan=15, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_roi_side, TOOLTIP_ROI_EDGE_DIST)
        self.entry_roi_side = EntryPlaceholder(self, "11", width=INPUT_SMALL)
        self.entry_roi_side.grid(row=12, column=15, columnspan=5)

        label_inter_roi = tk.Label(self, text="Minimum spacing between ROIs", bg='white', font=FONT_LABEL)
        label_inter_roi.grid(row=12, column=20, columnspan=15, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_inter_roi, TOOLTIP_ROI_INTER_DIST)
        self.entry_inter_roi = EntryPlaceholder(self, "6", width=INPUT_SMALL)
        self.entry_inter_roi.grid(row=12, column=35, columnspan=5)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=13, column=0, rowspan=1, columnspan=40, sticky='we')

        label_other_settings = tk.Label(self, text="Other settings", font=FONT_SUBHEADER, bg='white')
        label_other_settings.grid(row=14, column=0, columnspan=40, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_other_settings, TOOLTIP_ROI_OTHER_SETTINGS)

        label_all_figures = tk.Label(self, text="Create individual figures", font=FONT_LABEL, bg='white')
        label_all_figures.grid(row=15, column=0, columnspan=15, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_all_figures, TOOLTIP_ROI_INDIVIDUAL_FIGURES)
        self.variable_all_figures = tk.StringVar(self, value=False)
        check_figures = ttk.Checkbutton(self, variable=self.variable_all_figures, onvalue=True, offvalue=False)
        check_figures.grid(row=15, column=15, columnspan=5, sticky='EW', padx=PAD_SMALL)

        label_filter_size = tk.Label(self, text="Background median filter size", bg='white', font=FONT_LABEL)
        label_filter_size.grid(row=15, column=20, columnspan=15, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_filter_size, TOOLTIP_ROI_BACKGROUND_FILTER)
        self.entry_filter_size = EntryPlaceholder(self, "9", width=INPUT_SMALL)
        self.entry_filter_size.grid(row=15, column=35, columnspan=5)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=16, column=0, rowspan=1, columnspan=40, sticky='we')

        button_find_rois = ttk.Button(self, text="Find ROIs", command=lambda: self.fit_rois())
        button_find_rois.grid(row=17, column=0, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.label_number_of_rois = NormalLabel(self, text="TBD", row=17, column=8, columnspan=8, font=FONT_LABEL,
                                                padx=PAD_SMALL, sticky='EW')

        button_restore = ttk.Button(self, text="Restore default",
                                    command=lambda: self.restore_default())
        button_restore.grid(row=17, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)
        self.button_restore_saved = NormalButton(self, text="Restore saved",
                                                 state='disabled',
                                                 command=lambda: self.restore_saved(),
                                                 row=17, column=24,
                                                 columnspan=8, sticky='EW', padx=PAD_SMALL)

        button_save = ttk.Button(self, text="Save", command=lambda: self.save_roi_settings())
        button_save.grid(row=17, column=32, columnspan=8, sticky='EW', padx=PAD_SMALL)

        self.figure = FigureFrame(self, height=GUI_WIDTH * 0.4, width=GUI_WIDTH * 0.4, dpi=DPI)
        self.figure.grid(row=1, column=40, columnspan=8, rowspan=17, sticky='EW', padx=PAD_SMALL)

        button_accept = ttk.Button(self, text="Accept & Continue", command=lambda: self.accept())
        button_accept.grid(row=18, column=44, columnspan=4, rowspan=2, sticky='EW', padx=PAD_SMALL)

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
        # find variable to make histogram of
        if variable == "min_int" or variable == "max_int":
            self.to_hist = self.experiment.roi_finder.main(return_int=True)
        elif variable == "corr_min":
            self.to_hist = self.experiment.roi_finder.main(return_corr=True)
        else:
            self.to_hist = self.experiment.roi_finder.main(return_sigmas=True)

        # make histogram
        self.histogram_fig = plt.figure(figsize=(6.4 * 1.2, 4.8 * 1.2))
        self.histogram_fig.canvas.mpl_connect('close_event', self.histogram_close)
        self.histogram_open(variable)
        self.histogram_make(variable)

    def histogram_make(self, variable):
        """
        Makes histograms of parameters.
        Parameters
        ----------
        variable : The parameter to make a histogram of
        Returns
        -------
        None, output figure
        """
        # add sub figure to histogram and create histogram
        fig_sub = self.histogram_fig.add_subplot(111)
        hist, bins, _ = fig_sub.hist(self.to_hist, bins='auto')

        # get sliders
        min_int = self.slider_min_int.get()
        max_int = self.slider_max_int.get()
        min_sigma = self.slider_min_sigma.get()
        max_sigma = self.slider_max_sigma.get()
        min_corr = self.slider_min_corr.get()

        # draw sliders
        if variable == "min_sigma" or variable == "max_sigma":
            fig_sub.set_title("Sigma. Use graph select to change threshold")
            fig_sub.axvline(x=min_sigma, color='red', linestyle='--')
            fig_sub.axvline(x=max_sigma, color='red', linestyle='--')
            fig_sub.set_xlabel('Sigma (pixels)')
            fig_sub.set_ylabel('Occurrence')
        else:

            # reset bins
            self.histogram_fig.clear()
            del fig_sub
            logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
            fig_sub = self.histogram_fig.add_subplot(111)
            fig_sub.hist(self.to_hist, bins=logbins)
            if variable == "min_int" or variable == "max_int":
                # intensity
                fig_sub.set_title("Intensity. Use graph select to change threshold")
                fig_sub.axvline(x=min_int, color='red', linestyle='--')
                fig_sub.axvline(x=max_int, color='red', linestyle='--')
                fig_sub.set_xlabel('Peak Intensity (counts)')
            else:
                # correlation
                fig_sub.set_title("Correlation values for ROIs. Use graph select to change threshold")
                fig_sub.axvline(x=min_corr, color='red', linestyle='--')
                fig_sub.set_xlabel('Correlation (arb. un.)')
            fig_sub.set_xscale('log')
            fig_sub.set_ylabel('Occurrence')

        # show figure
        self.histogram_fig.show()

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
        # find click on figure
        try:
            click = self.histogram_fig.ginput(1)
        except AttributeError:
            tk.messagebox.showerror("You cannot do that", "You cannot click on a figure without having a figure open")
            return

        # set slider to click
        if variable == "min_int":
            self.slider_min_int.updater(start=int(click[0][0]))
        elif variable == 'max_int':
            self.slider_max_int.updater(start=int(click[0][0]))
        elif variable == "min_sigma":
            self.slider_min_sigma.updater(start=click[0][0])
        elif variable == "max_sigma":
            self.slider_max_sigma.updater(start=click[0][0])
        else:
            self.slider_min_corr.updater(start=click[0][0])

        # clear histogram and make new figure
        self.histogram_fig.clear()
        self.histogram_fig.canvas.mpl_connect('close_event', self.histogram_close)
        self.histogram_make(variable)

    def histogram_open(self, variable):
        """
        Called when new histogram is opened.
        :param variable: Decides which buttons are enabled and disabled depending on which variables is plotted
        :return: None. Changes GUI
        """
        # enable the right button to select from graph with
        if variable == "min_int" or variable == "max_int":
            self.button_max_int_histogram_select.updater(state='enabled',
                                                         command=lambda: self.histogram_select("max_int"))
            self.button_min_int_histogram_select.updater(state='enabled',
                                                         command=lambda: self.histogram_select("min_int"))
        elif variable == "corr_min":
            self.button_min_corr_histogram_select.updater(state='enabled',
                                                          command=lambda: self.histogram_select("corr_min"))
        else:
            self.button_max_sigma_histogram_select.updater(state='enabled',
                                                           command=lambda: self.histogram_select("max_sigma"))
            self.button_min_sigma_histogram_select.updater(state='enabled',
                                                           command=lambda: self.histogram_select("min_sigma"))

        # disable opening other graphs
        self.button_min_corr_histogram.updater(state='disabled')
        self.button_min_int_histogram.updater(state='disabled')
        self.button_min_sigma_histogram.updater(state='disabled')

    def histogram_close(self, _):
        """
        Called when histogram is closed.
        :param _: Event is discarded
        :return: Changes GUI
        """
        # enable opening other graphs
        self.button_min_corr_histogram.updater(state='enabled', command=lambda: self.fun_histogram("corr_min"))
        self.button_min_int_histogram.updater(state='enabled', command=lambda: self.fun_histogram("min_int"))
        self.button_min_sigma_histogram.updater(state='enabled', command=lambda: self.fun_histogram("min_sigma"))

        # disable select from graphs
        self.button_max_int_histogram_select.updater(state='disabled')
        self.button_min_int_histogram_select.updater(state='disabled')
        self.button_min_corr_histogram_select.updater(state='disabled')
        self.button_max_sigma_histogram_select.updater(state='disabled')
        self.button_min_sigma_histogram_select.updater(state='disabled')

    def read_out_settings(self):
        """
        Function that reads out all the settings, and saves it to a dict which it returns
        Returns
        -------
        settings: a dictionary of all read-out settings
        """
        int_min = self.slider_min_int.get()
        int_max = self.slider_max_int.get()
        sigma_min = self.slider_min_sigma.get()
        sigma_max = self.slider_max_sigma.get()
        corr_min = self.slider_min_corr.get()
        roi_size = 7  # hard-coded as 7 for ROI finding
        all_figures = bool(int(self.variable_all_figures.get()))

        # try to read out inputs and check integer
        try:
            filter_size = int(self.entry_filter_size.get())
            roi_side = int(self.entry_roi_side.get())
            inter_roi = int(self.entry_inter_roi.get())
        except:
            tk.messagebox.showerror("ERROR", "Filter size, side spacing, and ROI spacing must all be integers")
            return {}, False

        settings = {'int_max': int_max, 'int_min': int_min,
                    'sigma_min': sigma_min, 'sigma_max': sigma_max,
                    'corr_min': corr_min, 'roi_size': roi_size, 'filter_size': filter_size,
                    'roi_side': roi_side, 'inter_roi': inter_roi, 'all_figures': all_figures}

        # return settings and success boolean
        return settings, True

    def fit_rois(self):
        """
        Reads out settings and correlates frame with old frame to find active ROIs
        :return: alters active ROIs within dataset
        """
        settings, success = self.read_out_settings()
        if success is False:
            return False

        # check settings
        if settings['roi_side'] < int((settings['roi_size'] - 1) / 2):
            tk.messagebox.showerror("ERROR", "Distance to size cannot be smaller than 1D ROI size")
            return False
        if settings['filter_size'] % 2 != 1:
            tk.messagebox.showerror("ERROR", "Filter size should be odd")
            return False

        # change settings and show new ROIs
        self.experiment.change_rois(settings)
        self.experiment.show_rois("Experiment", figure=self.figure, overwrite=True)
        self.label_number_of_rois.updater(text="{} ROIs found".format(len(self.experiment.rois)))

        return True

    def restore_default(self):
        """
        Restores default settings
        :return: Changes GUI
        """
        # find default settings if need be
        if self.default_settings is None:
            self.default_settings = self.experiment.roi_finder.get_settings()
        else:
            pass
        # change all sliders and such
        self.slider_min_int.updater(from_=0, to=self.default_settings['int_max'] / 4,
                                    start=self.default_settings['int_min'])
        self.slider_max_int.updater(from_=0, to=self.default_settings['int_max'],
                                    start=self.default_settings['int_max'])
        self.slider_min_sigma.updater(from_=0, to=self.default_settings['sigma_max'],
                                      start=self.default_settings['sigma_min'])
        self.slider_max_sigma.updater(from_=0, to=self.default_settings['sigma_max'],
                                      start=self.default_settings['sigma_max'])
        self.slider_min_corr.updater(from_=0, to=1, start=self.default_settings['corr_min'])

        self.variable_all_figures.set(False)

        self.entry_filter_size.updater()
        self.entry_roi_side.updater()
        self.entry_inter_roi.updater()

        self.update()
        self.fit_rois()

    def restore_saved(self):
        """
        Restores saved settings to sliders etc.
        Returns
        -------
        None, updates GUI
        """
        # gets saved settings and changes sliders and such
        settings = self.saved_settings

        self.slider_min_int.updater(from_=0, to=self.default_settings['int_max'] / 4,
                                    start=settings['int_min'])
        self.slider_max_int.updater(from_=0, to=self.default_settings['int_max'],
                                    start=settings['int_max'])
        self.slider_min_sigma.updater(from_=0, to=self.default_settings['sigma_max'],
                                      start=settings['sigma_min'])
        self.slider_max_sigma.updater(from_=0, to=self.default_settings['sigma_max'],
                                      start=settings['sigma_max'])
        self.slider_min_corr.updater(from_=0, to=1, start=settings['corr_min'])

        self.entry_filter_size.updater(settings['filter_size'])
        self.entry_roi_side.updater(settings['roi_side'])
        self.entry_inter_roi.updater(settings['inter_roi'])

        self.variable_all_figures.set(settings['all_figures'])

        self.experiment.change_rois(settings)
        self.experiment.show_rois("Experiment", figure=self.figure)
        self.label_number_of_rois.updater(text="{} ROIs found".format(len(self.experiment.rois)))
        self.update()

    def save_roi_settings(self):
        """
        Save settings for later re-use
        """
        if self.fit_rois() is False:
            return

        self.saved_settings, _ = self.read_out_settings()
        self.button_restore_saved.updater()

    def clear_page(self):
        """
        Clear the memory of the page
        """
        self.experiment = None
        self.default_settings = None
        self.saved_settings = None
        self.histogram_fig = None
        self.to_hist = None
        self.figure.fig.clear()

    def accept(self):
        """
        Accept ROI settings and move to analysis page
        """
        # check
        if self.controller.proceed_question("Are you sure?", "You cannot change settings later.") is False:
            return
        settings, success = self.read_out_settings()
        if success is False:
            return
        # change ROI settings to latest set
        self.experiment.change_rois(settings)

        # get name and all figures variable
        name = self.entry_name.get()
        settings_experiment = {'All Figures': settings['all_figures']}
        self.experiment.finalize_rois(name, settings_experiment)

        # show correct analysis page
        if self.experiment.created_by == "TT":
            self.controller.show_page(TTPage, experiment=self.experiment)
        else:
            self.controller.show_page(HSMPage, experiment=self.experiment)
        # empty memory
        self.clear_page()

    def update_page(self, experiment=None):
        """
        Update page by changing name and showing figure
        :param experiment: experiment to show
        """
        self.experiment = experiment

        self.entry_name.updater(placeholder=self.experiment.datasets[-1].name)  # take name for only dataset in exp
        self.figure.updater(self.experiment.frame_for_rois)

        self.restore_default()

# %% Analysis Page Template


class AnalysisPageTemplate(BasePage):
    """
    Base analysis page. Inherits from BasePage and adds the correlation area and some initialisations
    """
    def __init__(self, container, controller):
        """
        Initializer of AnalysisPageTemplate. Sets a load of buttons and and experiment member.
        :param container: container
        :param controller: controller
        """
        super().__init__(container, controller)

        # experiment working on
        self.experiment = None

        label_loaded_video = tk.Label(self, text="Loaded nd2 name:", font=FONT_SUBHEADER, bg='white', anchor='e')
        label_loaded_video.grid(row=0, column=0, columnspan=16, rowspan=1, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_loaded_video, TOOLTIP_ANALYSIS_NAME_LOADED)
        self.label_loaded_video_status = NormalLabel(self, text="XX", row=0, column=16, columnspan=32, rowspan=1,
                                                     sticky="w", font=FONT_LABEL)

        label_name = tk.Label(self, text="Name dataset:", font=FONT_SUBHEADER, bg='white', anchor='e')
        label_name.grid(row=1, column=0, columnspan=16, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_name, TOOLTIP_ANALYSIS_NAME_DATASET)
        self.entry_name = EntryPlaceholder(self, "TBD", small=False)
        self.entry_name.grid(row=1, column=16, columnspan=32, sticky='EW', padx=PAD_SMALL)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=2, column=0, rowspan=1, columnspan=16, sticky='we')

        label_precrop = tk.Label(self, text="ROI finding for dataset", font=FONT_SUBHEADER, bg='white')
        label_precrop.grid(row=3, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_precrop, TOOLTIP_ANALYSIS_PRECROP)

        label_x_min = tk.Label(self, text="minimum x-value", font=FONT_LABEL, bg='white')
        label_x_min.grid(row=4, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_x_min, TOOLTIP_ANALYSIS_MIN_X)
        self.entry_x_min = EntryPlaceholder(self, "Leave empty for start")
        self.entry_x_min.grid(row=5, column=0, columnspan=16, padx=PAD_SMALL)
        label_x_max = tk.Label(self, text="maximum x-value", font=FONT_LABEL, bg='white')
        label_x_max.grid(row=6, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_x_max, TOOLTIP_ANALYSIS_MAX_X)
        self.entry_x_max = EntryPlaceholder(self, "Leave empty for end")
        self.entry_x_max.grid(row=7, column=0, columnspan=16, padx=PAD_SMALL)

        label_y_min = tk.Label(self, text="minimum y-value", font=FONT_LABEL, bg='white')
        label_y_min.grid(row=8, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_y_min, TOOLTIP_ANALYSIS_MIN_Y)
        self.entry_y_min = EntryPlaceholder(self, "Leave empty for start")
        self.entry_y_min.grid(row=9, column=0, columnspan=16, padx=PAD_SMALL)
        label_y_max = tk.Label(self, text="maximum y-value", font=FONT_LABEL, bg='white')
        label_y_max.grid(row=10, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_y_max, TOOLTIP_ANALYSIS_MAX_Y)
        self.entry_y_max = EntryPlaceholder(self, "Leave empty for end")
        self.entry_y_max.grid(row=11, column=0, columnspan=16, padx=PAD_SMALL)

        button_find_rois = ttk.Button(self, text="Find ROIs", command=lambda: self.fit_rois())
        button_find_rois.grid(row=12, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)

        self.figure_dataset = FigureFrame(self, height=GUI_WIDTH * 0.35, width=GUI_WIDTH * 0.35, dpi=DPI)
        self.figure_dataset.grid(row=2, column=16, columnspan=16, rowspan=11, sticky='EW', padx=PAD_SMALL)

        self.figure_experiment = FigureFrame(self, height=GUI_WIDTH * 0.35, width=GUI_WIDTH * 0.35, dpi=DPI)
        self.figure_experiment.grid(row=2, column=32, columnspan=16, rowspan=11, sticky='EW', padx=PAD_SMALL)

        line = ttk.Separator(self, orient='horizontal')
        line.grid(row=13, column=0, rowspan=1, columnspan=48, sticky='we')

        self.button_add_to_queue = NormalButton(self, text="Add to queue", state='disabled',
                                                row=19, column=42, columnspan=6, rowspan=1, sticky='EW', padx=PAD_SMALL)

    @staticmethod
    def check_invalid_input(input_string, start):
        """
        Check whether or not given string input is invalid input
        :param input_string: Input string
        :param start: Whether or not start or end
        :return: Boolean if valid
        """
        def is_int(to_check):
            # check if int
            try:
                int(to_check)
                return True
            except:
                return False

        # if start, check if base value or valid int
        if start:
            if input_string == "Leave empty for start" or is_int(input_string):
                return False
        else:  # otherwise, check for other base value or still int
            if input_string == "Leave empty for end" or is_int(input_string):
                return False
        return True

    def fit_rois(self):
        """
        Correlate frames to find active ROIs
        :return: Changes active ROIs
        """
        # get values from GUI
        x_min = self.entry_x_min.get()
        x_max = self.entry_x_max.get()
        y_min = self.entry_y_min.get()
        y_max = self.entry_y_max.get()
        # check invalid
        if self.check_invalid_input(x_min, True) or self.check_invalid_input(y_max, False) or \
                self.check_invalid_input(y_min, True) or self.check_invalid_input(y_max, False):
            tk.messagebox.showerror("ERROR", "x min and max and y min and max must all be integers")
            return

        # set settings and sent to dataset
        settings_correlation = {'x_min': x_min, 'x_max': x_max, 'y_min': y_min, 'y_max': y_max}
        self.experiment.find_rois_dataset(settings_correlation)
        self.experiment.show_rois("Dataset", self.figure_dataset)
        self.button_add_to_queue.updater(command=lambda: self.add_to_queue())

    def clear_page(self):
        """
        Clears memory of page
        """
        self.experiment = None
        self.figure_experiment.fig.clear()
        self.figure_dataset.fig.clear()

    def add_to_queue(self):
        """
        Empty add_to_queue. Changed by inheritance
        """
        pass

    def update_page(self, experiment=None):
        """
        Update page by adding experiment data to page.
        :param experiment: experiment to analyze
        """
        self.experiment = experiment
        experiment.show_rois("Experiment", self.figure_experiment)
        experiment.show_rois("Dataset", self.figure_dataset)

        self.label_loaded_video_status.updater(text=self.experiment.datasets[-1].name)
        self.entry_name.updater(placeholder=self.experiment.datasets[-1].name)

        self.entry_x_min.updater()
        self.entry_y_min.updater()
        self.entry_x_max.updater()
        self.entry_y_max.updater()

# %% TTPage


class TTPage(AnalysisPageTemplate):
    """
    TT Page. Inherits of AnalysisPageTemplate
    """
    def __init__(self, container, controller):
        """
        Initializer of TTPage. Adds even more buttons specifically for TT analysis
        :param container: container
        :param controller: controller
        """
        super().__init__(container, controller)

        label_tt = tk.Label(self, text="TT settings", font=FONT_SUBHEADER, bg='white')
        label_tt.grid(row=14, column=0, columnspan=48, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_tt, TOOLTIP_TT_MAIN)

        label_method = tk.Label(self, text="Method", font=FONT_LABEL, bg='white')
        label_method.grid(row=15, column=0, columnspan=16, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_method, TOOLTIP_TT_METHOD)
        self.variable_method = tk.StringVar(self)
        drop_method = ttk.OptionMenu(self, self.variable_method, fit_options[1], *fit_options)
        drop_method.grid(row=16, column=0, columnspan=16, sticky="ew")

        label_rejection = tk.Label(self, text="Rejection", bg='white', font=FONT_LABEL)
        label_rejection.grid(row=15, column=16, columnspan=8, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_rejection, TOOLTIP_TT_REJECTION)
        self.variable_rejection = tk.BooleanVar(self, value=False)
        check_rejection = ttk.Checkbutton(self, variable=self.variable_rejection, onvalue=True, offvalue=False)
        check_rejection.grid(row=16, column=16, columnspan=8, padx=PAD_SMALL)

        label_cores = tk.Label(self, text="#cores", font=FONT_LABEL, bg='white')
        label_cores.grid(row=15, column=24, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_cores, TOOLTIP_TT_CORES)
        total_cores = mp.cpu_count()
        cores_options = [1, int(total_cores / 2), int(total_cores * 3 / 4), int(total_cores)]
        self.variable_cores = tk.IntVar(self)
        drop_cores = ttk.OptionMenu(self, self.variable_cores, cores_options[0], *cores_options)
        drop_cores.grid(row=16, column=24, columnspan=8, sticky='EW', padx=PAD_BIG)

        label_dimensions = tk.Label(self, text="output in pixels or nm", font=FONT_LABEL, bg='white')
        label_dimensions.grid(row=15, column=32, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_dimensions, TOOLTIP_TT_PIXELS_OR_NM)
        self.variable_dimensions = tk.StringVar(self)
        drop_dimension = ttk.OptionMenu(self, self.variable_dimensions, dimension_options[0], *dimension_options)
        drop_dimension.grid(row=16, column=32, columnspan=8, sticky='EW', padx=PAD_BIG)

        label_correlation_interval = tk.Label(self, text="Correlation interval", font=FONT_LABEL, bg='white')
        label_correlation_interval.grid(row=15, column=40, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_correlation_interval, TOOLTIP_TT_CORRELATION_INTERVAL)
        self.entry_correlation_interval = EntryPlaceholder(self, "Never")
        self.entry_correlation_interval.grid(row=16, column=40, rowspan=1, columnspan=8, padx=PAD_SMALL)

        label_used_roi_spacing = tk.Label(self, text="Used ROI spacing:", bg='white', font=FONT_LABEL)
        label_used_roi_spacing.grid(row=18, column=0, rowspan=1, columnspan=10, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_used_roi_spacing, TOOLTIP_TT_USED_ROI_SPACING)
        self.label_roi_spacing_status = NormalLabel(self, text="TBD", row=18, column=10, rowspan=1, columnspan=6,
                                                    sticky='EW', padx=PAD_SMALL, font=FONT_LABEL)

        label_roi_size = tk.Label(self, text="ROI size", bg='white', font=FONT_LABEL)
        label_roi_size.grid(row=19, column=0, columnspan=10, rowspan=1, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_roi_size, TOOLTIP_TT_ROI_SIZE)
        self.variable_roi_size = tk.StringVar(self)
        drop_roi_size = ttk.OptionMenu(self, self.variable_roi_size, roi_size_options[0], *roi_size_options)
        drop_roi_size.grid(row=19, column=10, columnspan=6, rowspan=1, sticky='EW', padx=PAD_SMALL)

        label_begin_frame = tk.Label(self, text="First frame number", font=FONT_LABEL, bg='white')
        label_begin_frame.grid(row=18, column=16, rowspan=1, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_begin_frame, TOOLTIP_TT_FIRST_FRAME)
        self.entry_begin_frame = EntryPlaceholder(self, "Leave empty for start")
        self.entry_begin_frame.grid(row=19, column=16, rowspan=1, columnspan=8, padx=PAD_SMALL)

        label_end_frame = tk.Label(self, text="Last frame number", font=FONT_LABEL, bg='white')
        label_end_frame.grid(row=18, column=24, rowspan=1, columnspan=8, sticky='EW', padx=PAD_BIG)
        create_tooltip(label_end_frame, TOOLTIP_TT_LAST_FRAME)
        self.entry_end_frame = EntryPlaceholder(self, "Leave empty for end")
        self.entry_end_frame.grid(row=19, column=24, rowspan=1, columnspan=8, padx=PAD_SMALL)

    def add_to_queue(self):
        """
        Add to queue specific for TT analysis
        :return: None. Edits objects
        """
        # get all inputs
        name = self.entry_name.get()
        method = self.variable_method.get()
        rejection_type = self.variable_rejection.get()
        n_processes = self.variable_cores.get()
        dimension = self.variable_dimensions.get()
        frame_begin = self.entry_begin_frame.get()
        frame_end = self.entry_end_frame.get()
        roi_size = int(self.variable_roi_size.get()[0])
        corr_int = self.entry_correlation_interval.get()

        # check validity inputs
        if self.check_invalid_input(frame_begin, True) or self.check_invalid_input(frame_end, False):
            tk.messagebox.showerror("ERROR", "Frame begin and frame end must be integers")
            return

        # make settings dict and set to input
        settings_runtime = {'method': method, 'rejection': rejection_type, '#cores': n_processes,
                            'roi_size': roi_size, "pixels_or_nm": dimension, 'name': name,
                            'frame_begin': frame_begin, 'frame_end': frame_end, 'correlation_interval': corr_int}

        if self.experiment.add_to_queue(settings_runtime) is False:
            return

        self.controller.show_page(MainPage)
        # clear memory
        self.clear_page()

    def update_page(self, experiment=None):
        """
        Update page by setting standard values and showing experiment
        :param experiment: experiment to show
        :return: changes GUI
        """
        super().update_page(experiment=experiment)
        self.label_roi_spacing_status.updater(text=self.experiment.roi_finder.get_settings()['inter_roi'])

        self.variable_method.set(fit_options[1])
        self.variable_rejection.set(True)
        self.variable_cores.set(1)
        self.variable_dimensions.set(dimension_options[0])
        self.variable_roi_size.set(roi_size_options[0])
        self.entry_begin_frame.updater()
        self.entry_end_frame.updater()
        self.entry_correlation_interval.updater()

        self.button_add_to_queue.updater(state='disabled')

# %% HSMPage


class HSMPage(AnalysisPageTemplate):
    """
    HSM Page. Inherits from AnalysisPageTemplate.
    """
    def __init__(self, container, controller):
        """
        Sets more buttons specifically for HSM analysis
        -------------------
        :param container: container
        :param controller: controller
        """
        super().__init__(container, controller)

        label_hsm = tk.Label(self, text="HSM settings", font=FONT_SUBHEADER, bg='white')
        label_hsm.grid(row=14, column=0, columnspan=48, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_hsm, TOOLTIP_HSM_MAIN)

        label_hsm_correct = tk.Label(self, text="Correction file:", font=FONT_LABEL, bg='white', anchor='e')
        label_hsm_correct.grid(row=15, column=0, columnspan=8, rowspan=2, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_hsm_correct, TOOLTIP_HSM_CORRECTION_FILE)
        path_hsm_correct = getcwd() + "/spectral_corrections"
        hsm_correct_options = listdir(path_hsm_correct)
        hsm_correct_options = [option[:-4] for option in hsm_correct_options]  # remove .mat in name
        self.variable_hsm_correct = tk.StringVar(self)
        drop_hsm_correct = ttk.OptionMenu(self, self.variable_hsm_correct, [], *hsm_correct_options)
        drop_hsm_correct.grid(row=15, column=8, columnspan=24, rowspan=2, sticky="ew")

        label_hsm_wavelength = tk.Label(self, text="Wavelengths:", font=FONT_LABEL, bg='white', anchor='e')
        label_hsm_wavelength.grid(row=18, column=0, columnspan=8, rowspan=2, sticky='EW', padx=PAD_SMALL)
        create_tooltip(label_hsm_wavelength, TOOLTIP_HSM_WAVELENGTHS)

        self.entry_wavelength = EntryPlaceholder(self,
                                                 "Use MATLAB-like array notation. "
                                                 "Example: [500:10:520, 532, 540:10:800]", width=INPUT_BIG)
        self.entry_wavelength.grid(row=18, column=8, columnspan=24, rowspan=2, sticky='EW')

    def update_page(self, experiment=None):
        """
        Set base values and show experiment
        :param experiment: experiment to show
        :return: Changes GUI
        """
        super().update_page(experiment=experiment)

        self.variable_hsm_correct.set("")
        self.entry_wavelength.updater()

        self.button_add_to_queue.updater(state='disabled')

    def add_to_queue(self):
        """
        Add HSM dataset to queue and clears memory
        """
        # get input values
        hsm_correction = self.variable_hsm_correct.get()
        wavelengths = self.entry_wavelength.get()
        name = self.entry_name.get()

        # check input values
        if wavelengths == "Use MATLAB-like array notation. Example: [500:10:520, 532, 540:10:800]":
            tk.messagebox.showerror("Check again", "Wavelengths not given.")
            return
        if hsm_correction == "":
            tk.messagebox.showerror("Check again", "Select a correction file.")
            return

        # make settings dictionary and add to queue
        settings_runtime_hsm = {'correction_file': hsm_correction, 'wavelengths': wavelengths, 'name': name}
        if self.experiment.add_to_queue(settings_runtime_hsm) is False:
            return

        self.controller.show_page(MainPage)
        # clear memory
        self.clear_page()

# %% START GUI


if __name__ == '__main__':
    mp.freeze_support()
    # setup logging
    logger, formatter = logging_setup()
    stream_handler = GUIStream()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.WARNING)
    logger.addHandler(stream_handler)
    logger.info("Started new run\n-----------------------------\n")
    divertor = DivertorErrorsGUI(logger)

    gui = PLASMON(proceed_question=proceed_question)

    tk.Tk.report_callback_exception = divertor.error
    plt.ioff()
    gui.mainloop()
