# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
MBx Python Data Analysis

v1.0: first version GUI: 13/07/2020
v1.1: Adding ROI size, wavelength to GUI: 13/07/2020
v1.2: Small additions, 7x7 ROI
v1.3: Disabling buttons, removing from grid
v1.4: Save and change wavelength per dataset

"""

## GENERAL IMPORTS
import os # to get standard usage
from sys import exit
import time # for timekeeping
from win32api import GetSystemMetrics # Get sys info



## Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("TkAgg") #set back end to TK

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import style
import matplotlib.patches as patches

## GUI
import tkinter as tk # for GUI
from tkinter import ttk # GUI styling
from tkinter.filedialog import askopenfilenames # for popup that asks to select .nd2's

## ND2 related
from pims import ND2_Reader # reader of ND2 files
#from pims import Frame # for converting ND2 generator to one numpy array

## Own code
import _code.roi_finding as roi_finding
import _code.fitters as fitting
import _code.tools as tools

#%% Inputs
ROI_SIZE = 9 
THRESHOLD = 5 # X*sigma

#%% Initializations

FILETYPES = [("ND2", ".nd2")]

LARGE_FONT = ("Verdana", 12)
#style.use("ggplot")

fit_options = ["Phasor with intensity", "Phasor without intensity",
                   "Gaussian with background", "Gaussian without background"]

roi_size_options = ["7x7", "9x9"]

ttk_style = ttk.Style()
#ttk_style.configure('my.TButton', font=('Verdana', 12))
ttk_style.configure('my.TButton', font=('Verdana', 1000))

width = GetSystemMetrics(0)
height = GetSystemMetrics(1)
GUI_WIDTH = int(width/2)
GUI_HEIGHT = int(height/2)
GUI_WIDTH_START = int(width/4)
GUI_HEIGHT_START = int(height/4)
DPI = 100
#%% Container

class mbx_python(tk.Tk):
    
    def __init__(self, *args, **kwargs):
    
        tk.Tk.__init__(self, *args, **kwargs)
        
        #tk.Tk.iconbitmap(self, default = "")
        tk.Tk.wm_title(self, "MBx Python")
        container = tk.Frame(self)
        
        container.pack(side="top", fill="both", expand = True)
        
        container.grid_rowconfigure(0, weight = 1)
        container.grid_columnconfigure(0, weight = 1)
        
        self.frames = {}
        
        frame_tuple = (load_page, fitting_page)
        
        for to_load_frame in frame_tuple:
            frame = to_load_frame(container, self)
            self.frames[to_load_frame] = frame
            frame.grid(row = 0, column = 0, sticky = "nsew")
        
        self.show_load_frame(load_page)
        
    def show_load_frame(self, cont):
        
        frame = self.frames[cont]
        frame.tkraise()
    
    def show_frame(self, cont, filenames):
        
        frame = self.frames[cont]
        frame.tkraise()
        frame.update_init(filenames)

#%% Quit
    
def quit_program():
    global gui
    gui.destroy()
    exit(0)
    
#%% Big button
        
class my_button(ttk.Frame):
    def __init__(self, parent, height=None, width=None, text="", command=None, style=None, state = 'enabled'):
        ttk.Frame.__init__(self, parent, height=height, width=width, style='my.TButton')

        self.pack_propagate(0)
        self._btn = ttk.Button(self, text=text, command=command, style=style, state=state)
        self._btn.pack(fill=tk.BOTH, expand=1)
        
#%% Entry with placeholder
        
class entry_placeholder(tk.Entry):
    def __init__(self, master=None, placeholder="PLACEHOLDER", color='grey'):
        super().__init__(master)

        self.placeholder = placeholder
        self.placeholder_color = color
        self.default_fg_color = self['fg']

        self.bind("<FocusIn>", self.foc_in)
        self.bind("<FocusOut>", self.foc_out)

        self.put_placeholder()

    def put_placeholder(self):
        self.insert(0, self.placeholder)
        self['fg'] = self.placeholder_color

    def foc_in(self, *args):
        if self['fg'] == self.placeholder_color:
            self.delete('0', 'end')
            self['fg'] = self.default_fg_color

    def foc_out(self, *args):
        if not self.get():
            self.put_placeholder()
            
#%% Loading page

class load_page(tk.Frame):
    
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)
        
        button1 = my_button(self, text = "LOAD", height = int(GUI_HEIGHT/4), 
                             width = int(GUI_WIDTH/4),# style= 'my.TButton',
                             command= lambda: self.load_nd2(controller))
        button1.grid(row = 0, column = 0)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
    def load_nd2(self, controller):
        
        global filenames
        tk.Tk().withdraw()
        filenames = askopenfilenames(filetypes = FILETYPES, 
                                     title = "Select file", 
                                     initialdir =os.getcwd())
        
        if len(filenames) == 0:
            return
        
        controller.show_frame(fitting_page, filenames)
        
#%% Fitting page, initial update (after loading)
        
class fitting_page(tk.Frame):
    
    def update_init(self, filenames):
        
        self.filenames = filenames
        self.roi_locations = {}
        self.dataset_index = 0
        
        self.dataset_roi_status.grid_forget()
        self.dataset_roi_status = tk.Label(self, text = "Dataset " + str(self.dataset_index + 1) + " of " + str(len(filenames)))
        self.dataset_roi_status.grid(row = 9, column = 1, columnspan = 4) 
        
        self.roi_status.grid_forget()
        self.roi_status = tk.Label(self, text = "0 of " + str(len(filenames)) + " have settings")
        self.roi_status.grid(row = 21, column = 0, columnspan = 6)
        
        self.ND2 = ND2_Reader(filenames[self.dataset_index])
        self.frames = self.ND2
        self.metadata = self.ND2.metadata
        
        self.fig.clear()        
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(self.frames[0], extent=[0,self.frames[0].shape[1],self.frames[0].shape[0],0],
                   aspect='auto')
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 16, sticky = 'E')
        
        self.roi_finder = roi_finding.roi_finder(9, self.frames[0])
                
        self.min_int_slider.grid_forget()
        self.min_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max/4, 
                                  orient = 'horizontal')
        self.min_int_slider.set(self.roi_finder.intensity_min)
        self.min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        self.max_int_slider.grid_forget()
        self.max_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max, orient = 'horizontal')
        self.max_int_slider.set(self.roi_finder.intensity_max)
        self.max_int_slider.grid(row = 3, column = 0, columnspan = 3)
                
        self.min_sigma_slider.grid_forget()
        self.min_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.min_sigma_slider.set(self.roi_finder.sigma_min)
        self.min_sigma_slider.grid(row = 1, column = 3, columnspan = 3)
                
        self.max_sigma_slider.grid_forget()
        self.max_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.max_sigma_slider.set(self.roi_finder.sigma_max)
        self.max_sigma_slider.grid(row = 3, column = 3, columnspan = 3)
        
        if self.dataset_index == len(filenames)-1:
            self.button_right.grid_forget()
            self.button_right = ttk.Button(self, text = ">>", state = 'disabled')
            self.button_right.grid(row = 9, column = 5)
        else:
            self.button_right.grid_forget()
            self.button_right = ttk.Button(self, text = ">>", command = lambda: self.next_dataset())
            self.button_right.grid(row = 9, column = 5)
            
        if self.dataset_index == 0:
            self.button_left.grid_forget()
            self.button_left = ttk.Button(self, text = "<<", state = 'disabled')
            self.button_left.grid(row = 9, column = 0)
        else:
            self.button_left.grid_forget()
            self.button_left = ttk.Button(self, text = "<<", command = lambda: self.previous_dataset())
            self.button_left.grid(row = 9, column = 0)
                        
    #%% Fitting page, fit ROIs
        
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
        
        self.roi_finder.change_settings(intensity_min = min_int, 
                        intensity_max = max_int,
                        sigma_min = min_sigma, sigma_max = max_sigma)
        
        fitter = fitting.scipy_last_fit_guess(self.metadata, ROI_SIZE,
                                                  wavelength, THRESHOLD, 
                                                  "ScipyLastFitGuess", 5)
        
        self.temp_roi_locations = self.roi_finder.main(self.frames[0], fitter)
        
        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)      
        
        fig_sub.imshow(self.frames[0], extent=[0,self.frames[0].shape[1],self.frames[0].shape[0],0],
                   aspect='auto')
        
        roi_locations_temp = self.temp_roi_locations - self.roi_finder.roi_size_1d
        roi_size = self.roi_finder.roi_size
        
        for roi in roi_locations_temp:
            rect = patches.Rectangle((roi[1], roi[0]), roi_size, roi_size,
                                 linewidth=0.5, edgecolor='r', facecolor='none')
            fig_sub.add_patch(rect)
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 16, sticky = 'E')
        
    #%% Fitting page, save ROI setings
        
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
        
        settings = {}
        settings['max_int'] = max_int
        settings['min_int'] = min_int
        settings['min_sigma'] = min_sigma
        settings['max_sigma'] = max_sigma
        settings['wavelength'] = wavelength
        
        self.saved_settings[self.dataset_index] = settings
        
        self.roi_status.grid_forget()
        self.roi_status = tk.Label(self, text = str(len(self.roi_locations)) + 
                              " of " + str(len(self.filenames)) + " have settings")
        self.roi_status.grid(row = 21, column = 0, columnspan = 6)
        
        self.button_restore_saved.grid_forget()
        self.button_restore_saved = ttk.Button(self, text = "Restore saved",
                                 command = lambda: self.restore_saved())
        self.button_restore_saved.grid(row = 12, column = 3, columnspan = 2)
        
        if len(self.roi_locations) == len(self.filenames):
            self.button_fit.grid_forget()
            self.button_fit = my_button(self, text = "FIT", height = int(GUI_HEIGHT/8), 
                                 width = int(GUI_WIDTH/8), 
                                 command = lambda: self.start_fitting())#, style= 'my.TButton')
            self.button_fit.grid(row = 24, column = 6, columnspan = 3, rowspan = 5)
        
    #%% Fitting page, restore default settings
    
    def restore_default(self):
        
        self.roi_finder = roi_finding.roi_finder(9, self.frames[0])
                
        self.min_int_slider.grid_forget()
        self.min_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max/4, 
                                  orient = 'horizontal')
        self.min_int_slider.set(self.roi_finder.intensity_min)
        self.min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        self.max_int_slider.grid_forget()
        self.max_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max, orient = 'horizontal')
        self.max_int_slider.set(self.roi_finder.intensity_max)
        self.max_int_slider.grid(row = 3, column = 0, columnspan = 3)
                
        self.min_sigma_slider.grid_forget()
        self.min_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max,
                                         orient = 'horizontal', resolution = 0.01)
        self.min_sigma_slider.set(self.roi_finder.sigma_min)
        self.min_sigma_slider.grid(row = 1, column = 3, columnspan = 3)
          
        self.max_sigma_slider.grid_forget()          
        self.max_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.max_sigma_slider.set(self.roi_finder.sigma_max)
        self.max_sigma_slider.grid(row = 3, column = 3, columnspan = 3)
        
    #%% Fitting page, switch between datasets
        
    def next_dataset(self):
            
        self.dataset_index +=1
        
        self.dataset_roi_status.grid_forget()
        self.dataset_roi_status = tk.Label(self, text = "Dataset " + str(self.dataset_index + 1) + " of " + str(len(self.filenames)))
        self.dataset_roi_status.grid(row = 9, column = 1, columnspan = 4) 
        
        self.ND2 = ND2_Reader(self.filenames[self.dataset_index])
        self.frames = self.ND2
        self.metadata = self.ND2.metadata
        
        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(self.frames[0], extent=[0,self.frames[0].shape[1],self.frames[0].shape[0],0],
                   aspect='auto')
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 16, sticky = 'E')
        
        self.roi_finder = roi_finding.roi_finder(9, self.frames[0])
                
        self.min_int_slider.grid_forget()
        self.min_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max/4, 
                                  orient = 'horizontal')
        self.min_int_slider.set(self.roi_finder.intensity_min)
        self.min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        self.max_int_slider.grid_forget()        
        self.max_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max, orient = 'horizontal')
        self.max_int_slider.set(self.roi_finder.intensity_max)
        self.max_int_slider.grid(row = 3, column = 0, columnspan = 3)
             
        self.min_sigma_slider.grid_forget()
        self.min_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.min_sigma_slider.set(self.roi_finder.sigma_min)
        self.min_sigma_slider.grid(row = 1, column = 3, columnspan = 3)
                    
        self.max_sigma_slider.grid_forget()
        self.max_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.max_sigma_slider.set(self.roi_finder.sigma_max)
        self.max_sigma_slider.grid(row = 3, column = 3, columnspan = 3)
        
        if self.dataset_index == len(self.filenames)-1:
            self.button_right.grid_forget()
            self.button_right = ttk.Button(self, text = ">>", state = 'disabled')
            self.button_right.grid(row = 9, column = 5)
        else:
            self.button_right.grid_forget()
            self.button_right = ttk.Button(self, text = ">>", command = lambda: self.next_dataset())
            self.button_right.grid(row = 9, column = 5)
            
        if self.dataset_index == 0:
            self.button_left.grid_forget()
            self.button_left = ttk.Button(self, text = "<<", state = 'disabled')
            self.button_left.grid(row = 9, column = 0)
        else:
            self.button_left.grid_forget()
            self.button_left = ttk.Button(self, text = "<<", command = lambda: self.previous_dataset())
            self.button_left.grid(row = 9, column = 0)
            
        if self.dataset_index in self.saved_settings:
            self.button_restore_saved.grid_forget()
            self.button_restore_saved = ttk.Button(self, text = "Restore saved", 
                                 command = lambda: self.restore_saved())
            self.button_restore_saved.grid(row = 12, column = 3, columnspan = 2)
            
            wavelength = self.saved_settings[self.dataset_index]['wavelength']
            self.wavelength_input = entry_placeholder(self, "wavelength in nm")
            self.wavelength_input.grid(row = 12, column = 6, columnspan = 3)
            self.wavelength_input.foc_in()
            self.wavelength_input.insert(0, str(wavelength))
        else:
            self.button_restore_saved.grid_forget()
            self.button_restore_saved = ttk.Button(self, text = "Restore saved", 
                                                   state = 'disabled',
                                 command = lambda: self.restore_saved())
            self.button_restore_saved.grid(row = 12, column = 3, columnspan = 2)   
            
            self.wavelength_input = entry_placeholder(self, "wavelength in nm")
            self.wavelength_input.grid(row = 12, column = 6, columnspan = 3)
        
    def previous_dataset(self):
        
        self.dataset_index -=1
        
        self.dataset_roi_status.grid_forget()
        self.dataset_roi_status = tk.Label(self, text = "Dataset " + str(self.dataset_index + 1) + " of " + str(len(filenames)))
        self.dataset_roi_status.grid(row = 9, column = 1, columnspan = 4) 
        
        self.ND2 = ND2_Reader(filenames[self.dataset_index])
        self.frames = self.ND2
        self.metadata = self.ND2.metadata
        
        self.fig.clear()
        fig_sub = self.fig.add_subplot(111)
        fig_sub.imshow(self.frames[0], extent=[0,self.frames[0].shape[1],self.frames[0].shape[0],0],
                   aspect='auto')
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 16, sticky = 'E')
        
        self.roi_finder = roi_finding.roi_finder(9, self.frames[0])
                
        self.min_int_slider.grid_forget()
        self.min_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max/4, 
                                  orient = 'horizontal')
        self.min_int_slider.set(self.roi_finder.intensity_min)
        self.min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        self.max_int_slider.grid_forget()        
        self.max_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max, orient = 'horizontal')
        self.max_int_slider.set(self.roi_finder.intensity_max)
        self.max_int_slider.grid(row = 3, column = 0, columnspan = 3)
             
        self.min_sigma_slider.grid_forget()
        self.min_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.min_sigma_slider.set(self.roi_finder.sigma_min)
        self.min_sigma_slider.grid(row = 1, column = 3, columnspan = 3)
                    
        self.max_sigma_slider.grid_forget()
        self.max_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.max_sigma_slider.set(self.roi_finder.sigma_max)
        self.max_sigma_slider.grid(row = 3, column = 3, columnspan = 3)
        
        if self.dataset_index == len(self.filenames)-1:
            self.button_right.grid_forget()
            self.button_right = ttk.Button(self, text = ">>", state = 'disabled')
            self.button_right.grid(row = 9, column = 5)
        else:
            self.button_right.grid_forget()
            self.button_right = ttk.Button(self, text = ">>", command = lambda: self.next_dataset())
            self.button_right.grid(row = 9, column = 5)
        
        if self.dataset_index == 0:
            self.button_left.grid_forget()
            self.button_left = ttk.Button(self, text = "<<", state = 'disabled')
            self.button_left.grid(row = 9, column = 0)
        else:
            self.button_left.grid_forget()
            self.button_left = ttk.Button(self, text = "<<", command = lambda: self.previous_dataset())
            self.button_left.grid(row = 9, column = 0)
            
        if self.dataset_index in self.saved_settings:
            self.button_restore_saved.grid_forget()
            self.button_restore_saved = ttk.Button(self, text = "Restore saved", 
                                 command = lambda: self.restore_saved())
            self.button_restore_saved.grid(row = 12, column = 3, columnspan = 2)
            
            wavelength = self.saved_settings[self.dataset_index]['wavelength']
            self.wavelength_input = entry_placeholder(self, "wavelength in nm")
            self.wavelength_input.grid(row = 12, column = 6, columnspan = 3)
            self.wavelength_input.foc_in()
            self.wavelength_input.insert(0, str(wavelength))
        else:
            self.button_restore_saved.grid_forget()
            self.button_restore_saved = ttk.Button(self, text = "Restore saved", 
                                                   state = 'disabled',
                                 command = lambda: self.restore_saved())
            self.button_restore_saved.grid(row = 12, column = 3, columnspan = 2)
            
            self.wavelength_input = entry_placeholder(self, "wavelength in nm")
            self.wavelength_input.grid(row = 12, column = 6, columnspan = 3)
            
    #%% Fitting page, start fitting
            
    def start_fitting(self):
        
        if len(self.roi_locations) != len(self.filenames):
            tk.messagebox.showerror("ERROR", "Not all datasets have settings yet, cannot start")
            return
            
        check = tk.messagebox.askokcancel("Are you sure?", "Fitting may take a while. Are you sure everything is set up correctly?" )
        if check == False:
            return
        
        basedir = os.getcwd()
        start_frame = self.frame_begin_input.get()
        end_frame = self.frame_end_input.get()
        method = self.method_var.get()
        roi_size = int(self.roi_var.get()[0])
        
        self.start_time = time.time()
            
        for dataset_index, filename in enumerate(self.filenames):
            self.ND2 = ND2_Reader(filenames[self.dataset_index])
            self.frames = self.ND2
            self.metadata = self.ND2.metadata
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
                fitter = fitting.scipy_last_fit_guess_background(self.metadata, roi_size, wavelength, min_intensity, method, 6)
                
            if start_frame == "Leave empty for start" and end_frame == "Leave empty for end":
                to_fit = self.frames
            elif start_frame == "Leave empty for start" and end_frame != "Leave empty for end":
                to_fit = self.frames[:int(end_frame)]
            elif start_frame != "Leave empty for start" and end_frame == "Leave empty for end":
                to_fit = self.frames[int(start_frame):]
            elif start_frame != "Leave empty for start" and end_frame != "Leave empty for end":
                to_fit = self.frames[int(start_frame):int(end_frame)]
            
            results = fitter.main(to_fit, self.metadata, roi_locations, 
                                  gui = self) 
            
            ## create folder for output
            
            directory = filename.split(".")[0].split("/")[-1]
            path = os.path.join(basedir, directory)
            try:
                os.mkdir(path)
            except:
                pass
            
            ## Save
            
            metadata_filtered = {k:v for k, v in self.metadata.items() if v is not None}
            del metadata_filtered['time_start']
            del metadata_filtered['time_start_utc']
    
            ## ROI_locations dict
            ROI_locations_dict = dict(zip(['x', 'y'], roi_locations.T))
    
            ## Localization dict
            results_dict = {'Localizations': results}
            
        tk.messagebox.showinfo("Done!")

        tools.save_to_csv_mat('metadata', metadata_filtered, directory)
        tools.save_to_csv_mat('ROI_locations', ROI_locations_dict, directory)
        tools.save_to_csv_mat('Localizations', results_dict, directory)
        
    #%% Fitting page, update the status
        
    def update_status(self, progress, comparator):
        
        num_files = len(self.filenames)
        base_progress = self.dataset_index / num_files * 100
        
        file_progress = progress/comparator*100/num_files
        progress = base_progress + file_progress
        
        current_time = time.time()
        time_taken = current_time - self.start_time
        time_done_estimate = time_taken * 100 / progress + self.start_time
        tr = time.localtime(time_done_estimate)
        
        progress_text = str(round(progress, 1)) + "% done"
        time_text = str(tr[3])+":"+str(tr[4])+":"+str(tr[5])+" "+str(tr[2])+"/"+str(tr[1])
        
        
        self.progress_status_label.grid_forget()      
        self.progress_status_label = tk.Label(self, text = progress_text, 
                                         bd = 1, relief = 'sunken')
        self.progress_status_label.grid(row = 25, column = 10)
        
        self.time_status_label.grid_forget()      
        self.time_status_label = tk.Label(self, text = time_text, 
                                         bd = 1, relief = 'sunken')
        self.time_status_label.grid(row = 28, column = 10)
        self.update()
        
    #%% Fitting page, return to load page
        
    def load_new(self, controller):
        
        controller.show_load_frame(load_page)
        
    #%% Fitting page, restore saved settings
    
    def restore_saved(self):
        
        settings = self.saved_settings[self.dataset_index]
        
        self.min_int_slider.grid_forget()       
        self.min_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max/4, 
                                  orient = 'horizontal')
        self.min_int_slider.set(settings['min_int'])
        self.min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        self.max_int_slider.grid_forget()  
        self.max_int_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.intensity_max, orient = 'horizontal')
        self.max_int_slider.set(settings['max_int'])
        self.max_int_slider.grid(row = 3, column = 0, columnspan = 3)
                  
        self.min_sigma_slider.grid_forget()  
        self.min_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max,
                                         orient = 'horizontal', resolution = 0.01)
        self.min_sigma_slider.set(settings['min_sigma'])
        self.min_sigma_slider.grid(row = 1, column = 3, columnspan = 3)
              
        self.max_sigma_slider.grid_forget()  
        self.max_sigma_slider = tk.Scale(self, from_ = 0, to = self.roi_finder.sigma_max, 
                                         orient = 'horizontal', resolution = 0.01)
        self.max_sigma_slider.set(settings['max_sigma'])
        self.max_sigma_slider.grid(row = 3, column = 3, columnspan = 3)
        
        wavelength = settings['wavelength']
        self.wavelength_input = entry_placeholder(self, "wavelength in nm")
        self.wavelength_input.grid(row = 12, column = 6, columnspan = 3)
        self.wavelength_input.foc_in()
        self.wavelength_input.insert(0, str(wavelength))
        
    #%% Fitting page, initial decleration
        
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)
        self.saved_settings = {}
        
        min_int_label = tk.Label(self, text = "Minimum Intensity", font = LARGE_FONT)
        min_int_label.grid(row = 0, column = 0, columnspan = 3)
        
        self.min_int_slider = tk.Scale(self, from_ = 0, to = 1000, orient = 'horizontal')
        self.min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        max_int_label = tk.Label(self, text = "Maximum Intensity", font = LARGE_FONT)
        max_int_label.grid(row = 2, column = 0, columnspan = 3)
        
        self.max_int_slider = tk.Scale(self, from_ = 0, to = 5000, orient = 'horizontal')
        self.max_int_slider.grid(row = 3, column = 0, columnspan = 3)
        
        min_sigma_label = tk.Label(self, text = "Minimum Sigma", font = LARGE_FONT)
        min_sigma_label.grid(row = 0, column = 3, columnspan = 3)
        
        self.min_sigma_slider = tk.Scale(self, from_ = 0, to = 5, 
                                         orient = 'horizontal', resolution = 0.01)
        self.min_sigma_slider.grid(row = 1, column = 3, columnspan = 3)
        
        max_sigma_label = tk.Label(self, text = "Maximum Sigma", font = LARGE_FONT)
        max_sigma_label.grid(row = 2, column = 3, columnspan = 3)
        
        self.max_sigma_slider = tk.Scale(self, from_ = 0, to = 10, 
                                         orient = 'horizontal', resolution = 0.01)
        self.max_sigma_slider.grid(row = 3, column = 3, columnspan = 3)
        
        self.button_left = ttk.Button(self, text = "<<")
        self.button_left.grid(row = 9, column = 0)
        
        self.dataset_roi_status = tk.Label(self, text = "TBD")
        self.dataset_roi_status.grid(row = 9, column = 1, columnspan = 4)
        
        self.button_right = ttk.Button(self, text = ">>")
        self.button_right.grid(row = 9, column = 5)
        
        button_fit = ttk.Button(self, text = "Fit", command= lambda: self.fit_rois())
        button_fit.grid(row = 12, column = 0, columnspan = 1)
        
        button_restore = ttk.Button(self, text = "Restore default", 
                                 command = lambda: self.restore_default())
        button_restore.grid(row = 12, column = 1, columnspan = 2)
        
        self.button_restore_saved = ttk.Button(self, text = "Restore saved", 
                                    state = 'disabled',
                                 command = lambda: self.restore_saved())
        self.button_restore_saved.grid(row = 12, column = 3, columnspan = 2)
        
        button_save = ttk.Button(self, text = "Save", 
                                 command = lambda: self.save_roi_settings())
        button_save.grid(row = 12, column = 5, columnspan = 1)
                
        self.fig = Figure(figsize = (GUI_HEIGHT/DPI*0.7,GUI_HEIGHT/DPI*0.7), dpi = DPI)
        
        canvas = FigureCanvasTkAgg(self.fig, self)
        canvas.draw()
        canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 14, sticky = 'E')
        
        #toolbar = NavigationToolbar2Tk(canvas, self)
        #toolbar.update()
        #canvas._tkcanvas.grid(row = 15, column = 10)
        
        wavelength_label = tk.Label(self, text = "Laser Wavelength", font = LARGE_FONT)
        wavelength_label.grid(row = 10, column = 6, columnspan = 3)
        
        self.wavelength_input = entry_placeholder(self, "wavelength in nm")
        self.wavelength_input.grid(row = 12, column = 6, columnspan = 3)
        
        line = ttk.Separator(self, orient='horizontal')
        line.grid(column=0, row=18, rowspan = 2,  columnspan = 10, sticky='we')
        
        self.roi_status = tk.Label(self, text = "TBD")
        self.roi_status.grid(row = 21, column = 0, columnspan = 6)
        
        method_label = tk.Label(self, text = "Method", font = LARGE_FONT)
        method_label.grid(column = 0, row = 23, columnspan = 3)
        
        self.method_var = tk.StringVar(self)
        self.method_var.set(fit_options[2])
        
        method_drop = tk.OptionMenu(self, self.method_var, *fit_options)
        method_drop.grid(column = 0, row = 24, columnspan = 3)
        
        roi_size_label = tk.Label(self, text = "ROI size", font = LARGE_FONT)
        roi_size_label.grid(column = 3, row = 23, columnspan = 3)
        
        self.roi_var = tk.StringVar(self)
        self.roi_var.set(roi_size_options[0])
        
        roi_drop = tk.OptionMenu(self, self.roi_var, *roi_size_options)
        roi_drop.grid(column = 3, row = 24, columnspan = 3)
        
        frame_begin_label = tk.Label(self, text = "Begin frame", font = LARGE_FONT)
        frame_begin_label.grid(row = 27, column = 0, columnspan = 3)
        
        self.frame_begin_input = entry_placeholder(self, "Leave empty for start")
        self.frame_begin_input.grid(row = 28, column = 0, columnspan = 3)
        
        frame_end_label = tk.Label(self, text = "End frame", font = LARGE_FONT)
        frame_end_label.grid(row = 27, column = 3, columnspan = 3)
        
        self.frame_end_input = entry_placeholder(self, "Leave empty for end")
        self.frame_end_input.grid(row = 28, column = 3, columnspan = 3)
        
        self.button_fit = my_button(self, text = "FIT", height = int(GUI_HEIGHT/8), 
                             width = int(GUI_WIDTH/8), state = 'disabled', 
                             command = lambda: self.start_fitting())#, style= 'my.TButton')
        self.button_fit.grid(row = 24, column = 6, columnspan = 3, rowspan = 5)
        
        
        progress_label = tk.Label(self, text = "Progress", font = LARGE_FONT)
        progress_label.grid(row = 24, column = 10)
        
        self.progress_status_label = tk.Label(self, text = "Not yet started", 
                                         bd = 1, relief = 'sunken')
        self.progress_status_label.grid(row = 25, column = 10)
        
        time_label = tk.Label(self, text = "Estimated time done", font = LARGE_FONT)
        time_label.grid(row = 27, column = 10)
        
        self.time_status_label = tk.Label(self, text = "Not yet started", 
                                         bd = 1, relief = 'sunken')
        self.time_status_label.grid(row = 28, column = 10)
        
        button_quit = ttk.Button(self, text = "Quit", command=quit_program)
        button_quit.grid(row = 50, column = 8)
        
        button_load_new = ttk.Button(self, text = "Load new", command= lambda: self.load_new(controller))
        button_load_new.grid(row = 50, column = 6)
        
        for i in range(0,10):
            self.columnconfigure(i, weight=1)

#%% START GUI
gui = mbx_python()
gui.geometry(str(GUI_WIDTH)+"x"+str(GUI_HEIGHT)+"+"+str(GUI_WIDTH_START)+"+"+str(GUI_HEIGHT_START))
gui.mainloop()