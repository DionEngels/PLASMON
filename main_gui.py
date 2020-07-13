# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
MBx Python Data Analysis

v0.1: first version GUI: 

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
WAVELENGTH = 637 #nm
THRESHOLD = 5 # X*sigma

#%% Initializations

FILETYPES = [("ND2", ".nd2")]

LARGE_FONT = ("Verdana", 12)
style.use("ggplot")

fit_options = ["Phasor with intensity", "Phasor without intensity",
                   "Gaussian with background", "Gaussian without background"]

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
#%% Declarations
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
        
def prints(param):
    print(param)
    
def quit_program():
    global gui
    gui.destroy()
    exit(0)
        
class my_button(ttk.Frame):
    def __init__(self, parent, height=None, width=None, text="", command=None, style=None):
        ttk.Frame.__init__(self, parent, height=height, width=width, style='my.TButton')

        self.pack_propagate(0)
        self._btn = ttk.Button(self, text=text, command=command, style=style)
        self._btn.pack(fill=tk.BOTH, expand=1)
        
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
        
class fitting_page(tk.Frame):
    
    def update_init(self, filenames):
        
        self.filenames = filenames
        dataset_roi_status = tk.Label(self, text = "Dataset 1 of " + str(len(filenames)))
        dataset_roi_status.grid(row = 9, column = 1, columnspan = 4) 
        
        roi_status = tk.Label(self, text = "0 of " + str(len(filenames)) + " have settings")
        roi_status.grid(row = 21, column = 0, columnspan = 6)
        
        
        with ND2_Reader(filenames[0]) as ND2:
            self.frames = ND2
            self.metadata = ND2.metadata
            self.fig_sub = plt.imshow(self.frames[0], extent=[0,self.frames[0].shape[1],self.frames[0].shape[0],0],
                       aspect='auto')
            
            canvas = FigureCanvasTkAgg(self.fig, self)
            canvas.draw()
            canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 16, sticky = 'E')
        
    
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)
        
        min_int_label = tk.Label(self, text = "Minimum Intensity", font = LARGE_FONT)
        min_int_label.grid(row = 0, column = 0, columnspan = 3)
        
        min_int_slider = tk.Scale(self, from_ = 0, to = 1000, orient = 'horizontal')
        min_int_slider.grid(row = 1, column = 0, columnspan = 3)
        
        max_int_label = tk.Label(self, text = "Maximum Intensity", font = LARGE_FONT)
        max_int_label.grid(row = 2, column = 0, columnspan = 3)
        
        max_int_slider = tk.Scale(self, from_ = 0, to = 5000, orient = 'horizontal')
        max_int_slider.grid(row = 3, column = 0, columnspan = 3)
        
        min_int_label = tk.Label(self, text = "Minimum Sigma", font = LARGE_FONT)
        min_int_label.grid(row = 0, column = 3, columnspan = 3)
        
        min_int_slider = tk.Scale(self, from_ = 0, to = 5, orient = 'horizontal')
        min_int_slider.grid(row = 1, column = 3, columnspan = 3)
        
        max_int_label = tk.Label(self, text = "Maximum Sigma", font = LARGE_FONT)
        max_int_label.grid(row = 2, column = 3, columnspan = 3)
        
        max_int_slider = tk.Scale(self, from_ = 0, to = 10, orient = 'horizontal')
        max_int_slider.grid(row = 3, column = 3, columnspan = 3)
        
        button_left = ttk.Button(self, text = "<<")
        button_left.grid(row = 9, column = 0)
        
        dataset_roi_status = tk.Label(self, text = "TBD")
        dataset_roi_status.grid(row = 9, column = 1, columnspan = 4)
        
        button_right = ttk.Button(self, text = ">>")
        button_right.grid(row = 9, column = 5)
        
        button_fit = ttk.Button(self, text = "Fit", command= lambda: prints("Joejoe"))
        button_fit.grid(row = 12, column = 0, columnspan = 2)
        
        button_save = ttk.Button(self, text = "Save these settings")
        button_save.grid(row = 12, column = 2, columnspan = 4)
                
        self.fig = Figure(figsize = (GUI_HEIGHT/DPI*0.7,GUI_HEIGHT/DPI*0.7), dpi = DPI)
        self.fig_sub = self.fig.add_subplot(111)
        # #fig_sub.plot([1,2,3,4,5], [3,5,1,5,2])
        
        canvas = FigureCanvasTkAgg(self.fig, self)
        canvas.draw()
        canvas.get_tk_widget().grid(row = 0, column = 10, rowspan = 16, sticky = 'E')
        
        # toolbar = NavigationToolbar2Tk(canvas, self)
        # toolbar.update()
        # # canvas._tkcanvas.grid(row = 1, column = 9)
        
        line = ttk.Separator(self, orient='horizontal')
        line.grid(column=0, row=18, rowspan = 2,  columnspan = 10, sticky='we')
        
        roi_status = tk.Label(self, text = "TBD")
        roi_status.grid(row = 21, column = 0, columnspan = 6)
        
        method_var = tk.StringVar(self)
        method_var.set(fit_options[0])
        
        method_drop = tk.OptionMenu(self, method_var, *fit_options)
        method_drop.grid(column = 0, row = 24, columnspan = 6)
        
        frame_begin_label = tk.Label(self, text = "Begin frame", font = LARGE_FONT)
        frame_begin_label.grid(row = 27, column = 0, columnspan = 3)
        
        frame_begin_input = entry_placeholder(self, "Leave empty for start")
        frame_begin_input.grid(row = 28, column = 0, columnspan = 3)
        
        frame_end_label = tk.Label(self, text = "End frame", font = LARGE_FONT)
        frame_end_label.grid(row = 27, column = 3, columnspan = 3)
        
        frame_end_input = entry_placeholder(self, "Leave empty for end")
        frame_end_input.grid(row = 28, column = 3, columnspan = 3)
        
        button_fit = my_button(self, text = "FIT", height = int(GUI_HEIGHT/8), 
                             width = int(GUI_WIDTH/8))#, style= 'my.TButton')
        button_fit.grid(row = 24, column = 6, columnspan = 3, rowspan = 5)
        
        
        progress_label = tk.Label(self, text = "Progress", font = LARGE_FONT)
        progress_label.grid(row = 24, column = 10)
        
        progress_status_label = tk.Label(self, text = "Not yet started", 
                                         bd = 1, relief = 'sunken')
        progress_status_label.grid(row = 25, column = 10)
        
        time_label = tk.Label(self, text = "Estimated time done", font = LARGE_FONT)
        time_label.grid(row = 27, column = 10)
        
        time_status_label = tk.Label(self, text = "Not yet started", 
                                         bd = 1, relief = 'sunken')
        time_status_label.grid(row = 28, column = 10)
        
        button_quit = ttk.Button(self, text = "Quit", command=quit_program)
        button_quit.grid(row = 50, column = 6, columnspan = 3)
        
        for i in range(0,10):
            self.columnconfigure(i, weight=1)

        

#%% START GUI
gui = mbx_python()
gui.geometry(str(GUI_WIDTH)+"x"+str(GUI_HEIGHT)+"+"+str(GUI_WIDTH_START)+"+"+str(GUI_HEIGHT_START))
gui.mainloop()