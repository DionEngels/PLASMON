# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:14:16 2020

@author: Dion Engels
MBx Python Data Analysis

v0.1: first version GUI: 

"""

## GENERAL IMPORTS
import os # to get standard usage
#import math # for generic math
import time # for timekeeping

## Numpy and matplotlib, for linear algebra and plotting respectively
import numpy as np
import matplotlib.pyplot as plt

## GUI
import tkinter as tk # for GUI
from tkinter.filedialog import askopenfilenames # for popup that asks to select .nd2's

## ND2 related
from pims import ND2_Reader # reader of ND2 files
#from pims import Frame # for converting ND2 generator to one numpy array

## Multiprocessing
import multiprocessing as mp

## comparison to other fitters
import scipy.io

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

#%% Declarations
class gui(tk.Tk):
    
    def __init__(self, *args, **kwargs):
    
        tk.Tk.__init__(self, *args, **kwargs)
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
        
        self.show_frame(load_page)
        
    def show_frame(self, cont):
        
        frame = self.frames[cont]
        frame.tkraise()
        
def prints(param):
    print(param)
        
class load_page(tk.Frame):
    
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text = "Start Page", font = LARGE_FONT)
        label.pack(pady = 10, padx = 10)
        
        button1 = tk.Button(self, text = "LOAD", command= lambda: self.load_nd2(controller))
        button1.pack()
        
    def load_nd2(self, controller):
        
        filenames = askopenfilenames(filetypes = FILETYPES, 
                                     title = "Select file", 
                                     initialdir =os.getcwd())
        
        controller.show_frame(fitting_page)
        
class fitting_page(tk.Frame):
    
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text = "Fitting Page", font = LARGE_FONT)
        label.pack(pady = 10, padx = 10)
        
        button1 = tk.Button(self, text = "DO FITTING THINGIES", command= lambda: prints("Joejoe"))
        button1.pack()

#%% START GUI
gui = gui()
gui.mainloop()
# filenames = askopenfilenames(filetypes = FILETYPES,title = "Select file", initialdir =os.getcwd()) 
# show an "Open" dialog of only .nd2's starting at code location