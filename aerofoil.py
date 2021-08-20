#!/usr/bin/env python3
"""
Aerofoil pressure distribution visualisation tool for QUT Young Accelerators "Aerodynamics" workshop
"""
import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

__author__ = "Jarod Lam"
__copyright__ = "Copyright 2021 Queensland University of Techonology"
__verson__ = "1.0"
__email__ = "accelerators@qut.edu.au"

class MainFrame(ttk.Frame):
    """Master Tk frame"""
    
    def __init__(self, master):
         self.master = master
         ttk.Frame.__init__(self, master)
         self.setup_gui()
         self.pack()

    def setup_gui(self):
        self.topFrame = ttk.Frame(self)
        self.bottomFrame = ttk.Frame(self)
        
        self.params = ParametersFrame(self.topFrame)
        self.params.pack(side=tk.LEFT, fill="both", padx=10, ipadx=10, pady=10)
        
        self.pressure = PressureFrame(self.topFrame)
        self.pressure.pack(side=tk.LEFT, fill="both")
        
        self.airfoil = AirfoilFrame(self.bottomFrame)
        self.airfoil.pack(fill="x")
        
        self.topFrame.pack(fill="both")
        self.bottomFrame.pack(fill="both", expand=True)

class ParametersFrame(ttk.Frame):
    """Parameters selection in top-left of window"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
         
    def setup_gui(self):
        s = ttk.Style()
        s.configure('Calculate.TButton', foreground='green')
        
        self.label1 = ttk.Label(self, text="Parameters")
        self.label1.grid(row=0, column=0)
        
        self.fieldsFrame = ttk.Frame(self)
        self.fieldChordLen = ParametersField(self.fieldsFrame, "Chord length", "m")
        self.fieldChordLen.pack(expand=True)
        self.fieldThickness = ParametersField(self.fieldsFrame, "Thickness", "%")
        self.fieldThickness.pack(expand=True)
        self.fieldMaxCamber = ParametersField(self.fieldsFrame, "Maximum camber", "%")
        self.fieldMaxCamber.pack(expand=True)
        self.fieldCamberPos = ParametersField(self.fieldsFrame, "Camber position", "%")
        self.fieldCamberPos.pack(expand=True)
        self.fieldAngleAttack = ParametersField(self.fieldsFrame, "Angle of attack", "°")
        self.fieldAngleAttack.pack(expand=True)
        self.fieldVelocity = ParametersField(self.fieldsFrame, "Velocity", "km/h")
        self.fieldVelocity.pack(expand=True)
        self.fieldDataPoints = ParametersField(self.fieldsFrame, "# Data points", "")
        self.fieldDataPoints.pack(expand=True)
        self.fieldsFrame.grid(row=1, column=0, sticky=tk.NSEW)
        
        self.buttonsFrame = ttk.Frame(self)
        self.buttonCalculate = ttk.Button(self.buttonsFrame, text="Calculate", default="active", takefocus = 0)
        self.buttonCalculate.pack(side=tk.LEFT)
        self.buttonReset = ttk.Button(self.buttonsFrame, text="Reset", takefocus = 0)
        self.buttonReset.pack(side=tk.LEFT)
        self.buttonReset = ttk.Button(self.buttonsFrame, text="Output Data", takefocus = 0)
        self.buttonReset.pack(side=tk.LEFT)
        self.buttonsFrame.grid(row=2, column=0)
        
        self.totalForce = ttk.Label(self, text="Total force: ")
        self.totalForce.grid(row=3, column=0, sticky="nsew")
        
class ParametersField(ttk.Frame):
    """A label, text field, and unit label"""
    
    def __init__(self, master, label, unit):
        self.master = master
        ttk.Frame.__init__(self, master)
        
        self.label = label
        self.unit = unit
        self.setup_gui()
    
    def setup_gui(self):
        s = ttk.Style()
        s.configure('TEntry', background='white')
        
        self.label = ttk.Label(self, text=self.label, width=20)
        self.label.pack(side=tk.LEFT)
        
        self.entry = ttk.Entry(self, width=5)
        self.entry.pack(side=tk.LEFT)
        
        self.unit = ttk.Label(self, text=self.unit, width=4)
        self.unit.pack(side=tk.LEFT)

class PressureFrame(ttk.Frame):
    """Pressure distribution graph in top-right"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
        
    def setup_gui(self):
        self.figure = plt.Figure(figsize=(5, 3.5))
        t = np.arange(0, 3, .01)
        ax = self.figure.add_subplot()
        ax.plot(t, 2 * np.sin(2 * np.pi * t))
        ax.set_title("Pressure Distribution")
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class AirfoilFrame(ttk.Frame):
    """Airfoil profile graph at bottom"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
        
    def setup_gui(self):
        self.figure = plt.Figure()
        t = np.arange(0, 3, .01)
        ax = self.figure.add_subplot()
        ax.plot(t, 2 * np.sin(2 * np.pi * t))
        ax.set_title("Airfoil Design")
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

if __name__ == "__main__":
    # Set up Tk instance
    root = tk.Tk()
    root.title("Airfoil GUI")
    root.geometry('800x600')
    root.minsize(800, 600)
    root.maxsize(800, 600)

    # Create main frame
    mainFrame = MainFrame(root)
    mainFrame.pack(fill="both", expand=True)

    # Run
    root.mainloop()