#!/usr/bin/env python3
"""Aerofoil pressure distribution visualisation tool

Authors:
- Stuart Johnston, 2013: Developed initial MATLAB version
- Joshua Smith, 2017: Revisions to MATLAB version
- Jarod Lam, 2020: Ported to Python
Version: 2.0
Copyright 2021 Queensland University of Technology

Developed for the QUT Young Accelerators "Aerodynamics and modelling lift" 
workshop.

Calculates the shape of an NACA 4 series airfoil and the corresponding pressure 
distribution according to the input using a vortex panel method approach. 
"""

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

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
        
        self.params = ParametersFrame(self.topFrame, self.calculate)
        self.params.pack(side=tk.LEFT, fill="both", padx=10, ipadx=10, pady=10)
        
        self.pressure = PressureFrame(self.topFrame)
        self.pressure.pack(side=tk.LEFT, fill="both")
        
        self.aerofoil = AerofoilFrame(self.bottomFrame)
        self.aerofoil.pack(fill="x")
        
        self.topFrame.pack(fill="both")
        self.bottomFrame.pack(fill="both", expand=True)
        
        self.reset()
        
    def calculate(self):
        # Get all parameters from text boxes
        c = self.params.fieldChordLen.get_value()
        t =  self.params.fieldThickness.get_value()
        m = self.params.fieldMaxCamber.get_value()
        p = self.params.fieldCamberPos.get_value()
        aoa = self.params.fieldAngleAttack.get_value()
        vinf = self.params.fieldVelocity.get_value()
        
        xa, ya, xc, yc, xmid, pd, L = lift(c, t, m, p, aoa, vinf)
        
    def reset(self):
        self.params.fieldChordLen.set_value(10)
        self.params.fieldThickness.set_value(20)
        self.params.fieldMaxCamber.set_value(20)
        self.params.fieldCamberPos.set_value(60)
        self.params.fieldAngleAttack.set_value(5)
        self.params.fieldVelocity.set_value(200)

class ParametersFrame(ttk.Frame):
    """Parameters selection in top-left of window"""
    
    def __init__(self, master, calculateCallback):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui(calculateCallback)
         
    def setup_gui(self, calculateCallback):
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
        self.fieldsFrame.grid(row=1, column=0, sticky=tk.NSEW)
        
        self.buttonsFrame = ttk.Frame(self)
        self.buttonCalculate = ttk.Button(self.buttonsFrame, text="Calculate", default="active", takefocus = 0, command=calculateCallback)
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
    
    def get_value(self):
        return int(self.entry.get())
        
    def set_value(self, value):
        self.entry.delete(0, tk.END)
        self.entry.insert(tk.END, value)

class PressureFrame(ttk.Frame):
    """Pressure distribution graph in top-right"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
        
    def setup_gui(self):
        self.figure = plt.Figure(figsize=(5, 3.5))
        t = np.arange(0, 3, .01)
        self.ax = self.figure.add_subplot()
        self.ax.plot(t, 2 * np.sin(2 * np.pi * t))
        self.ax.set_title("Pressure Distribution")
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class AerofoilFrame(ttk.Frame):
    """Aerofoil profile graph at bottom"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
        
    def setup_gui(self):
        self.figure = plt.Figure()
        t = np.arange(0, 3, .01)
        ax = self.figure.add_subplot()
        ax.plot(t, 2 * np.sin(2 * np.pi * t))
        ax.set_title("Aerofoil Design")
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

def lift(chord, thickness, camber, camberposition, angleofattack, velocity):
    """
    Function is used to calculate the shape of an NACA 4 series airfoil
    and the corresponding pressure distribution according to the input
    using a vortex panel method approach.
    """
    
    ## Rescaling of parameters
    # Thickness is percentage of cord (0-100)
    # Camber is percentage of cord(0-100)
    # Camberposition is percentage of cord (0-100)
    # Angleofattack is angle between chord and x axis (0-90)
    # Velocity is stream velocity in km/h (0-1000)
    c = chord                   # Chord length
    m = camber/100              # Maximum distance between camber and chord as fraction of chord
    p = camberposition/100      # Point where m occurs as fraction of chord
    t = thickness/100           # Thickness of airfoil as fraction of chord
    pc = p*c                    # Actual point where m occurs
    aoa = angleofattack*np.pi/180  # Angle the airfoil is inclined in radians
    vinf = velocity*1000/3600   # Stream velocity in m/s

    rhoinf = 1.225/1000         # Density at sea level at fifeen degrees
    pinf = 101.324              # Atmospheric pressure
    
    ## Plotting aerofoil
    # Set up grid
    N = 51                      # Number of grid points
    x1 = np.linspace(0,pc,N)    # Set up grid between start and m
    x2 = np.linspace(pc,c,N)    # Set up grid between m and end
    x2 = x2[2:]                 # Delete overlapping points
    
    # Calculate camber line
    yc1 = m*x1/p**2.*(2*p-x1/c)
    yc2 = m*(c-x2)/(1-p)**2.*(1+x2/c-2*p)
    
    yc = np.concatenate([yc1, yc2])
    xc = np.concatenate([x1, x2])
    
    # Calculate thickness
    yt = (t / 0.2) * c * (
        0.2969 * np.power((xc / c), 0.5) 
      - 0.126 * (xc / c) 
      - 0.3516 * np.power((xc / c), 2) 
      + 0.2843 * np.power((xc / c), 3) 
      - 0.1015 * np.power((xc / c), 4) 
    )
    
    # Angle of camber line
    theta = np.arctan2(yc, xc)
    
    xc = xc - c/2               # Rescale x to be centred at 0
    
    xu = xc - yt * np.sin(theta)   # Upper x values
    yu = yc + yt * np.cos(theta)   # Upper y values
    xl = xc + yt * np.sin(theta)   # Lower x values
    yl = yc - yt * np.cos(theta)   # Lower y values
    xa = np.concatenate([np.flip(xu), xl]) # Combine x values
    ya = np.concatenate([np.flip(yu), yl]) # Combine y values
    x = np.concatenate([xa[0:2*N-1], xa[2*N:]])  # Remove overlap
    y = np.concatenate([ya[0:2*N-1], ya[2*N:]])  # Remove overlap
    
    ## Calculate pressure distribution
    # n = len(x) - 1
    # a = np.zeros([n + 1, n + 1]);
    #
    # # Segment length
    # r = np.sqrt(
    #     (x[1:] - x[0:-1])**2
    #   + (y[1:]-y[0:-1])**2
    # )
    #
    # for j in range (0, n):
    #     a[j, n] = 1
    #     for i in range(0, n):
    #         if i == j:
    #             a[i, i] = r[i] / (2*np.pi) * (np.log(0.5*r[i]) - 1)
    #         else:
    #             xm1 = 0.5 * (x[j] + x[j+1])    # Segment x midpoint
    #             ym1 = 0.5 * (y[j] + y[j+1])    # Segment y midpoint
    #             dx = (x[i+1] - x[i]) / r[i]    # Rate of change of x over segment
    #             dy = (y[i+1] - y[i]) / r[i]    # Rate of change of y over segment
    #
    #             t1 = x[i] - xm1
    #             t2 = y[i] - ym1
    #             t3 = x[i+1] - xm1
    #             t7 = y[i+1] - ym1
    #             t4 = t1 * dx + t2 * dy
    #             t5 = t3 * dx + t7 * dy
    #             t6 = t2 * dx - t1 * dy
    #             t1 = t5 * np.log(t5**2 + t6**2) - t4 * np.log(t4**2 + t6**2)
    #             t2 = np.arctan2(t6,t4) - np.arctan2(t6,t5)
    #
    #             a[j, i] = (0.5*t1-t5+t4+t6*t2) / (2*np.pi)    # Influence coefficient matrix
    #
    #     a[n, 0] = 1
    #     a[n, n-1] = 1
    #
    # print(a)
    # print(a.shape)
    #
    # xmid = 0.5 * (x[0:-1] + x[1:]).transpose()    # Midpoint values of x
    # ymid = 0.5 * (y[0:-1] + y[1:]).transpose()    # Midpoint values of y
    # rhs = np.array([[ymid * np.cos(aoa) - xmid * np.sin(aoa)], [0]])  # Right hand side of matrix
    # gamma = np.linalg.solve(a, rhs)                                   # Solve linear system
    # cp = 1 - gamma[0:-1]**2                       # Calculate pressure coefficient
    # p = cp*0.5*rhoinf*vinf^2+pinf                 # Solve for pressure distribution
    # Ll = 0
    # Lu = 0
    #
    # for i in range(0, 99):
    #     Lu = Lu + 0.5 * (xmid[i+1+100] - xmid[i+100]) * (p[i+100] + p[i+1+100])
    #     Ll = Ll + 0.5 * (xmid[i] - xmid[i+1]) * (p[i] + p[i+1])
    #
    # Lift = Lu - Ll
    # Lift = Lift * 1000
    #
    # # Rescale to start at 0
    # xa = xa + c/2
    # xc = xc + c/2
    # xmid = xmid + c/2
    
    xmid = 0
    p = 0
    Lift = 0
    
    return (xa, ya, xc, yc, xmid, p, Lift)

if __name__ == "__main__":
    # Set up Tk instance
    root = tk.Tk()
    root.title("Aerofoil GUI")
    root.geometry('800x600')
    root.minsize(800, 600)
    root.maxsize(800, 600)

    # Create main frame
    mainFrame = MainFrame(root)
    mainFrame.pack(fill="both", expand=True)

    # Run
    root.mainloop()