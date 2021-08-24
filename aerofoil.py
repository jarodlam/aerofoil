 #!/usr/bin/env python3
"""Aerofoil pressure distribution visualisation tool

Authors:
- Stuart Johnston, 2013: Developed initial MATLAB version
- Joshua Smith, 2017: Revisions to MATLAB version
- Jarod Lam, 2021: Ported to Python
Version: 2.0
Copyright 2021 Queensland University of Technology

Developed for the QUT Young Accelerators "Aerodynamics and modelling lift" 
workshop.

Calculates the shape of an NACA 4 series airfoil and the corresponding pressure 
distribution according to the input using a vortex panel method approach. 

Requires Tkinter, NumPy, Matplotlib, OpenPyXL and PyInstaller

To package into a single executable, run:
pyinstaller --onefile --windowed --icon icon.ico aerofoil.py
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from openpyxl import Workbook

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
        
        self.params = ParametersFrame(self.topFrame, self.calculate, self.reset, self.output)
        self.params.pack(side=tk.LEFT, anchor=tk.N, fill="x", padx=10, pady=10)
        
        self.pressure = PressureFrame(self.topFrame)
        self.pressure.pack(side=tk.LEFT, fill="both", expand=True)
        
        self.aerofoil = AerofoilFrame(self.bottomFrame)
        self.aerofoil.pack(fill="x")
        
        self.topFrame.pack(fill="both")
        self.bottomFrame.pack(fill="both", expand=True)
        
        self.reset()
        
    def calculate(self):
        if not self.params.all_fields_are_valid():
            return
        
        # Get all parameters from text boxes
        c = self.params.fieldChordLen.get_value()
        t =  self.params.fieldThickness.get_value()
        m = self.params.fieldMaxCamber.get_value()
        p = self.params.fieldCamberPos.get_value()
        aoa = self.params.fieldAngleAttack.get_value()
        vinf = self.params.fieldVelocity.get_value()
        
        # Use lift function to calculate points
        xa, ya, xc, yc, xmid, pd, L = lift(c, t, m, p, aoa, vinf)
        
        # Plot values in their respective frames
        self.aerofoil.update_plot(xa, ya, xc, yc)
        self.pressure.update_plot(xmid, pd)
        
        # Output total force
        self.params.set_total_force(L)
        
        # Save pressure for outputting
        self.xmid = xmid
        self.pd = pd
        
        # Enable output button so we can export data (if not already enabled)
        self.params.set_output_button_enable(True)
        
    def reset(self):
        """Reset all values to default state"""
        
        self.params.fieldChordLen.set_value(10)
        self.params.fieldThickness.set_value(12)
        self.params.fieldMaxCamber.set_value(2)
        self.params.fieldCamberPos.set_value(40)
        self.params.fieldAngleAttack.set_value(5)
        self.params.fieldVelocity.set_value(200)
        self.params.set_output_button_enable(False)
        self.params.set_total_force(None)
        self.aerofoil.reset_plot()
        self.pressure.reset_plot()
    
    def output(self):
        """Save pressure distribution to Excel workbook
        Format is: x | upper pressure | x | lower pressure
        """
        
        # Open a "Save as" dialog box to get file name
        filepath = tk.filedialog.asksaveasfilename(
            confirmoverwrite=True,
            defaultextension=".xlsx",
            filetypes=[("Excel workbook", ".xlsx")],
            initialfile="aerofoil.xlsx",
            parent=self.master
        )
        
        # User pressed Cancel
        if not filepath:
            return
        
        # Set up Excell workbook
        wb = Workbook()
        ws = wb.active
        
        # Write data
        ws.append(["x", "Under aerofoil pressure (kPa)", "x", "Above aerofoil pressure (kPa)"])
        data = self.get_data_average()
        for row in data:
            ws.append(row.tolist())
        
        # Save the file
        wb.save(filename=filepath)
    
    def get_data_full(self):
        mid = int(np.ceil(len(self.xmid)/2))
        xu, xa = self.xmid[mid:], self.xmid[0:mid]
        pu, pa = self.pd[mid:], self.pd[0:mid]
        data = np.transpose([xu, pu, xa, pa])
        return data
    
    def get_data_average(self):
        data_full = self.get_data_full()
        
        n = 10    # Number of data points
        data_split = np.array_split(data_full, n)
        data = []
        for d in data_split:
            data.append(np.mean(d, axis=0))
        
        data = np.array(data)
        return data

class ParametersFrame(ttk.Frame):
    """Parameters selection in top-left of window"""
    
    def __init__(self, master, calculateCallback, resetCallback, outputCallback):
        self.master = master
        ttk.Frame.__init__(self, master)
        
        self.setup_gui(calculateCallback, resetCallback, outputCallback)
         
    def setup_gui(self, calculateCallback, resetCallback, outputCallback):
        self.fieldsFrame = ttk.LabelFrame(self, text="Parameters")
        self.fieldChordLen = ParametersField(self.fieldsFrame, "Chord length (1-100)", "m", 1, 100)
        self.fieldChordLen.pack(expand=True, padx=5, pady=5)
        self.fieldThickness = ParametersField(self.fieldsFrame, "Thickness (1-40)", "%", 1, 40)
        self.fieldThickness.pack(expand=True, padx=5, pady=5)
        self.fieldMaxCamber = ParametersField(self.fieldsFrame, "Maximum camber (0-10)", "%", 0, 10)
        self.fieldMaxCamber.pack(expand=True, padx=5, pady=5)
        self.fieldCamberPos = ParametersField(self.fieldsFrame, "Camber position (1-90)", "%", 1, 90)
        self.fieldCamberPos.pack(expand=True, padx=5, pady=5)
        self.fieldAngleAttack = ParametersField(self.fieldsFrame, "Angle of attack (0-10)", "°", 0, 10)
        self.fieldAngleAttack.pack(expand=True, padx=5, pady=5)
        self.fieldVelocity = ParametersField(self.fieldsFrame, "Velocity (1-10000)", "km/h", 1, 10000)
        self.fieldVelocity.pack(expand=True, padx=5, pady=5)
        self.fieldsFrame.grid(row=1, column=0, sticky=tk.NSEW)
        
        self.buttonsFrame = ttk.Frame(self)
        self.buttonCalculate = ttk.Button(self.buttonsFrame, text="Calculate", state=tk.NORMAL, default=tk.ACTIVE, takefocus=False, command=calculateCallback)
        self.buttonCalculate.pack(side=tk.LEFT, pady=5, expand=True)
        self.buttonReset = ttk.Button(self.buttonsFrame, text="Reset", state=tk.NORMAL, takefocus = 0, command=resetCallback)
        self.buttonReset.pack(side=tk.LEFT, pady=5, expand=True)
        self.buttonOutput = ttk.Button(self.buttonsFrame, text="Save", state=tk.DISABLED, takefocus = 0, command=outputCallback)
        self.buttonOutput.pack(side=tk.LEFT, pady=5, expand=True)
        self.buttonsFrame.grid(row=2, column=0, sticky="nsew")
        
        self.forceFrame = ttk.LabelFrame(self, text="Total force")
        self.forceFrame.grid(row=3, column=0, sticky="nsew", pady=5)
        self.totalForce = ttk.Entry(self.forceFrame, state="readonly", justify="right")
        self.totalForce.pack(side=tk.LEFT, expand=True, fill="both", padx=(5,0), pady=5)
        self.labelNewtons = ttk.Label(self.forceFrame, text="N", width=5)
        self.labelNewtons.pack(side=tk.LEFT, padx=(0,5), pady=5)
    
    def set_output_button_enable(self, enable):
        self.buttonOutput["state"] = tk.NORMAL if enable else tk.DISABLED
        
    def set_total_force(self, newtons):
        self.totalForce["state"] = "normal"
        self.totalForce.delete(0, tk.END)
        if newtons is not None:
            self.totalForce.insert(tk.END, "{:.0f}".format(newtons))
        self.totalForce["state"] = "readonly"
    
    def all_fields_are_valid(self):
        return \
            self.fieldChordLen.is_within_bounds() and \
            self.fieldThickness.is_within_bounds() and \
            self.fieldMaxCamber.is_within_bounds() and \
            self.fieldCamberPos.is_within_bounds() and \
            self.fieldAngleAttack.is_within_bounds() and \
            self.fieldVelocity.is_within_bounds()
        
class ParametersField(ttk.Frame):
    """A label, text field, and unit label"""
    
    def __init__(self, master, label, unit, minVal, maxVal):
        self.master = master
        ttk.Frame.__init__(self, master)
        
        self.label = label
        self.unit = unit
        self.minVal = minVal
        self.maxVal = maxVal
        
        self.setup_gui()
    
    def setup_gui(self):
        self.s = ttk.Style()
        self.styleName = self.label + ".TSpinbox"
        self.s.configure(self.styleName, background='white')
        
        self.label = ttk.Label(self, text=self.label, width=25)
        self.label.pack(side=tk.LEFT)
        
        vcmd = (self.register(self.validate), '%P')
        self.entry = ttk.Spinbox(self, width=5, style=self.styleName, from_=self.minVal, to=self.maxVal, validate="key", validatecommand=vcmd)
        self.entry.pack(side=tk.LEFT)
        
        self.unit = ttk.Label(self, text=self.unit, width=5)
        self.unit.pack(side=tk.LEFT)
    
    def get_value(self):
        try:
            return int(self.entry.get())
        except ValueError:
            return None
        
    def set_value(self, value):
        self.entry.delete(0, tk.END)
        self.entry.insert(tk.END, value)

    def validate(self, P):
        """Validate input while user is typing"""
        
        # Allow blank field
        if P == "":
            self.set_entry_display_error(False)
            return True
        
        # Check if is an integer
        try:
            value = int(P)
            self.set_entry_display_error(not self.is_within_bounds(value))
            return True
        except ValueError:
            return False
    
    def set_entry_display_error(self, error):
        colour = "red" if error else "black"
        self.s.configure(self.styleName, foreground=colour)
    
    def is_within_bounds(self, value=None):
        if value is None:
            value = self.get_value()
        
        if value is None or value == "":
            return False
        
        return value >= self.minVal and value <= self.maxVal

class PressureFrame(ttk.Frame):
    """Pressure distribution graph in top-right"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
        
    def setup_gui(self):
        self.figure, self.ax = plt.subplots()
        self.setup_plot()
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def setup_plot(self):
        self.ax.set_title("Pressure distribution")
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("Pressure (kPa)")
        self.ax.legend()
    
    def update_plot(self, xmid, pd):
        self.ax.clear()
        mid = int(np.ceil(len(xmid)/2))
        self.ax.plot(xmid[mid:], pd[mid:], "r", label="Under aerofoil pressure", linewidth=1)
        self.ax.plot(xmid[0:mid], pd[0:mid], "g", label="Above aerofoil pressure", linewidth=1)
        self.setup_plot()
        self.canvas.draw()
    
    def reset_plot(self):
        self.ax.clear()
        self.setup_plot()
        self.canvas.draw()

class AerofoilFrame(ttk.Frame):
    """Aerofoil profile graph at bottom"""
    
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master, borderwidth = 1)
        
        self.setup_gui()
        
    def setup_gui(self):
        self.figure, self.ax = plt.subplots()
        self.setup_plot()
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def setup_plot(self):
        self.ax.set_title("Aerofoil design")
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")
        self.figure.subplots_adjust(bottom=0.2)
        self.ax.legend()
        self.ax.set_aspect('equal')
    
    def update_plot(self, xa, ya, xc, yc):
        self.ax.clear()
        self.ax.plot(xa, ya, "b", label="Aerofoil surface", linewidth=1)
        self.ax.plot(xc, yc, "r", label="Mean camber line", linewidth=1)
        self.setup_plot()
        self.canvas.draw()
    
    def reset_plot(self):
        self.ax.clear()
        self.setup_plot()
        self.canvas.draw()

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
    x2 = x2[1:]                 # Delete overlapping points
    
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
    n = len(x) - 1
    a = np.zeros([n+1, n+1]);
    
    # Segment length
    r = np.sqrt(
        (x[1:] - x[0:-1])**2
      + (y[1:]-y[0:-1])**2
    )
    
    for j in range (0, n):
        a[j, n] = 1
        for i in range(0, n):
            if i == j:
                a[i, i] = r[i] / (2*np.pi) * (np.log(0.5*r[i]) - 1)
            else:
                xm1 = 0.5 * (x[j] + x[j+1])    # Segment x midpoint
                ym1 = 0.5 * (y[j] + y[j+1])    # Segment y midpoint
                dx = (x[i+1] - x[i]) / r[i]    # Rate of change of x over segment
                dy = (y[i+1] - y[i]) / r[i]    # Rate of change of y over segment
    
                t1 = x[i] - xm1
                t2 = y[i] - ym1
                t3 = x[i+1] - xm1
                t7 = y[i+1] - ym1
                t4 = t1 * dx + t2 * dy
                t5 = t3 * dx + t7 * dy
                t6 = t2 * dx - t1 * dy
                t1 = t5 * np.log(t5**2 + t6**2) - t4 * np.log(t4**2 + t6**2)
                t2 = np.arctan2(t6,t4) - np.arctan2(t6,t5)
    
                a[j, i] = (0.5*t1-t5+t4+t6*t2) / (2*np.pi)    # Influence coefficient matrix
    
        a[n, 0] = 1
        a[n, n-1] = 1
    
    xmid = 0.5 * (x[0:-1] + x[1:]).transpose()    # Midpoint values of x
    ymid = 0.5 * (y[0:-1] + y[1:]).transpose()    # Midpoint values of y
    rhs = ymid * np.cos(aoa) - xmid * np.sin(aoa) # Right hand side of matrix
    rhs = np.append(rhs, 0)
    gamma = np.linalg.solve(a, rhs)                                   # Solve linear system
    cp = 1 - gamma[0:-1]**2                       # Calculate pressure coefficient
    p = cp*0.5*rhoinf*vinf**2+pinf                 # Solve for pressure distribution
    Ll = 0
    Lu = 0
    
    for i in range(0, 99):
        Lu = Lu + 0.5 * (xmid[i+1+100] - xmid[i+100]) * (p[i+100] + p[i+1+100])
        Ll = Ll + 0.5 * (xmid[i] - xmid[i+1]) * (p[i] + p[i+1])
    
    Lift = Lu - Ll
    Lift = Lift * 1000
    
    # Rescale to start at 0
    xa = xa + c/2
    xc = xc + c/2
    xmid = xmid + c/2
    
    return (xa, ya, xc, yc, xmid, p, Lift)

if __name__ == "__main__":
    # Set up Tk instance
    root = tk.Tk()
    root.title("QUT Young Accelerators aerofoil")
    root.geometry('800x800')
    root.minsize(800, 800)

    # Create main frame
    mainFrame = MainFrame(root)
    mainFrame.pack(fill="both", expand=True)

    # Run
    root.mainloop()
