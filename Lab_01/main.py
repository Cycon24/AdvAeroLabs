# Adv. Aero Lab 1 Plotting
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# constants
rho = 1.225
R = 287
T = 300.372
y = 1.4
Ps = 98397.16

# import data from files
presDF = pd.read_csv("Lab_01/NACA_4412_Data.csv")
travDF = pd.read_csv("Lab_01/traverse.csv")

# converting to numpy arrays
# [0][:] returns dynamic press (psi), [1][:] returns surface pres (psi)
pressAOA0 = np.array([presDF["Unnamed: 1"][3:21], presDF["Unnamed: 2"][3:21]], dtype="float64")
pressAOA4 = np.array([presDF["Unnamed: 1"][25:43], presDF["Unnamed: 2"][25:43]], dtype="float64")
pressAOA6 = np.array([presDF["Unnamed: 1"][47:65], presDF["Unnamed: 2"][47:65]], dtype="float64")
pressAOA8 = np.array([presDF["Unnamed: 1"][69:87], presDF["Unnamed: 2"][69:87]], dtype="float64")
pressAOA10 = np.array([presDF["Unnamed: 1"][91:109], presDF["Unnamed: 2"][91:109]], dtype="float64")
pressAOA12 = np.array([presDF["Unnamed: 1"][113:131], presDF["Unnamed: 2"][113:131]], dtype="float64")
pressAOA14 = np.array([presDF["Unnamed: 1"][135:153], presDF["Unnamed: 2"][135:153]], dtype="float64")
# [0][:] returns traverse locations, [1][:] returns velocity (ft/s)
velsAOA0 = np.array([travDF["AOA"][1:], travDF["0"][1:]], dtype="float64")
velsAOA8 = np.array([travDF["AOA"][1:], travDF["8"][1:]], dtype="float64")
velsAOA12 = np.array([travDF["AOA"][1:], travDF["12"][1:]], dtype="float64")

# plot scale factor
scale = 3

# pressure tap coordinates
xs = np.array([0., 2.5, 5., 10., 20., 30., 40., 50., 70., 90., 85.8, 70., 50., 30., 20., 10., 6., 4.])/100*scale
ys = np.array([0., 3.4, 4.7, 6.6, 8.8, 9.8, 9.8, 9.2, 6.7, 2.7, -0.3, -0.7, -1.4, -2.3, -2.7, -2.9, -2.6, -2.3])/100*scale

# airfoil plotting
xps = np.array([1.0000, 0.9500, 0.9000, 0.8000, 0.7000, 0.6000, 0.5000, 0.4000, 0.3000, 0.2500, 
              0.2000, 0.1500, 0.1000, 0.0750, 0.0500, 0.0250, 0.0125, 0.0000, 0.0125, 0.0250, 
              0.0500, 0.0750, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.4000, 0.5000, 0.6000, 
              0.7000, 0.8000, 0.9000, 0.9500, 1.0000, 1.])*scale
yps = np.array([0.0013, 0.0147, 0.0271, 0.0489, 0.0669, 0.0814, 0.0919, 0.0980, 0.0976, 0.0941, 
              0.0880, 0.0789, 0.0659, 0.0576, 0.0473, 0.0339, 0.0244, 0.0000, -0.0143, -0.0195, 
              -0.0249, -0.0274, -0.0286, -0.0288, -0.0274, -0.0250, -0.0226, -0.0180, -0.0140, 
              -0.0100, -0.0065, -0.0039, -0.0022, -0.0016, -0.0013, 0.0013])*scale

# function for plotting pressure dist.
def presDist(press, aoa):
    # determine q
    q = np.mean(press[0][:])
    
    # calculate Cp for each location
    CPs = press[1][:]/q
    
    # plot it
    plt.figure(f"CP Dist. AOA={aoa:.0f} deg (2)")
    plt.plot(xs, CPs)
    plt.plot(xps, -yps, "k-")
    plt.xlim(0, 2); plt.xlabel("$x/c$")
    plt.xticks(np.linspace(0, scale, 11), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.ylim(1, -2.5); plt.ylabel("$C_P$")
    plt.title(f"AOA={aoa:.0f}$^\circ$", fontsize=14)
    plt.gca().set_aspect("equal")
    plt.grid(); plt.tight_layout()

# do the plotting (part 2)
presDist(pressAOA0, 0)
presDist(pressAOA8, 8)
presDist(pressAOA10, 10)
presDist(pressAOA12, 12)

plt.show()