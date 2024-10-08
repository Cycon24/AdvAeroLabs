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
mu = 1.73e-5
l = 6/39.27

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
xs = np.array([0., 2.5, 5., 10., 20., 30., 40., 50., 70., 90., 85.8, 70., 50., 30., 20., 10., 6., 4.])
ys = np.array([0., 3.4, 4.7, 6.6, 8.8, 9.8, 9.8, 9.2, 6.7, 2.7, -0.3, -0.7, -1.4, -2.3, -2.7, -2.9, -2.6, -2.3])

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
    plt.plot(xs/100*scale, CPs)
    plt.plot(xps, -yps, "k-")
    plt.xlim(0, 2); plt.xlabel("$x/c$")
    plt.xticks(np.linspace(0, scale, 11), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.ylim(1, -2.5); plt.ylabel("$C_P$")
    plt.title(f"AOA={aoa:.0f}$^\circ$", fontsize=14)
    plt.gca().set_aspect("equal")
    plt.grid(); plt.tight_layout()

# # do the plotting (part 2)
# presDist(pressAOA0, 0)
# presDist(pressAOA8, 8)
# presDist(pressAOA10, 10)
# presDist(pressAOA12, 12)

def clcdCalc(press, aoa):
    # determine q
    q = np.mean(press[0][:])
    
    # split pressure data to upper and lower surfaces
    CPsu = np.append(press[1][:10]/q, 0)
    CPsl = np.flip(np.append(0, press[1][10:18]/q))
    
    # split tap locations into upper and lower as well
    xsu = np.append(xs[:10]/100, 1)
    xsl = np.flip(np.append(1, xs[10:18]/100))
    ysu = np.append(ys[:10]/100, 0)
    ysl = np.flip(np.append(0, ys[10:18]/100))

    # function to determine force coefficients (since all the summations are essentially the same)
    def coef(CPs, xs):
        coef = 0.
        
        # loop to do numerical integration (summation)
        for i in range(len(CPs) - 1):
            coef += 0.5*(CPs[i + 1] + CPs[i])*(xs[i + 1] - xs[i])

        return coef
    
    # calculate Cn for the upper and lower surface using force coeff fxn
    Cnu = coef(CPsu, xsu)
    Cnl = coef(CPsl, xsl)
    
    # calculate Cc for the upper and lower surface using force coeff fxn
    Ccu = coef(CPsu, ysu)
    Ccl = coef(CPsl, ysl)
    
    # determine actual Cn and Cc
    Cn = Cnl - Cnu
    Cc = Ccu - Ccl
    print(f"aoa: {aoa:.0f}")
    print(f"Cn: {Cn:.3f}, Cc: {Cc:.3f}")
    print(f"Cnu: {Cnu:.3f}, Cnl: {Cnl:.3f}, Ccu: {Ccu:.3f}, Ccl: {Ccl:.3f}")
    # calculate lift and drag coeffs
    Cl = Cn*np.cos(np.deg2rad(aoa)) - Cc*np.sin(np.deg2rad(aoa))
    Cd = Cc*np.cos(np.deg2rad(aoa)) + Cn*np.sin(np.deg2rad(aoa))
    print(f"Cl: {Cl:.3f}, Cd: {Cd:.3f}")
    print()
    return Cl, Cd

# # calc cl/cd for each airfoil
# Cl0, Cd0 = clcdCalc(pressAOA0, 0)
# Cl4, Cd4 = clcdCalc(pressAOA4, 4)
# Cl6, Cd6 = clcdCalc(pressAOA6, 6)
# Cl8, Cd8 = clcdCalc(pressAOA8, 8)
# Cl10, Cd10 = clcdCalc(pressAOA10, 10)
# Cl12, Cd12 = clcdCalc(pressAOA12, 12)
# Cl14, Cd14 = clcdCalc(pressAOA14, 14)
# Cls = np.array([Cl0, Cl4, Cl6, Cl8, Cl10, Cl12, Cl14])
# Cds = np.array([Cd0, Cd4, Cd6, Cd8, Cd10, Cd12, Cd14])
# aoas = np.array([0, 4, 6, 8, 10, 12, 14])

# # Cl/Cd plots
# plt.figure()
# plt.plot(aoas, Cls)
# plt.plot(aoas, Cds)
# # plt.plot(aoas, Cls/Cds)
# plt.xlim(0, 14); plt.xlabel("AOA (deg)")


# calculate U for each case
q0 = np.mean(pressAOA0[0][:])
q4 = np.mean(pressAOA4[0][:])
q6 = np.mean(pressAOA6[0][:])
q8 = np.mean(pressAOA8[0][:])
q10 = np.mean(pressAOA10[0][:])
q12 = np.mean(pressAOA12[0][:])
q14 = np.mean(pressAOA14[0][:])
qs = np.array([q0, q4, q6, q8, q10, q12, q14])*248.84
Us = np.sqrt(2*qs/rho)

# calculate Re for each case
Res = rho*Us*l/mu

# determine stagnation


# plt.show()