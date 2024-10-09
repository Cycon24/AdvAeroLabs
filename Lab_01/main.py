# Adv. Aero Lab 1 Plotting
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# set constants
R = 287             # R of air (SI)
T = 300.372         # tunnel temp (K)
y = 1.4             # gamma baby
mu = 1.73e-5        # viscosity of air (SI)
l = 6/39.27         # airfoil length (m)
pinf = 98397.16     # barometric pressure (Pa)
rho = pinf/R/T      # density (kg/m^3)

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate table of pressure coefficient values
def presDist(press):
    # determine q
    q = np.mean(press[0][:])

    # calculate Cp for each location
    CPs = press[1][:]/q

    return CPs

CPdistAOA0 = presDist(pressAOA0)
CPdistAOA4 = presDist(pressAOA4)
CPdistAOA6 = presDist(pressAOA6)
CPdistAOA8 = presDist(pressAOA8)
CPdistAOA10 = presDist(pressAOA10)
CPdistAOA12 = presDist(pressAOA12)
CPdistAOA14 = presDist(pressAOA14)

pressDistVals = {"AOA 0": CPdistAOA0, 
                 "AOA 4": CPdistAOA4, 
                 "AOA 6": CPdistAOA6,
                 "AOA 8": CPdistAOA8,
                 "AOA 10": CPdistAOA10,
                 "AOA 12": CPdistAOA12,
                 "AOA 14": CPdistAOA14,
                 }
pressDistVals = pd.DataFrame(pressDistVals)
print(pressDistVals)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function for plotting pressure distribution for AOAs 0, 8, 12
def presDistPlot(CPs, aoa):

    plt.figure(f"CP Dist. AOA={aoa:.0f} deg (2)")
    plt.plot(xs/100*scale, CPs)
    plt.plot(xps, -yps, "k-")
    plt.xlim(0, 2); plt.xlabel("$x/c$")
    plt.xticks(np.linspace(0, scale, 11), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.ylim(1, -2.5); plt.ylabel("$C_P$")
    plt.title(f"AOA={aoa:.0f}$^\circ$", fontsize=14)
    plt.gca().set_aspect("equal")
    plt.grid(); plt.tight_layout()

# do the plotting
presDistPlot(CPdistAOA0, 0)
presDistPlot(CPdistAOA8, 8)
presDistPlot(CPdistAOA10, 10)
presDistPlot(CPdistAOA12, 12)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def clcdCalc(press, aoa):
    # determine q
    q = np.mean(press[0][:])
    
    # split pressure data to upper and lower surfaces
    CPsu = np.append(press[1][:10]/q, 0)
    CPsl = np.flip(np.append(0, press[1][10:18]/q))
    
    # split tap locations into upper and lower as well, also change lower surface direction to LE -> TE
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
    
    # determine total Cn and Cc
    Cn = Cnl - Cnu
    Cc = Ccu - Ccl

    # calculate lift and drag coeffs
    Cl = Cn*np.cos(np.deg2rad(aoa)) - Cc*np.sin(np.deg2rad(aoa))
    Cd = Cc*np.cos(np.deg2rad(aoa)) + Cn*np.sin(np.deg2rad(aoa))

    return Cl, Cd

# calc cl/cd for each airfoil
Cl0, Cd0 = clcdCalc(pressAOA0, 0)
Cl4, Cd4 = clcdCalc(pressAOA4, 4)
Cl6, Cd6 = clcdCalc(pressAOA6, 6)
Cl8, Cd8 = clcdCalc(pressAOA8, 8)
Cl10, Cd10 = clcdCalc(pressAOA10, 10)
Cl12, Cd12 = clcdCalc(pressAOA12, 12)
Cl14, Cd14 = clcdCalc(pressAOA14, 14)
Cls = np.array([Cl0, Cl4, Cl6, Cl8, Cl10, Cl12, Cl14])
Cds = np.array([Cd0, Cd4, Cd6, Cd8, Cd10, Cd12, Cd14])
aoas = np.array([0, 4, 6, 8, 10, 12, 14])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate Cl/Cd plots
plt.figure("CL vs AOA (4)")
plt.plot(aoas, Cls)
plt.xlim(0, 14); plt.xlabel("AOA (deg)")
plt.ylim(0.2, 1); plt.ylabel("$C_L$")
plt.grid(); plt.title("CL vs AOA"); plt.tight_layout()
plt.figure("CD vs AOA (4)")
plt.plot(aoas, Cds)
plt.xlim(0, 14); plt.xlabel("AOA (deg)")
plt.ylim(-0.01, 0.05); plt.ylabel("$C_D$")
plt.grid(); plt.title("CD vs AOA"); plt.tight_layout()
plt.figure("CL vs CD (4)")
plt.plot(Cds, Cls)
plt.xlim(-0.01, 0.05); plt.xlabel("$C_D$")
plt.ylim(0.2, 1); plt.ylabel("$C_L$")
plt.grid(); plt.title("CL vs. CD"); plt.tight_layout()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 5
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 6
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# slope of Cl curve
coefs = np.polyfit(aoas[0:5], Cls[0:5], 1)
plt.figure("Curve Fit CL vs AOA (6)")
plt.plot(aoas, Cls)
plt.plot(aoas[0:5], coefs[0]*aoas[0:5] + coefs[1], "k--")
plt.xlim(0, 14); plt.xlabel("AOA (deg)")
plt.ylim(0.2, 1); plt.ylabel("$C_L$")
plt.title("Curve Fit CL vs AOA (6)")
plt.grid(); plt.tight_layout()
print(f"CL vs AOA slope: {np.rad2deg(coefs[0]):.2f}/rad")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Traverse Data Analysis: Parts 1 & 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get vars for cd0 calc
Ms = Us/np.sqrt(y*R*T)
p0infs = pinf*(1 + (y -1)/2*Ms**2)**(y/(y - 1))

# loop to do numerical integration for cd0
dycs = np.flip(velsAOA0[0][:]/39.27/l)

def traverse(vels, aoa, Cd):
    # determine which U inf and poinf to use
    if aoa == 0:
        j = 0
    elif aoa == 8:
        j = 3
    elif aoa == 12:
        j = 5
        
    # do numerical integration
    Cd0 = 0.
    vels = np.flip(vels[1])
    for k in range(len(vels) - 1):
        Cd0 += (rho*(Us[j]**2 - vels[k]**2))*(dycs[k + 1] - dycs[k])
    
    Cd0 *= 1/(p0infs[j] - pinf)
    
    # calculate viscous drag
    Cvisc = Cd0 - Cd

    Cd_per = (Cd/Cd0)*100
    Cvisc_per = (Cvisc/Cd0)*100
    
    return Cd0, Cvisc, Cd_per, Cvisc_per

# calculate 
Cd0AOA0, CviscAOA0, Cd_perAOA0, Cvisc_perAOA0 = traverse(velsAOA0/3.281, 0, Cd0)
Cd0AOA8, CviscAOA8, Cd_perAOA8, Cvisc_perAOA8 = traverse(velsAOA8/3.281, 8, Cd8)
Cd0AOA12, CviscAOA12, Cd_perAOA12, Cvisc_perAOA12 = traverse(velsAOA12/3.281, 12, Cd12)

# tabulate data
ReCd0d = {"AOA": [0, 8, 12], 
          "Re": np.array([Res[0], Res[1], Res[2]]).round(0), 
          "Cd0": np.array([Cd0AOA0, Cd0AOA8, Cd0AOA12]).round(3),
          "Cdvisc": np.array([CviscAOA0, CviscAOA8, CviscAOA12]).round(3),
          "Cdvisc_per": np.array([Cvisc_perAOA0, Cvisc_perAOA8, Cvisc_perAOA12]).round(1),
          "Cd": np.array([Cd0, Cd8, Cd12]).round(3),
          "Cd_per": np.array([Cd_perAOA0, Cd_perAOA8, Cd_perAOA12]).round(1)}
ReCd0df = pd.DataFrame(ReCd0d)
print(ReCd0df)

# plt.show()