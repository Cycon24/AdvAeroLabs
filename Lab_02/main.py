# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 17:20:32 2024

@author: cycon
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# =============================================================================
# Retrieve and organize data
# =============================================================================
df_1 = pd.read_csv(r'Lab_02/AeroLab2_DataSheet.csv').to_numpy()
df_2 = pd.read_csv(r'Lab_02/AeroLab2_DataSheet_02.csv').to_numpy()

ProbePosition = np.array(df_1[1:,0], dtype=int)
ProbePressure = np.array(df_1[1:,1:], dtype=float)

BackPressure    = np.array(df_2[0,1:], dtype=float) # psig
InletTemp       = np.array(df_2[1,1:], dtype=float) # C
OutletTemp      = np.array(df_2[2,1:], dtype=float) # C
deltaP          = np.array(df_2[3,1:], dtype=float) # inH2O - Orifice dP
BarometricP     = float(df_2[4,1])                # inHg

# =============================================================================
# Convert Values
# =============================================================================
# Convert into SI Units and absolute units (if needed)



# =============================================================================
# Make Calculations
# =============================================================================
# Needed Variables
A_th  = 0.# [m^2] - Throat Area 
F_a   = 1.                  # - Thermal Expansion Factor
C_d   = 0.623               # - Discharge Coefficient
Beta  = 1.6255/3.1875       # d/D - Diameter ratio
g_c   = 1.                  # [kg*m/N*s^2] - Gravitational Constant
# rho_f = # [kg/m^3] - Fluid Density
# dP   =
Prats = BarometricP/ProbePressure # TODO: is ProbePressure total? is this eqtn right idfk


# =============================================================================
# # Plot Results
# =============================================================================\
# TODO: can use locations or probe positions, unclear which he wants
locs = np.arange(-0.5, 2.31, 0.1)
def normPressProfile(Pb, Poss, Prats):
    plt.figure(f"Normalized Pressure Profile (1a): Pb={Pb:.0f} psig")
    plt.title(f"Normalized Pressure Profile, Pb={Pb:.0f} psig")
    plt.plot(Poss, Prats)
    plt.xlim(-0.5, 2.5); plt.ylim(0, 6)
    plt.xlabel("Position (in.)"); plt.ylabel(r"$\dfrac{P}{P_0}$")
    plt.grid()


# TODO: check this section, unsure but makes plot
plt.figure("Normalized Exit Pressure vs Normalized Back Pressure (1b)")
erats = ProbePressure[25,:]/ProbePressure[0,:]
brats = BackPressure/ProbePressure[0,:]
plt.plot(erats, brats)
plt.grid()

# can generate multiple plots here (comment out if working on dif. section)
for i in range(8):
    normPressProfile(BackPressure[i], locs, Prats[:,i])

plt.show()
