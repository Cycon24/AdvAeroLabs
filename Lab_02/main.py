# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 17:20:32 2024

@author: cycon
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import GasDynamics as GD



# =============================================================================
# Retrieve and organize data
# ============================================================================+
# r'Lab_02/
df_1 = pd.read_csv('AeroLab2_DataSheet.csv').to_numpy()
df_2 = pd.read_csv('AeroLab2_DataSheet_02.csv').to_numpy()

ProbePosition = np.array(df_1[1:,0], dtype=int)     # Integer
ProbePressure = np.array(df_1[1:,1:], dtype=float)  # kPa

BackPressure    = np.array(df_2[0,1:], dtype=float) # psig
InletTemp       = np.array(df_2[1,1:], dtype=float) # C
OutletTemp      = np.array(df_2[2,1:], dtype=float) # C
deltaP          = np.array(df_2[3,1:], dtype=float) # inH2O - Orifice dP
BarometricP     = float(df_2[4,1])                # inHg

# =============================================================================
# Convert Values
# =============================================================================
# Convert into Imperial Units and absolute units (if needed)
# BackPressure                          # Keep in psig
InletTemp =  InletTemp *9/5 + 491.67    # Covnert to R
OutletTemp = OutletTemp*9/5 + 491.67    # Covnert to R
deltaP *= 0.0360912                     # Convert to psi [lbf/in^2]
deltaP_ft2 = deltaP*(12/1)**2           # Convert to [lbf/ft^2]
BarometricP *= 0.4911543448             # Convert t0 psia

# =============================================================================
# Make Calculations
# =============================================================================
# Needed Variables
d_th  = 0.2504              # [in] - Diameter of Nozzle Throat
d_probe = 0.130             # [in] - Diameter of Probe
d_e   = 0.2869              # [in] - Diameter of Nozzle Exit

A_th  = np.pi*((d_th/2)**2 - (d_probe/2)**2)*(1/12)**2   # [ft^2] - Throat Area 
A_e  =  np.pi*((d_e/2)**2  - (d_probe/2)**2)*(1/12)**2    # [ft^2] - Exit Area 

F_a   = 1.0                 # [NonDim] - Thermal Expansion Factor
C_d   = 0.623               # [NonDim] - Discharge Coefficient
Beta  = 1.6255/3.1875       # [NonDim] - d/D - Diameter ratio
g_c   = 32.174              # [𝑙𝑏𝑚−𝑓𝑡/𝑙𝑏𝑓−𝑠^2] - Gravitational Constant

R_univ = 10.731577089016    # [psi⋅ft3/lbmol⋅°R]  - Universal Gas Constant
MW_air = 28.966             # [lbm/lbmol]
R      = R_univ / MW_air    # [psi⋅ft3/lbm⋅°R] Gas Constant if  air

# dP   =

# Calculate density
P_f = (ProbePressure[28,:]+BarometricP) # lbf/in^2
rho_f = P_f / (R*OutletTemp)# [lbm/ft^3] - Fluid Density (needs calculated)

# Calculate Pressure Ratios
Po = 85 # Total pressure was 85 psig
PRs = ProbePressure/Po # TODO: is ProbePressure total? is this eqtn right idfk
PRs_e = ProbePressure[25,:]/Po
PRs_b = ProbePressure[28,:]/Po

# Calculate mass flow rate
paran = (C_d/(np.sqrt(1-Beta**4)))
sqrt  = np.sqrt(2*g_c*rho_f*deltaP_ft2)
mdots = A_th*F_a*paran*sqrt


# =============================================================================
# # Plot Results
plt.clf()
# =============================================================================\

    # TODO: can use locations or probe positions, unclear which he wants
locs = np.arange(-0.5, 2.31, 0.1)
def normPressProfile_single(Pb, Poss, PRs, inParent=False):
    if inParent:
        plt.plot(Poss, PRs, label="Pb={:.0f} psig".format(Pb))
        plt.xlabel("Position (in.)"); plt.ylabel(r"$\dfrac{P}{P_0}$")
    else:
        plt.figure(f"Normalized Pressure Profile (1a): Pb={Pb:.0f} psig")
        plt.title(f"Normalized Pressure Profile, Pb={Pb:.0f} psig")
        plt.plot(Poss, PRs)
        plt.xlim(-0.5, 2.5); plt.ylim(0, 6)
        plt.xlabel("Position (in.)"); plt.ylabel(r"$\dfrac{P}{P_0}$")
        plt.grid()
    
def normPressProfile_all(BackPressureArray, locs, PRsArray):
    plt.figure(f"Normalized Pressure Profile (1a)")
    plt.title(f"Normalized Pressure Profile")
    for i in range(8):
        normPressProfile_single(BackPressureArray[i], locs, PRsArray[:,i],True)
    # Plot Critical Pressure
    plt.plot([min(locs), max(locs)], [0.5283, 0.5283], label="Pb=P_crit")
    plt.legend()
    plt.grid()

# TODO: check this section, unsure but makes plot
plt.figure("Normalized Exit Pressure vs Normalized Back Pressure (1b)")
plt.plot(PRs_e, PRs_b,'*-',label='data')
plt.plot([0,1],[0,1], label='Pe=Pb')
plt.xlabel('Pe/Po')
plt.ylabel('Pb/Po')
plt.legend()
plt.grid()

# can generate multiple plots here (comment out if working on dif. section)
# for i in range(8):
#     normPressProfile(BackPressure[i], locs, PRs[:,i])
normPressProfile_all(BackPressure, locs, PRs)
plt.show()
