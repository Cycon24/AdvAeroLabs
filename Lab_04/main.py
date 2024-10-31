# Adv. Aero Lab 4 Plotting

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from sklearn.linear_model import LinearRegression

# Set constants
R = 287             # R of air (SI)
T = 297.594         # tunnel temp (K)
y = 1.4             # gamma
mu = 1.73e-5        # viscosity of air (SI)
pinf = 99187.324    # barometric pressure (Pa)
rho = pinf/R/T      # density (kg/m^3)
Sref = 64/1550      # wing area (m^2)
MAC = 4/3.281       # mean aerodynamic chord (m)
Fusl = 10.1/3.281   # fuselage length (m)
AR = 2.             # aspect ratio

# Import data from file
cwd = os.getcwd()
data_loc = str(cwd + '/Lab_04/Lab4_AllData.csv')
alldataDF = pd.read_csv(data_loc)
alldataDF.columns = alldataDF.columns.str.strip()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 1 - Pre-processing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Average all data points for each case
casesDF = pd.DataFrame(alldataDF.groupby('Case ID', as_index=False).mean())

# Switch signs for normal force, side force, and pitching moment
casesDF['Normal Force'] = casesDF['Normal Force']*-1
casesDF['Side Force'] = casesDF['Side Force']*-1
casesDF['Pitching Moment'] = casesDF['Pitching Moment']*-1

# Convert all data values to SI
casesDF['Side Force'] = casesDF['Side Force']/4.448                 # lbf to N
casesDF['Pitching Moment'] = casesDF['Pitching Moment']*0.119825    # lbf*in to N*m
casesDF['Normal Force'] = casesDF['Normal Force']/4.448             # lbf to N
casesDF['Yawing Moment'] = casesDF['Yawing Moment']*0.119825        # lbf*in to N*m
casesDF['Rolling Moment'] = casesDF['Rolling Moment']*0.119825      # lbf*in to N*m
casesDF['Axial Force'] = casesDF['Axial Force']/4.448               # lbf to N
casesDF['Dynamic Pressure'] = casesDF['Dynamic Pressure']*6894.76   # psi to Pa
casesDF['Static Pressure'] = casesDF['Static Pressure']*6894.76     # psi to Pa
casesDF['Air Speed'] = casesDF['Air Speed']/3.281                   # ft/s to m/s

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 2 - Parasitic Drag
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Velocity sweeps - 0 AoA, cases 1-5

Cd0s = []
ReMACs = []
ReFUSEs = []
Aspeeds = []
CLs = []
CYs = []

for i in range(5):
    D = casesDF['Axial Force'].iloc[i]
    L = casesDF['Normal Force'].iloc[i]
    Y = casesDF['Side Force'].iloc[i]
    q = casesDF['Dynamic Pressure'].iloc[i]
    Aspeed = np.sqrt((2*q)/rho)
    Aspeeds.append(Aspeed)

    # Calculate parasitic drag coefficient and Reynold's number (2a)
    Cd0s.append(D/(q*Sref)) # for AoA = 0, Cd0 = CD
    ReMACs.append((rho*Aspeed*MAC)/mu)
    ReFUSEs.append((rho*Aspeed*Fusl)/mu)

    # Calculate lift and side force coefficients (2b)
    CL = L/(q*Sref)
    CY = Y/(q*Sref)
    CLs.append(CL)
    CYs.append(CY)

# Tabulate calculated values with data to verify 0 AoA and yaw angle (2c)
ReCs = {"Case": [1, 2, 3, 4, 5], 
        "Airspeed, m/s": Aspeeds,
        "ReMAC": ReMACs, 
        "ReFus": ReFUSEs,
        "Axial Force, N": casesDF['Axial Force'][0:5].tolist(),
        "Normal Force, N": casesDF['Normal Force'][0:5].tolist(),
        "Side Force, N": casesDF['Side Force'][0:5].tolist(),
        "Cd0": Cd0s,
        "CL": CLs,
        "CY": CYs
        }

ReCsDF = pd.DataFrame(ReCs)
print('\nPart 2c: \n', ReCsDF)

# Generate plots for Cd0 vs ReMAC and Cd0 vs ReFus (2d)
plt.figure('Cd0 vs ReMAC')
plt.plot(ReMACs, Cd0s)
plt.xlabel('$Re_{MAC}$'); plt.xlim(0,7000000)
plt.ylabel('$C_{d0}$'); plt.ylim(0.0006,0.0012)
plt.title('$C_{d0}$ vs $Re_{MAC}$')
plt.grid(); plt.tight_layout()

plt.figure('Cd0 vs ReFus')
plt.plot(ReFUSEs, Cd0s)
plt.xlabel('$Re_{fus}$'); plt.xlim(0,16000000)
plt.ylabel('$C_{d0}$'); plt.ylim(0.0006,0.0012)
plt.title('$C_{d0}$ vs $Re_{fus}$')
plt.grid(); plt.tight_layout()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Data Analysis: Part 3 - Lift, Drag and Pitching Moment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# AoA sweeps - constant velocity, cases 3 and 6-10
indexList = [2, 5, 6, 7, 8, 9]
AoAs = [0., 4., 6., 8., 10., 12.]

CLs = []
CDs = []
CMs = []
Cd0s = []

AFs = []
NFs = []
PMs = []

j = 0
for i in range(11):

    if i in indexList:
        
        a = np.deg2rad(AoAs[j])
        j += 1

        AF = casesDF['Axial Force'].iloc[i]; AFs.append(AF)
        NF = casesDF['Normal Force'].iloc[i]; NFs.append(NF)
        PM = casesDF['Pitching Moment'].iloc[i]; PMs.append(PM)
        q = casesDF['Dynamic Pressure'].iloc[i]

        L = NF*np.cos(a) - AF*np.sin(a)
        D = AF*np.cos(a) + NF*np.sin(a)

        # Calculate CL, CD and CM (3a)
        CL = L/(q*Sref)
        CD = D/(q*Sref)
        CM = PM/(q*Sref*MAC)
        CLs.append(CL)
        CDs.append(CD)
        CMs.append(CM)

        # Cd0s.append(D - (1/(np.pi*AR*0.8))*CL**2)

    else:
        i += 1

# Tabulate calculated values with data (3b)
AoACs = {"Case": [3, 6, 7, 8, 9, 10], 
         "AoA, degrees": AoAs,
         "Axial Force, N": AFs,
         "Normal Force, N": NFs,
         "Pitching Moment, N*m": PMs,
         "CL": CLs,
         "CD": CDs,
         "CM": CMs
        }
AoACsDF = pd.DataFrame(AoACs)
print('\nPart 3b: \n', AoACsDF)

# Plot CL vs. AoA, plus generate lift curve slope (3c)
coefs = np.polyfit(AoAs[0:3], CLs[0:3], 1)
slope = []
for i in range(len(AoAs)):
    slope.append(AoAs[i]*coefs[0] + coefs[1])

plt.figure("CL vs AOA")
plt.plot(AoAs, CLs, label="$C_L$")
plt.plot(AoAs[0:3], slope[0:3], "k--", label="Lift Curve Slope = 0.0045/deg")
plt.xlim(0, 12); plt.xlabel("AoA (deg)")
plt.ylim(0,0.04); plt.ylabel("$C_L$")
plt.title("$C_L$ vs AoA")
plt.legend(loc='lower right')
plt.grid(); plt.tight_layout()

print(f"\nPart 3c: CL vs AOA slope: {(coefs[0]):.4f}/deg or {np.rad2deg(coefs[0]):.4f}/rad")

# Plot CL vs CD (3d)
plt.figure("CL vs CD")
plt.plot(CDs, CLs)
plt.xlim(0,0.009); plt.xlabel("$C_D$")
plt.ylim(0,0.04); plt.ylabel("$C_L$")
plt.title("$C_L$ vs $C_D$")
plt.grid(); plt.tight_layout()

# Plot CL vs CD and estimate K (3e)
Cd0 = CDs[0] # Assuming Cd0 from AoA = 0 case

# CD - Cd0
CDs_est = np.array([CDs - Cd0]).reshape(-1,1)
CDs_est = CDs_est[:-2]  # Removing AoAs 10 and 12 deg (stall)

# CL^2
CLs_copy = CLs.copy()
CLs_copy[0] = 0.
CLs_sq = np.array([CLs])**2
CLs_sq = CLs_sq.reshape(-1,1)
CLs_sq = CLs_sq[:-2]    # Removing AoAs 10 and 12 deg (stall)

# Carry out linear regression to determine K
model = LinearRegression(fit_intercept=False)
model.fit(CLs_sq, CDs_est)
K = model.coef_[0]

print(f'\nPart 3e: Estimated drag due to lift factor, K = {K[0]:.2f}')

# Re-calculate CLs across span of CD values with new K
CDpred = np.linspace(np.array(CDs[:-2]).min(), np.array(CDs[:-2]).max(), 100)
CLpred = []
for i in range(len(CDpred)):
    CLpred.append(np.sqrt((CDpred[i] - Cd0) / K)) 

plt.figure("Drag due to lift factor, K")
plt.plot(CDs, CLs, label="$C_L$ vs $C_D$ curve")
plt.plot(CDpred, CLpred, label="Fitted Curve, K = 1.40")
plt.xlim(0,0.009); plt.xlabel("$C_D$")
plt.ylim(0,0.04); plt.ylabel("$C_L$")
plt.title("$C_L$ vs $C_D$ with K estimation")
plt.legend(loc='lower right')
plt.grid(); plt.tight_layout()

# Plot CM vs AoA and estimate the pitching moment derivative (3f)
CMa = (CMs[-3] - CMs[0])/(AoAs[-3] - AoAs[0])  # Removing AoAs 10 and 12 deg (stall)
print(f'\nPart 3f: Pitching Moment Derivative, CMa = {CMa:.6f}\n')

plt.figure("CM vs AOA")
plt.plot(AoAs, CMs, label="$C_M$")
plt.plot([0., 8.,], [CMs[0], CMs[-3]], "k--", label="Pitching Moment Derivative = 0.000861/deg")
plt.xlim(0, 12); plt.xlabel("AoA (deg)")
plt.ylim(0,0.008); plt.ylabel("$C_M$")
plt.title("$C_M$ vs AoA")
plt.legend(loc='lower right')
plt.grid(); plt.tight_layout()

# plt.show()


