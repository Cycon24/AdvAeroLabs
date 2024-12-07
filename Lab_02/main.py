# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 17:20:32 2024

@author: cycon
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import GasDynamics as GD

isSladeorBricWho = True # True for bric, False for swade

# =============================================================================
# Retrieve and organize data
# ============================================================================+
# r'Lab_02/
root = '' if isSladeorBricWho else 'Lab_02/'
df_1 = pd.read_csv(root+'AeroLab2_DataSheet.csv').to_numpy()
df_2 = pd.read_csv(root+'AeroLab2_DataSheet_02.csv').to_numpy()

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
InletTemp =  InletTemp *9/5 + 491.67    # Covnert to R
OutletTemp = OutletTemp*9/5 + 491.67    # Covnert to R
deltaP *= 0.0360912                     # Convert to psi [lbf/in^2]
deltaP_ft2 = deltaP*(12/1)**2           # Convert to [lbf/ft^2]
BarometricP *= 0.4911543448             # Convert t0 psia
BackPressure_abs  = BackPressure + BarometricP # Convert to psia
ProbePressure_abs = ProbePressure+ BarometricP # Convert to psia

# =============================================================================
# Make Calculations
# =============================================================================
# Needed Variables
d_th  = 0.2504/12              # [ft] - Diameter of Nozzle Throat
d_probe = 0.130/12             # [ft] - Diameter of Probe
d_e   = 0.2869/12              # [ft] - Diameter of Nozzle Exit

A_p = np.pi*(d_th/2)**2
A_th  = np.pi*((d_th/2)**2 - (d_probe/2)**2)   # [ft^2] - Throat Area 
A_e  =  np.pi*((d_e/2)**2  - (d_probe/2)**2)   # [ft^2] - Exit Area 
A_orf = np.pi*((1.6255/2)**2)*(1/12)**2        # [ft^2] - Orifice opening Area

F_a   = 1.0                 # [NonDim] - Thermal Expansion Factor
C_d   = 0.623               # [NonDim] - Discharge Coefficient
Beta  = 1.6255/3.1875       # [NonDim] - d/D - Diameter ratio
g_c   = 32.174              # [𝑙𝑏𝑚−𝑓𝑡/𝑙𝑏𝑓−𝑠^2] - Gravitational Constant
R      = 53.3523            # [lbf⋅ft/lbm⋅°R] Gas Constant if  air
R_prime= R/g_c              # [lbf^2⋅s^2 / lbm^2⋅°R] Gas Constant divided by gc


# Calculate density (NOT SURE IF RIGHT)
P_f = (BarometricP+deltaP)*(12/1)**2  # [lbf/ft^2]
rho_f = P_f / (R*OutletTemp)          # [lbm/ft^3] - Fluid Density (Needs verified)

# Calculate Pressure Ratios
Po = 85+BarometricP     # psia, Total pressure was 85 psig 
Po_ft2 = Po*(12/1)**2   # [lbf/ft^2]
PRs =   ProbePressure_abs/Po # TODO: is ProbePressure total? is this eqtn right idfk
PRs_e = ProbePressure_abs[25,:]/Po
PRs_b = ProbePressure_abs[28,:]/Po

# Define important pressure values
# fe = Fully Expanded, ue = Under Expanded
M_ue_sub, M_fe_sup = GD.Mach_at_A(A_e/A_th, Gamma=1.4)
PR_fe = 1/GD.Po_P_ratio(M_fe_sup, Gamma=1.4) # P/Po
PR_ue = 1/GD.Po_P_ratio(M_ue_sub, Gamma=1.4) # P/Po
P_fe = PR_fe * Po
P_ue = PR_ue * Po
Pcrit = 0.5283*Po


# =============================================================================
#   Mass Flow Rates
# =============================================================================
# Experimental mass flow rate
paran = (C_d/(np.sqrt(1-Beta**4)))
sqrt  = np.sqrt(2*g_c*rho_f*deltaP_ft2)
mdots = A_orf*F_a*paran*sqrt # [lbm/s]


# Calculate theoretical mass flow
mdot_the = np.zeros(9)
# Finding it based on the exit pressure aka position 27 (It is same as back pressure, here Pe=Pb)
for i, Ps in enumerate(ProbePressure_abs[26,:]):
    # _ue -> Subsonic isentropic solution
    # _fe -> Supersonic isentropic solution
    if Ps > P_ue:
        # print(i,' - Not Choked {:0.2f}'.format(Ps))
        # Not Choked Flow, find Mach number at exit plane with back pressure
        M = GD.Mach_at_PR(Po/Ps, Gamma=1.4)
        mdot_the[i] = GD.mdot(Po=Po_ft2, To=InletTemp[i], A=A_e,  Mach=M, Gamma=1.4, R=R_prime) # [lbm/s]
    elif Ps > P_fe:
        # choked flow, But shock occurs in nozzle (not ful expand)
        # print(i,' - Choked(ue) {:0.2f}'.format(Ps))
        mdot_the[i] = GD.mdot(Po=Po_ft2, To=InletTemp[i], A=A_th, Mach=1, Gamma=1.4, R=R_prime) # [lbm/s]
    else:
        # Choked and fully expanded flow
        # print(i,' - Choked (fe) {:0.2f}'.format(Ps))
        mdot_the[i] = GD.mdot(Po=Po_ft2, To=InletTemp[i], A=A_th, Mach=1, Gamma=1.4, R=R_prime) # [lbm/s]
    

# =============================================================================
#     Throat Pressures
# =============================================================================
'''
Maybe we can start by finding the theoretical Me based on Pb=Pe,
and the Pb/Po ratio. 
We then find the area ratio Ae/A* from this mach number
Next, we know the area ratio Ae/At, and so we can find the ratio
At/A*, this will give us Mach number at the throat and hence the
pressure ratio at the throat

the = theoretical
th  = throat
star = choked location
'''  

P_th = np.zeros(9)
for i, Pb in enumerate(BackPressure_abs):
    if Pb > P_ue:
        # print(i,' - Not Choked {:0.2f}'.format(Pb)) 
        # Not Choked Flow, find Mach number at exit plane with back pressure, (pe=pb)
        # Subsonic Cases (Assuming Pe=Pb)
        Me_the = GD.Mach_at_PR(Po/Pb, Gamma=1.4)
        A_star = A_e/GD.A_ratio(Me_the, Gamma=1.4) # Theoretical throat area (aka area of throat if flow became choked)

        M_th, M_sup = GD.Mach_at_A(A_th/A_star, Gamma=1.4)
        P_th[i] = Po/GD.Po_P_ratio(M_th, Gamma=1.4)
        
    elif Pb > P_fe:
        # choked flow, But shock occurs in nozzle (not ful expand, pb /= pe)
        # print(i,' - Choked(ue) {:0.2f}'.format(Pb))
        P_th[i] = Po/GD.Po_P_ratio(Mach=1, Gamma=1.4)
    else:
        # Choked and fully expanded flow (pe = pb)
        # print(i,' - Choked (fe) {:0.2f}'.format(Pb))
        P_th[i] = Po/GD.Po_P_ratio(Mach=1, Gamma=1.4)
   
  
# =============================================================================
#   Normal Shock Position
# =============================================================================
'''
Find the position of the shock in relation to the throat position from the experimental
normalized pressure profile graph. Then calculate the area at this shock position using the
equation given below. The area through the nozzle, accounting for the area of the search
tube, is prescribed as: 𝐴(𝑥) = [(0.02246 ⋅ 𝑋) + 0.2504]^2⋅(𝜋/4) − 0.01327 where 𝑋 =
[𝑠ℎ𝑜𝑐𝑘 𝑝𝑜𝑠𝑖𝑡𝑖𝑜𝑛] − [𝑡ℎ𝑟𝑜𝑎𝑡 𝑝𝑜𝑠𝑖𝑡𝑖𝑜𝑛] where X = inches and A = inches2
'''
# Experimental
X_ShockPos_ex = np.array([None, 1.85,1.65,1.25,0.75,None,None,None,None]) # 0 = Nan
X_th = 0.4 # in

A_shock_ex = np.zeros(X_ShockPos_ex.shape)
for i, X_sp in enumerate(X_ShockPos_ex):
    if X_sp == None:
        A_shock_ex[i] = None
    else:
        X_ex = X_sp - X_th
        A_shock_ex[i] = (np.pi/4)*(0.02246*X_ex + 0.2504)**(2) - 0.01327 # in^2


# Theoretical Shock Positions need an iterative solver to find

    
# =============================================================================
# # Plot Results
plt.close('all')
# =============================================================================\

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
        plt.grid(); plt.tight_layout()
    
def normPressProfile_all(BackPressureArray, locs, PRsArray):
    plt.figure(f"Normalized Pressure Profile (1a)")
    plt.title(f"Normalized Pressure Profile")
    for i in range(9):
        normPressProfile_single(BackPressureArray[i], locs, PRsArray[:,i],True)
    # Plot Critical Pressure
    plt.plot([min(locs), max(locs)], [0.5283, 0.5283],'--', label="Pb=P_crit")
    plt.plot([min(locs), max(locs)], [PR_fe, PR_fe], '--',label="Pb=P_sup")
    plt.plot([min(locs), max(locs)], [PR_ue, PR_ue], '--',label="Pb=P_sub")
    plt.plot([locs[5], locs[5]], [0,1], '--k',label='Inlet')
    plt.plot([locs[8], locs[8]], [0,1], '--k',label='Throat')
    plt.plot([locs[25], locs[25]], [0,1], '--k',label='Exit')
    
    handles, labels = plt.gca().get_legend_handles_labels()
    order = []
    for i in range(0,3):
        order.append(i)
        order.append(i+3)
        order.append(i+6)
        order.append(i+9)
        order.append(i+12)
    plt.xlim([-0.5,2.5])
    plt.ylim([-0.25,1.05])    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], ncols=3)
    plt.grid(); plt.tight_layout()


normPressProfile_all(BackPressure, locs, PRs)


plt.figure("Normalized Exit Pressure vs Normalized Back Pressure (1b)")
plt.title("Normalized Exit Pressure vs Normalized Back Pressure")
plt.plot(PRs_e, PRs_b,'*-',label='data')
plt.plot([0,1],[0,1], label='Pe=Pb')
plt.xlabel('Pe/Po')
plt.ylabel('Pb/Po')
plt.legend()
plt.grid(); plt.tight_layout()

plt.figure("Theoretical vs Measured Throat Pressure (1c)")
plt.title("Theoretical vs Measured Throat Pressure")
plt.plot(PRs_b, P_th, "k--", label="Theoretical")
plt.plot(PRs_b, ProbePressure_abs[8], "k-", label="Measured")
plt.xlim(0.1, 1); plt.ylim(0, 100)
plt.xlabel("Pb/Po"); plt.ylabel("Throat Pressure (psia)")
plt.legend(); plt.grid(); plt.tight_layout()


plt.figure("Theoretical vs Measured Mass Flowrate (1e)")
plt.title("Theoretical vs Measured Mass Flowrate")
plt.plot(PRs_b, mdot_the, "k--", label="Theoretical")
plt.plot(PRs_b, mdots, "k-", label="Measured")
plt.axvline(PR_ue, color="tab:blue", linestyle="--", label="Theoretical Critical PR")
plt.axvline(PRs_b[-3], color="tab:blue", linestyle="-", label="Experimental Critical PR")
# plt.xlim(0.1, 1); plt.ylim(0.0007, 0.0015)
plt.xlabel("Pb/Po"); plt.ylabel("Mass Flowrate (lbm/s)")
plt.legend(); plt.grid(); plt.tight_layout()

plt.figure("Area of Sock in Nozzle vs Normalized Back Pressure2")
plt.title("Area of Sock in Nozzle vs Normalized Back Pressure")
plt.plot(PRs_b, A_shock_ex, '*-')
plt.xlim([-.1,1.1])
plt.ylabel('Nozzle Area at Shock Location (in^2)')
plt.xlabel("Pb/Po")
plt.grid()


plt.show()
