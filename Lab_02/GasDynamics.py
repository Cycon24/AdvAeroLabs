# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:34:05 2023

@author: cycon
"""
# Imports
import numpy as np
import matplotlib.pyplot as plt
import math

# import _tools.RootAlgorithms as rt
# import _tools.Interpolator as intrp

# =============================================================================
# From root search  imports (copy-pasted so dont have to import shit)
# =============================================================================
def rootsearch(f,a,b,dx):
    '''
    x1,x2 = rootsearch(f,a,b,dx).
    Searches the interval (a,b) in increments dx for
    the bounds (x1,x2) of the smallest root of f(x).
    Returns x1 = x2 = None if no roots were detected.
    '''
    x1 = a; f1 = f(a)
    x2 = a + dx; f2 = f(x2)
    while np.sign(f1) == np.sign(f2):
        if x1 >= b: 
            print('Root Search not found between a={} b={}, with dx={}'.format(a, b, dx))
            return None, None
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2)
    return x1,x2

def newtonRaphson(f,df,a,b,tol=1.0e-9):
    fa = f(a)
    if fa == 0.0: return a
    fb = f(b)
    if fb == 0.0: return b
    if np.sign(fa) == np.sign(fb): 
        print('Root is not bracketed')
        return None
    
    x = 0.5*(a + b)
    
    for i in range(30):
        fx = f(x)
        if fx == 0.0: return x
        
        # Tighten the brackets on the root
        if np.sign(fa) != np.sign(fx): b = x
        else: a = x
        
        # Try a Newton-Raphson step
        dfx = df(x)
        
        # If division by zero, push x out of bounds
        try: dx = -fx/dfx
        except ZeroDivisionError: dx = b - a
        x = x + dx
        
        # If the result is outside the brackets, use bisection
        if (b - x)*(x - a) < 0.0:
            dx = 0.5*(b - a)
            x = a + dx
            
        # Check for convergence
        if abs(dx) < tol*max(abs(b),1.0): return x
    print('Too many iterations in Newton-Raphson')
    
    
# =============================================================================
#               Isentropic Exapnsion/Compression
# =============================================================================
# To/T
def To_T_ratio(Mach,Gamma=1.4):
    return 1 + (Gamma-1)*(Mach**2)/2

# Po/P
def Po_P_ratio(Mach, Gamma =1.4):
    return (1 + (Gamma-1)*(Mach**2)/2)**(Gamma/(Gamma-1))

# Po/P
def Po_P_ratio2(To, T, Gamma=1.4):
    return (To/T)**(Gamma/(Gamma-1))

def RHO_o_Rho(Mach, Gamma=1.4):
    return (1 + (Gamma-1)*(Mach**2)/2)**(1/(Gamma-1))


def mdot(Po, To, A, Mach=1, Gamma=1.4, R=287):
    first = (Po*A/(np.sqrt(R*To))) *Mach*np.sqrt(Gamma)
    second = 1 + ((Gamma-1)/2)*(Mach**2)
    power = (1+Gamma)/(2*(1-Gamma))
    return first*np.power(second, power)


# =============================================================================
#            Normal Shock Relations Properties 
# =============================================================================
# velocity ratio across normal shock: v2/v1
def vR_n(Mach1, Mach2, Temp1, Temp2):
    # Temperatures are static, in absolute units
    return (Mach2/Mach1)*np.sqrt(Temp2/Temp1)

# density ratio across normal shock: rho2/rho1
def rhoR_n(Mach1, Mach2, Temp1, Temp2):
    # Temperatures are static, in absolute units
    return (Mach1/Mach2)*np.sqrt(Temp1/Temp2)

# static temperature ratio across normal shock: T2/T1
def TR_n(Mach1, Mach2, Gamma=1.4):
    return (1+(Mach1**2)*(Gamma-1)/2)/(1+(Mach2**2)*(Gamma-1)/2)

# static pressure ratio across normal shock: p2/p1
def PR_n(Mach1, Mach2, Gamma=1.4):
    return (1+Gamma*Mach1**2)/(1+Gamma*Mach2**2)

# mach number after a shock wave
def Mach2_n(Mach1, Gamma=1.4):
    num = Mach1**2 + 2/(Gamma-1)
    den = -1 + (2*Gamma/(Gamma-1))*Mach1**2   
    return np.sqrt(num/den)

# Redefine equations all based on M1 and Prop1:
def T2_T1_n(Mach1, Gamma=1.4):
    Mach2 = Mach2_n(Mach1, Gamma)
    return TR_n(Mach1, Mach2, Gamma)

def P2_P1_n(Mach1, Gamma=1.4):
    Mach2 = Mach2_n(Mach1, Gamma)
    return PR_n(Mach1, Mach2, Gamma)

def Po2_Po1_n(Mach1, Gamma=1.4):
    p2_p1 = P2_P1_n(Mach1, Gamma)
    T2_T1 = T2_T1_n(Mach1, Gamma)
    return p2_p1*(1/T2_T1)**(Gamma/(Gamma-1))

# ---- Stagnation Properties ----
# stagnation pressure ratio across normal shock: po2/po1
def stag_PR_n(Pres1, Pres2, Temp1, Temp2, Gamma=1.4):
    return (Pres2/Pres1)*(Temp1/Temp2)**(Gamma/(Gamma-1))

# stagnation density ratio across normal shock: ρo2/ρo1
def stag_rhoR_n(stagPres1, stagPres2):
    return stagPres1/stagPres2

# throat area across normal shock? A2*/A1*
def throatAreaR_n(stagPres1, stagPres2):
    return stagPres1/stagPres2

#_________________________________________________________



# =============================================================================
#                Nozzle Relations (Isentropic Flow)
# =============================================================================

def Mach_at_A(Ai_At, Gamma=1.4):
    '''
    Calcultes the MachNumberat the area ratio A/A*
    Returns subsonic and supersonic results

    Parameters
    ----------
    Ai_At : Float
        Area ratio: A/A*.
    Gamma : Float, optional
        Gas Dependent. The default is 1.4 (air).
    AutoProceed: Bool, optional
        Toggles if the function will automatically proceed if a solution is not found.
        I

    Returns
    -------
    M_subsonic, M_supersonic : Float
        Subsonic and Supersonic Mach number solution.

    '''
    if Ai_At < 1: raise TypeError('Area Ratio MUST be >= 1')
    
    def FA_ratio(M):
        return -Ai_At + (1/M)*((1 + M**2*(Gamma-1)/2)/((Gamma+1)/2))**((Gamma+1)/(2*Gamma-2))
    
    def dAR_M(M):
        first_num = - (2 + (M**2)*(Gamma-1))**((Gamma+1)/(2*Gamma-2))
        first_den = (M**2)*(Gamma+1)**((Gamma+1)/(2*Gamma-2))
        second = ((2 + (M**2)*(Gamma-1))/(Gamma+1))**((-Gamma+3)/(2*Gamma-2))
        return first_num/first_den + second
    
    # Solve for first root (subsonic solution)
    low_r, high_r = rootsearch(FA_ratio, 0.001, 8, 0.01)
    if (low_r == None or high_r == None):
        raise TypeError('Root Search Solution not found')
    Msub = newtonRaphson(FA_ratio,dAR_M,low_r,high_r)
    
    # Solve for second root (sonic solution)
    low_r, high_r = rootsearch(FA_ratio, Msub+0.01, 15, 0.01)
    Msup = newtonRaphson(FA_ratio,dAR_M, low_r,high_r)
    
    return Msub, Msup

# Area ratio: A/A*
def A_ratio(M, Gamma=1.4):
    '''
    Calculates the ratio of Area / Throat Area from Mach number and Gamma

    Parameters
    ----------
    M : Float
        Mach Number.
    Gamma : Float, optional
        The default is 1.4.

    Returns
    -------
    A_At : Float
        A/A*: Area ratio.

    '''
    A_At = (1/M)*((1 + M**2*(Gamma-1)/2)/((Gamma+1)/2))**((Gamma+1)/(2*Gamma-2))
    return A_At

def A_throat(m_dot, Po, To, Gamma=1.4, R=287.1):
    '''
    Calculates the Throat Area [m^2] for choked flow.

    Parameters
    ----------
    m_dot : Float
        Mass flow rate in kg/s.
    Po : Float
        Total Pressure in Pa.
    To : Float
        Total Temperature in K.
    Gamma : Float, optional
        Gas dependent. The default is 1.4 (air).
    R : Float, optional
        Gas Constant J/kg*K. The default is 287.1 (air).

    Returns
    -------
    Float
        Area of the throat of choked flow.

    '''
    return m_dot/((Po/np.sqrt(R*To))*np.sqrt(Gamma)*np.power(1+(Gamma-1)/2, (Gamma+1)/(2-2*Gamma)))


def Mach_at_PR(Po_P, Gamma=1.4):
    '''
    Calculates the Mach number of isentropic flow
    At a pressure ratio of Po/P

    Parameters
    ----------
    Po_P : Float
        Pressure ratio.
    Gamma : Float, optional
        Gas dependent. The default is 1.4 (air).

    Returns
    -------
    Float
        Mach number at the point where P = Po / (Po/P).

    '''
    return np.sqrt(((Po_P)**((Gamma-1)/Gamma) - 1)*2/(Gamma-1))

def A_At_shock(Pb, Po, Ae_At, Gamma=1.4, R=287.1):
    return None

# __________________________________________


# =============================================================================
#           Main 
# =============================================================================

if __name__ == "__main__":
    print('bing bong')
    print(Mach_at_A(1.13122))
    