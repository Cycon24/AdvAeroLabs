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
df_1 = pd.read_csv('AeroLab2_DataSheet.csv').to_numpy()
df_2 = pd.read_csv('AeroLab2_DataSheet_02.csv').to_numpy()

ProbePosition = np.array(df_1[1:,0], dtype=int)
PropePressure = np.array(df_1[1:,1:], dtype=float)

BackPressure    = np.array(df_2[0,1:], dtype=int) # psig
InletTemp       = np.array(df_2[1,1:], dtype=int) # C
OutletTemp      = np.array(df_2[2,1:], dtype=int) # C
deltaP          = np.array(df_2[3,1:], dtype=int) # inH2O - Orifice dP
BarometricP     = float(df_2[4,1])                # inHg

# =============================================================================
# Convert Values
# =============================================================================




# =============================================================================
# Make Calculations
# =============================================================================




# =============================================================================
# # Plot Results
# =============================================================================

