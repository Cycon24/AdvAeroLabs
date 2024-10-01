# Adv. Aero Lab 1 Plotting
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import data from files
presDF = pd.read_csv("Lab_01/NACA_4412_Data.csv")
travDF = pd.read_csv("Lab_01/traverse.csv")

# converting to numpy arrays

# [0][:] returns traverse locations, [1][:] returns velocity (ft/s)
velsAOA0 = np.array([travDF["AOA"][1:], travDF["0"][1:]], dtype="float64")
velsAOA8 = np.array([travDF["AOA"][1:], travDF["8"][1:]], dtype="float64")
velsAOA12 = np.array([travDF["AOA"][1:], travDF["12"][1:]], dtype="float64")