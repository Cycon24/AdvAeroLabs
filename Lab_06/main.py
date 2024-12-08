'''
 AEEM 5081L 
 Adv. Aero Lab 6
 Group 11
 
 Supersonic wind tunnel conical shock analysis
 '''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from GasDynamics import Po_P_ratio
plt.close('all')
data_root = 'LabData'

mu = 1.789e-5 # N*s/m^2

# Import CSV data (image data and tunnel correlation)
shock_data = pd.read_csv(data_root+'\\ImageData.csv')

# Import txt data (wind tunnel)
all_files =  os.listdir(data_root)
txt_files = list(filter(lambda x: x[-4:] == '.txt', all_files))  # only text files
tunnel_data = {}
for i,file in enumerate(txt_files):
    if i == 0:
        # Need to do this becuase 2 tests are stored in first txt file
        units = pd.read_csv(data_root+'\\'+file,sep='\t').iloc[0]
        tunnel_data['Test ' + str(i+1)] = pd.read_csv(data_root+'\\'+file,sep='\t').drop(index=0).astype(np.float64) # This data will be chopped
        tunnel_data['Test ' + str(i+2)] = pd.read_csv(data_root+'\\'+file,sep='\t').drop(index=0).astype(np.float64) # This data will be chopped
    else:
        tunnel_data['Test ' + str(i+2)] = pd.read_csv(data_root+'\\'+file,sep='\t').drop(index=0).astype(np.float64) # This data will be chopped

# Get Column names for tunnel data
col_names = tunnel_data['Test 1'].columns.tolist()


# Will need to break Test 1 into 2 parts (at index 400)
tunnel_data['Test 1'] = tunnel_data['Test 1'].iloc[0:400]
tunnel_data['Test 2'] = tunnel_data['Test 2'].iloc[400:]

# Cut data to stable state data
mach_tol = 0.01
Stats = {'Mach mean':[], 'Mach std_dev':[],'Po mean':[], 'Po std_dev':[], 'mdot mean':[], 'mdot std_dev':[]}
for test in tunnel_data:
    # Will use mach number to determine started tunnel, 
    # find max and cut down to before/after the mach number is
    # within a tolerance of the max value
    
    # First get max mach number 
    Mach_array = tunnel_data[test]['Mach Number']#.to_numpy()
    max_M = np.max(Mach_array)
    i_start = None
    i_end = None
    for i,M in enumerate(Mach_array):
        if i_start == None:
            if abs(M-max_M) < mach_tol: # enters max mach range
                i_start = i
        else:
            if abs(M-max_M) > mach_tol: # exits max mach range
                i_end = i 
                break
    tunnel_data[test] = tunnel_data[test].iloc[i_start:i_end]
    
    
    # Calculate mean and standard deviage for Mach, stag pressure, and mass flow
    Stats['Mach mean'].append(np.average(tunnel_data[test]['Mach Number'].to_numpy()))
    Stats['Mach std_dev'].append(np.std(tunnel_data[test]['Mach Number'].to_numpy()))
    
    Stats['Velocity Mean'] = np.average(tunnel_data[test]['Velocity'].to_numpy())
    Stats['Velocity std_dev'] = np.std(tunnel_data[test]['Velocity'].to_numpy())
    
    Stats['Density Mean'] = np.average(tunnel_data[test]['Static Density'].to_numpy())
    Stats['Density std_dev'] = np.std(tunnel_data[test]['Static Density'].to_numpy())
    
    Stats['Po mean'].append(np.average(tunnel_data[test]['Stagnation Pressure'].to_numpy()))
    Stats['Po std_dev'].append(np.std(tunnel_data[test]['Stagnation Pressure'].to_numpy()))
    
    Stats['mdot mean'].append(np.average(tunnel_data[test]['Mass Flow'].to_numpy()))
    Stats['mdot std_dev'].append(np.std(tunnel_data[test]['Mass Flow'].to_numpy()))



# =============================================================================
# NEED TO COMPUTE REYNALDS NUMBER
# =============================================================================
Stats['Re'] = Stats['Density Mean']*Stats['Velocity Mean'] * (6*2.54/100) / mu


# =============================================================================
# Save tables
# =============================================================================
Stats = pd.DataFrame(Stats)
Stats.to_csv('TunnelStats.csv', index=False)

for test in tunnel_data:
    tunnel_data[test].to_csv(test+'_time_data.csv',index=False)

# =============================================================================
# Compute pressure after shock
# =============================================================================
shock_data['Top Pc']= Stats['Po mean']*(1/Po_P_ratio(Stats['Mach mean'])*shock_data['Pc/P1 Top'])
shock_data['Bot Pc']= Stats['Po mean']*(1/Po_P_ratio(Stats['Mach mean'])*shock_data['Pc/P1 Bot'])


'''
Post-Processing Test Data and Shock Wave Analysis
1. Create a table for each test run to report the requested information in the steps below.
2. For each test run, compute the mean and standard deviation of the Mach number, stagnation
pressure, and mass flow rate. Only average the data when the tunnel is started and the Mach
number reaches a threshold near the average Mach number. Compute the Reynolds number
based on tunnel height.

3. Plot an example of the run data in Step 2 for each of the three variables. Label axes. Discuss
the stability of the wind tunnel.

4. e. Plot the results from the test campaign compared to the theoretical shock angles
and pressures based on the cone AoA’s. You can plot on the conical shock chart
(you can digitize the chart or plot it to scale overlaying your data on top) or use the
theoretical data from the calculator.


5. Pick one of the test cases and select a series of images (3-6 frames) showing the startup
sequence and discuss what is happening during the startup sequence. Describe some of the
flow features you observe and point them out (annotating in powerpoint is fine).


'''



# Plot all the properties of tunnel over time for each test
# plt.figure('Mach Numbers')
# =============================================================================
# Mach Numbers vs Time for each test
# =============================================================================
# for test in tunnel_data:
#     t = np.arange(0,len(tunnel_data[test].index)/4,0.25)
#     data = tunnel_data[test]
#     plt.figure('Mach Numbers vs Time')
#     plt.plot(t, data['Mach Number'], label=test)
#     plt.title('Mach Number for Each Test')
#     plt.ylim([1.5,2])
#     plt.grid()
#     plt.legend()

# =============================================================================
# 3. Plot each var for a test over time
# =============================================================================
test = 'Test 11'
t = np.arange(0,len(tunnel_data[test].index)/4,0.25)

fig = plt.figure('Mach Number Test 11')
plt.title('Mach Number')
plt.plot(t, tunnel_data[test]['Mach Number'])
plt.xlabel('Time [s]')
plt.ylabel('Mach Number')
plt.ylim([1.85,1.95])
plt.grid()
plt.savefig('Mach Number Test 11')

fig = plt.figure('Stagnation Pressure Test 11')
plt.title('Stagnation Pressure')
plt.plot(t, tunnel_data[test]['Stagnation Pressure'])
plt.xlabel('Time [s]')
plt.ylabel('Stagnation Pressure [psia]')
plt.ylim([15,30])
plt.grid()
plt.savefig('Stagnation Pressure Test 11')

fig = plt.figure('Mass Flow Rates Test 11')
plt.title('Mass Flow Rates')
plt.plot(t, tunnel_data[test]['Mass Flow'])
plt.xlabel('Time [s]')
plt.ylabel('Mass Flow Rate [kg/s]')
plt.ylim([5,7])
plt.grid()
plt.savefig('Mass Flow Rates Test 11')


# =============================================================================
#  4 Shock angle Plots
# =============================================================================
fig = plt.figure('Shock Angles vs AoA')
shock_data_sorted = shock_data.sort_values(by='AoA')
plt.plot(shock_data_sorted['AoA'], shock_data_sorted['TopShock'], '.-', label='Exp. Shock Angle Top')
plt.plot(shock_data_sorted['AoA'], shock_data_sorted['BotShock'], '.-', label='Exp. Shock Angle Top')
plt.plot(shock_data_sorted['AoA'], shock_data_sorted['CalcTopShock'],'.--', label='Th. Shock Angle Top')
plt.plot(shock_data_sorted['AoA'], shock_data_sorted['CalcBotShock'],'.--', label='Th. Shock Angle Top')
plt.xlabel('Angle of Attack [deg]')
plt.ylabel('Shock Angle [deg]')
plt.legend()
plt.grid()
plt.savefig('Shock Angles vs AoA')

fig = plt.figure('Shock Angles vs Cone Angle')
shock_data_sorted = shock_data.sort_values(by='AoA')
plt.plot(shock_data_sorted['TopAoA'], shock_data_sorted['TopShock'], '.-', label='Exp. Shock Angle Top')
plt.plot(shock_data_sorted['BotAoA'], shock_data_sorted['BotShock'], '.-', label='Exp. Shock Angle Top')
plt.plot(shock_data_sorted['TopAoA'], shock_data_sorted['CalcTopShock'],'.--', label='Th. Shock Angle Top')
plt.plot(shock_data_sorted['BotAoA'], shock_data_sorted['CalcBotShock'],'.--', label='Th. Shock Angle Top')
plt.xlabel('Cone Angle [deg]')
plt.ylabel('Shock Angle [deg]')
plt.legend()
plt.grid()
plt.savefig('Shock Angles vs Cone Angle')

# =============================================================================
# 4 Surface Pressure Plots
# =============================================================================
fig = plt.figure('Surface Pressure vs AoA')
plt.plot(shock_data_sorted['AoA'], shock_data_sorted['Top Pc'], '.-', label='Top')
plt.plot(shock_data_sorted['AoA'], shock_data_sorted['Bot Pc'], '.-', label='Bot')
plt.xlabel('Angle of Attack [deg]')
plt.ylabel('Cone Surface Pressure [psia]')
plt.title('Surface Pressures vs AoA')
plt.legend()
plt.grid()
plt.savefig('Surface Pressure vs AoA')


fig = plt.figure('Surface Pressure vs Cone Angle')
plt.plot(shock_data_sorted['TopAoA'], shock_data_sorted['Top Pc'], '.-', label='Top')
plt.plot(shock_data_sorted['BotAoA'], shock_data_sorted['Bot Pc'], '.-', label='Bot')
plt.xlabel('Relative Cone Angle [deg]')
plt.ylabel('Cone Surface Pressure [psia]')
plt.title('Surface Pressures vs Cone Angle')
plt.legend()
plt.grid()
plt.savefig('Surface Pressure vs Cone Angle')