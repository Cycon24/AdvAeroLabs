# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 22:48:46 2024

@author: piewi
"""

import scipy.io as sio
import numpy as np
import os
import matplotlib.pyplot as plt


# Path Containing 4 Folders
path = 'C:\\Users\\piewi\\Documents\\Python Scripts\\Lab 3 Data'

# Folder name of each test case
baseline_path = '\\Baseline(freestream)'
cylinder_path = '\\Cylinder'
aoa0_path = '\\Swept_0aoa'
aoa16_path = '\\Swept_16aoa'

# Empty array to hold each frame data dicitonary
baseline = []
cylinder = []
aoa0 = []
aoa16 = []

# Load in Each Dataframe Into the Respective List
for f in os.listdir(path+baseline_path):
    if f[1] == 'e':
        baseline.append(sio.loadmat(path + baseline_path + '\\' + f))
        
for f in os.listdir(path+cylinder_path):
    if f[1] == 'e':
        cylinder.append(sio.loadmat(path + cylinder_path + '\\' + f))
        
for f in os.listdir(path+aoa0_path):
    if f[1] == 'e':
        aoa0.append(sio.loadmat(path + aoa0_path + '\\' + f))
        
for f in os.listdir(path+aoa16_path):
    if f[1] == 'e':
        aoa16.append(sio.loadmat(path + aoa16_path + '\\' + f))
        
# Obtain measurement positions for each test case
def get_pos(dfs):
    x_pos = dfs[0]['positionX']
    y_pos = dfs[0]['positionY']
    
    X, Y = np.meshgrid(x_pos, y_pos)
    
    return X, Y

b_x_pos, b_y_pos = get_pos(baseline)
c_x_pos, c_y_pos = get_pos(cylinder)
aoa0_x_pos, aoa0_y_pos = get_pos(aoa0)
aoa16_x_pos, aoa16_y_pos = get_pos(aoa16)

# Calculate Average of Each Test Case
def avg_all_fields(dfs):
    # Empty list of velocities in each axis
    vX_lst = []
    vY_lst = []
    
    # Append the velocity matrices for each frame
    for d in dfs:
        vX_lst.append(d['vectorsX'])
        vY_lst.append(d['vectorsY'])
        
    # Convert velocity list to numpy array
    vX = np.array(vX_lst)
    vY = np.array(vY_lst)
    
    # Stack arrays
    all_vX = np.stack(vX)
    all_vY = np.stack(vY)
    
    # Average velocity fields while excluding NaNs
    avg_vX = np.nanmean(all_vX, axis = 0)
    avg_vY = np.nanmean(all_vY, axis = 0)
        
    return avg_vX, avg_vY

b_avg_vX, b_avg_vY = avg_all_fields(baseline)
c_avg_vX, c_avg_vY = avg_all_fields(cylinder)
aoa0_avg_vX, aoa0_avg_vY = avg_all_fields(aoa0)
aoa16_avg_vX, aoa16_avg_vY = avg_all_fields(aoa16)


# Calculate RMS Average of Each Test Case
def RMS_avg(dfs, mean_field_X, mean_field_Y):
    # Empty list of velocities in each axis
    vX_lst = []
    vY_lst = []
    
    # Append the velocity matrices for each frame
    for d in dfs:
        vX_lst.append(d['vectorsX']-mean_field_X)
        vY_lst.append(d['vectorsY']-mean_field_Y)
        
    # Convert velocity list to numpy array
    vX = np.array(vX_lst)
    vY = np.array(vY_lst)
    
    # Stack arrays
    all_vX = np.stack(vX)
    all_vY = np.stack(vY)
    
    # Average velocity fields while excluding NaNs
    rms_avg_vX = np.nanmean(all_vX, axis = 0)
    rms_avg_vY = np.nanmean(all_vY, axis = 0)
        
    return rms_avg_vX, rms_avg_vY

b_rms_vX, b_rms_vY = RMS_avg(baseline, b_avg_vX, b_avg_vY)
c_rms_vX, c_rms_vY = RMS_avg(cylinder, c_avg_vX, c_avg_vY)
aoa0_rms_vX, aoa0_rms_vY = RMS_avg(aoa0, aoa0_avg_vX, aoa0_avg_vY)
aoa16_rms_vX, aoa16_rms_vY = RMS_avg(aoa16, aoa16_avg_vX, aoa16_avg_vY)

# Set the DPI for Plots
plt.rcParams['figure.dpi'] = 500

# Plot Quiver Plots for Each Average
plt.quiver(b_x_pos[:,::2], b_y_pos[:,::2], b_avg_vX[:,::2], b_avg_vY[:,::2], scale=8000,headwidth=5,headlength=4, minlength=0.7, width=0.0009)
plt.title('Baseline Mean Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.quiver(c_x_pos[:,::2], c_y_pos[:,::2], c_avg_vX[:,::2], c_avg_vY[:,::2], scale=8000,headwidth=5,headlength=4, minlength=0.7, width=0.0009)
plt.title('Cylinder Mean Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.quiver(aoa0_x_pos[:,::2], aoa0_y_pos[:,::2], aoa0_avg_vX[:,::2], aoa0_avg_vY[:,::2], scale=8000,headwidth=5,headlength=4, minlength=0.7, width=0.0009)
plt.title('AoA 0° Mean Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.quiver(aoa16_x_pos[:,::2], aoa16_y_pos[:,::2], aoa16_avg_vX[:,::2], aoa16_avg_vY[:,::2], scale=8000,headwidth=5,headlength=4, minlength=0.7, width=0.0009)
plt.title('AoA 16° Mean Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()

# Plot Quiver Plots for Each RMS Avgerage
plt.quiver(b_x_pos[:,::2], b_y_pos[:,::2], b_rms_vX[:,::2], b_rms_vY[:,::2],headwidth=5,headlength=4, minlength=0.7, width=0.0015)
plt.title('Baseline RMS Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.quiver(c_x_pos[:,::2], c_y_pos[:,::2], c_rms_vX[:,::2], c_rms_vY[:,::2],headwidth=5,headlength=4, minlength=0.7, width=0.0015)
plt.title('Cylinder RMS Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.quiver(aoa0_x_pos[:,::2], aoa0_y_pos[:,::2], aoa0_rms_vX[:,::2], aoa0_rms_vY[:,::2],headwidth=5,headlength=4, minlength=0.7, width=0.0015)
plt.title('AoA 0° RMS Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.quiver(aoa16_x_pos[:,::2], aoa16_y_pos[:,::2], aoa16_rms_vX[:,::2], aoa16_rms_vY[:,::2],headwidth=5,headlength=4, minlength=0.7, width=0.0015)
plt.title('AoA 16° RMS Vector Field')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()


# Plot Contour Plots for Each Average X Velocity
plt.contourf(b_x_pos, b_y_pos, b_avg_vX, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('Baseline Horizontal Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(c_x_pos, c_y_pos, c_avg_vX, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('Cylinder Horizontal Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(aoa0_x_pos, aoa0_y_pos, aoa0_avg_vX, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('AoA 0° Horizontal Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(aoa16_x_pos, aoa16_y_pos, aoa16_avg_vX, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('AoA 16° Horizontal Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()  
    
# Plot Contour Plots for Each Average Y Velocity
plt.contourf(b_x_pos, b_y_pos, b_avg_vY, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('Baseline Vertical Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(c_x_pos, c_y_pos, c_avg_vY, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('Cylinder Vertical Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(aoa0_x_pos, aoa0_y_pos, aoa0_avg_vY, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('AoA 0° Vertical Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(aoa16_x_pos, aoa16_y_pos, aoa16_avg_vY, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('AoA 16° Vertical Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show() 
    
# Calculate total velocity for each test case
b_avg_sum = np.sqrt(b_avg_vX**2 + b_avg_vY**2)
c_avg_sum = np.sqrt(c_avg_vX**2 + c_avg_vY**2)
aoa0_avg_sum = np.sqrt(aoa0_avg_vX**2 + aoa0_avg_vY**2)
aoa16_avg_sum = np.sqrt(aoa16_avg_vX**2 + aoa16_avg_vY**2)

# Plot Contour Plots for Each Average Total Velocity
plt.contourf(b_x_pos, b_y_pos, b_avg_sum, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('Baseline Total Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(c_x_pos, c_y_pos, c_avg_sum, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('Cylinder Total Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(aoa0_x_pos, aoa0_y_pos, aoa0_avg_sum, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('AoA 0° Total Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show()
plt.contourf(aoa16_x_pos, aoa16_y_pos, aoa16_avg_sum, 10)
plt.colorbar(orientation='vertical', pad=0.025, shrink=0.8, label='Velocity (mm/s)')
plt.title('AoA 16° Total Velocity')
plt.xlabel('X Position (mm)')
plt.ylabel('Y Position (mm)')
plt.show() 

# Calculate Velocity Profiles
def vel_prof(dfs, dfs_avg, x_index):
    y = dfs[0]['positionY'][0]
    y_offset = y[0]
    y -= y_offset
    vels = dfs_avg[:, x_index]
    return y, vels

# Calculate Upstream Profile
b_vel_prof_y, b_vel_prof_upstream = vel_prof(baseline, b_avg_vX, 0)
c_vel_prof_y, c_vel_prof_upstream = vel_prof(cylinder, c_avg_vX, 0)
aoa0_vel_prof_y, aoa0_vel_prof_upstream = vel_prof(aoa0, aoa0_avg_vX, 0)
aoa16_vel_prof_y, aoa16_vel_prof_upstream = vel_prof(aoa16, aoa16_avg_vX, 0)  

# Plot Upstream Velocity Profile
plt.plot(b_vel_prof_upstream, b_vel_prof_y, label='Baseline')
plt.plot(c_vel_prof_upstream, c_vel_prof_y, label='Cylinder')
plt.plot(aoa0_vel_prof_upstream, aoa0_vel_prof_y, label='AoA 0°')
plt.plot(aoa16_vel_prof_upstream, aoa16_vel_prof_y, label='AoA 16°')
plt.title('Upstream Velocity Profile')
plt.xlabel('Velocity (mm/s)')
plt.ylabel('Height (mm)')
plt.legend(loc='best')
plt.show()

# Last Index of X Array
x_len = 78

# Calculate First Downstream Profile
b_vel_prof_y, b_vel_prof_downstream1 = vel_prof(baseline, b_avg_vX, x_len-10)
c_vel_prof_y, c_vel_prof_downstream1 = vel_prof(cylinder, c_avg_vX, x_len-10)
aoa0_vel_prof_y, aoa0_vel_prof_downstream1 = vel_prof(aoa0, aoa0_avg_vX, x_len-10)
aoa16_vel_prof_y, aoa16_vel_prof_downstream1 = vel_prof(aoa16, aoa16_avg_vX, x_len-10)

# Plot First Downstream Velocity Profile
plt.plot(b_vel_prof_downstream1, b_vel_prof_y, label='Baseline')
plt.plot(c_vel_prof_downstream1, c_vel_prof_y, label='Cylinder')
plt.plot(aoa0_vel_prof_downstream1, aoa0_vel_prof_y, label='AoA 0°')
plt.plot(aoa16_vel_prof_downstream1, aoa16_vel_prof_y, label='AoA 16°')
plt.title('First Downstream Velocity Profile')
plt.xlabel('Velocity (mm/s)')
plt.ylabel('Height (mm)')
plt.legend(loc='best')
plt.show()

# Calculate Last Downstream Profile
b_vel_prof_y, b_vel_prof_downstream2 = vel_prof(baseline, b_avg_vX, x_len-5)
c_vel_prof_y, c_vel_prof_downstream2 = vel_prof(cylinder, c_avg_vX, x_len-5)
aoa0_vel_prof_y, aoa0_vel_prof_downstream2 = vel_prof(aoa0, aoa0_avg_vX, x_len-5)
aoa16_vel_prof_y, aoa16_vel_prof_downstream2 = vel_prof(aoa16, aoa16_avg_vX, x_len-5)

# Plot Last Downstream Velocity Profile
plt.plot(b_vel_prof_downstream2, b_vel_prof_y, label='Baseline')
plt.plot(c_vel_prof_downstream2, c_vel_prof_y, label='Cylinder')
plt.plot(aoa0_vel_prof_downstream2, aoa0_vel_prof_y, label='AoA 0°')
plt.plot(aoa16_vel_prof_downstream2, aoa16_vel_prof_y, label='AoA 16°')
plt.title('Last Downstream Velocity Profile')
plt.xlabel('Velocity (mm/s)')
plt.ylabel('Height (mm)')
plt.legend(loc='best')
plt.show()

# Percent Change in Control Volume
def dy_c(y):
    c = y[-1]
    dy = []
    for i in range(len(y)-1):
        dy.append(y[i+1]-y[i])
    return np.mean(dy/c)

# Get Upstream and Downstream Velocities and Percent Change in CV Height
def cv_vels(dfs, dfs_avg):
    upstream = dfs_avg[:, 0]
    downstream = dfs_avg[:, -1]
    
    y = dfs[0]['positionY'][0]
    y_offset = y[0]
    y -= y_offset
    dy = dy_c(y/1000)

    return dy, upstream/1000, downstream/1000   # Output in m/s

b_cv_dy, b_cv_upstream, b_cv_downstream = cv_vels(baseline, b_avg_vX)
c_cv_dy, c_cv_upstream, c_cv_downstream = cv_vels(cylinder, c_avg_vX)
aoa0_cv_dy, aoa0_cv_upstream, aoa0_cv_downstream = cv_vels(aoa0, aoa0_avg_vX)
aoa16_cv_dy, aoa16_cv_upstream, aoa16_cv_downstream = cv_vels(aoa16, aoa16_avg_vX)

# Calculate Dynamic Pressure for Each Case in SI Units
def q_inf(upstream_v):
    rho = 998.0 #kg/m^3
    u_avg = np.mean(upstream_v) #m/s
    
    return (rho/2)*(u_avg**2) #kg/(m*s^2)

b_q_inf = q_inf(b_cv_upstream)
c_q_inf = q_inf(c_cv_upstream)
aoa0_q_inf = q_inf(aoa0_cv_upstream)
aoa16_q_inf = q_inf(aoa16_cv_upstream)

# Calculate CD using Control Volume Analysis for Each Case
def CV_calc(q, upstream_v, downstream_v, dy):
    rho = 998.0
    
    sums = 0
    for i in range(len(upstream_v)):
        sums += ((upstream_v[i]**2) - (downstream_v[i]**2))*rho*dy
    
    return sums/q

b_CD = CV_calc(b_q_inf, b_cv_upstream, b_cv_downstream, b_cv_dy)
c_CD = CV_calc(c_q_inf, c_cv_upstream, c_cv_downstream, c_cv_dy)
aoa0_CD = CV_calc(aoa0_q_inf, aoa0_cv_upstream, aoa0_cv_downstream, aoa0_cv_dy)
aoa16_CD = CV_calc(aoa16_q_inf, aoa16_cv_upstream, aoa16_cv_downstream, aoa16_cv_dy)

# Print Drag Coefficients for Each Case
print('\nBaseline CD:', b_CD)
print('\nCylinder CD:', c_CD)
print('\nAoA 0 Degrees CD:', aoa0_CD)
print('\nAoA 16 Degrees CD:', aoa16_CD)

# Lengths for Reynold's Numbers
D = 1.5/39.37 #in to m
MAC = 0.64607 #m

# Cylinder Average Upstream Velocity
c_v_avg = np.mean(c_cv_upstream)

# Wing Avgerage Upstream Velocity 
aoa0_v_avg = np.mean(aoa0_cv_upstream)
aoa16_v_avg = np.mean(aoa16_cv_upstream)

# Calculate Reynold's Number in SI Units
def Re_calc(u, L):
    rho = 998.0
    mu = 1.002
    Re = rho*u*L/mu
    return Re

Re_cylinder = Re_calc(c_v_avg, D)
Re_aoa0 = Re_calc(aoa0_v_avg, MAC)
Re_aoa16 = Re_calc(aoa16_v_avg, MAC)

# Display Reynold's Numbers
print('\n----------------------------------------------------------')
print('\nReynold\'s Number for Cylinder:', Re_cylinder)
print('\nReynold\'s Number for Wing (AoA 0 Degrees):', Re_aoa0)
print('\nReynold\'s Number for Wing (AoA 16 Degrees):', Re_aoa16)
