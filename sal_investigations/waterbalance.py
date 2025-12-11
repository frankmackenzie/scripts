#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:52:28 2024

@author: mackenfr
"""

import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import os

#%% paths
path_to_data = "/Volumes/SanDisk/" 

#%%
lig_wais_file = path_to_data + "lig_wais_lig/280/tavg_cat.waterbal.9000.nc"
lig_wais = xr.open_dataset(lig_wais_file, decode_times=False)

pi_wais_file = path_to_data + "pi_wais_lig/280/tavg_cat.waterbal.9000.nc"
pi_wais = xr.open_dataset(pi_wais_file, decode_times=False)


pi_control_file = path_to_data + "picontrol/280/tavg_cat.waterbal.9000.nc"#"/Users/mackenfr/data/uvic_outputs/pi_wais_lig_wais/pi_wais_lig_wais_tavg.09000.01.01.nc"
pi_control = xr.open_dataset(pi_control_file, decode_times=False)

#if average from dif no of timesteps
#lig_wais, pi_control = xr.alig_waisn(lig_wais, pi_control, join='override')

lon = lig_wais["lon"]
lat = lig_wais["lat"]


#%%  global  water balance

for data in [lig_wais, pi_wais, pi_control]:
    precip = data["F_precip"][-1,:,:].sum()
    evap = data["F_evap"][-1,:,:].sum()
    water_bal = precip - evap
    print(water_bal)

#%% water bal over time

times = lig_wais["time"]
lig_wais_waterbal = []
for i in list(range(len(lig_wais["time"]))):
    precip = lig_wais["F_precip"][i,:,:].sum()
    evap = lig_wais["F_evap"][i,:,:].sum()
    water_bal = precip - evap
    lig_wais_waterbal.append(float(water_bal))
    
pi_wais_waterbal = []
for i in list(range(len(lig_wais["time"]))):
    precip = pi_wais["F_precip"][i,:,:].sum()
    evap = pi_wais["F_evap"][i,:,:].sum()
    water_bal = precip - evap
    pi_wais_waterbal.append(float(water_bal))
    

pi_control_waterbal = []
for i in list(range(len(lig_wais["time"]))):
    precip = pi_control["F_precip"][i,:,:].sum()
    evap = pi_control["F_evap"][i,:,:].sum()
    water_bal = precip - evap
    pi_control_waterbal.append(float(water_bal))
    
    
fig, ax = plt.subplots()
#ax2 = ax.twinx()

ax.plot(times, lig_wais_waterbal, label="lig_wais_lig")
ax.plot(times, pi_wais_waterbal, label="pi_wais_lig")
ax.plot(times, pi_control_waterbal, label="pi_control")
ax.set_title("total water balance over time")
ax.legend()

#%% zonal mean water bal??? for last time step


lats = lig_wais["lat"]
lig_wais_waterbal = []
for i in list(range(len(lig_wais["lat"]))):
    precip = lig_wais["F_precip"][-1,i,:].mean()
    evap = lig_wais["F_evap"][-1,i,:].mean()
    water_bal = precip - evap
    lig_wais_waterbal.append(float(water_bal))
    
pi_wais_waterbal = []
for i in list(range(len(lig_wais["lat"]))):
    precip = pi_wais["F_precip"][-1,i,:].mean()
    evap = pi_wais["F_evap"][-1,i,:].mean()
    water_bal = precip - evap
    pi_wais_waterbal.append(float(water_bal))
    

pi_control_waterbal = []
for i in list(range(len(lig_wais["lat"]))):
    precip = pi_control["F_precip"][-1,i,:].mean()
    evap = pi_control["F_evap"][-1,i,:].mean()
    water_bal = precip - evap
    pi_control_waterbal.append(float(water_bal))
    
fig, ax = plt.subplots()

ax.plot(lats, pi_control_waterbal, label="pi_control")
ax.plot(lats, lig_wais_waterbal, label="lig_wais_lig")
ax.plot(lats, pi_wais_waterbal, label="pi_wais_lig")


ax.set_title("Zonal mean water balance")
ax.legend()

