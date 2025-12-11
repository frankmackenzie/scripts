#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 17:28:56 2024

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
lig_wais_file = path_to_data + "lig_wais_lig/280/tavg_cat.rivdis.sal.9000.nc"
lig_wais = xr.open_dataset(lig_wais_file, decode_times=False)

pi_wais_file = path_to_data + "pi_wais_lig/280/tavg_cat.rivdis.sal.9000.nc"
pi_wais = xr.open_dataset(pi_wais_file, decode_times=False)


picontrol_file = path_to_data + "picontrol/280/tavg_cat.rivdis.sal.9000.nc"#"/Users/mackenfr/data/uvic_outputs/pi_wais_lig_wais/pi_wais_lig_wais_tavg.09000.01.01.nc"
picontrol = xr.open_dataset(picontrol_file, decode_times=False)

#if average from dif no of timesteps
#lig_wais, picontrol = xr.alig_waisn(lig_wais, picontrol, join='override')

lon = lig_wais["lon"]
lat = lig_wais["lat"]


#%%  global river discharge  and salinity

#river discharge
lig_wais_timesteps = lig_wais["time"]
lig_wais_rivdis = []
for i in list(range(len(lig_wais["time"]))):
    lig_wais_riv = lig_wais['F_rivdis'][i, :, :].sum()
    lig_wais_rivdis.append(float(lig_wais_riv))
    
pi_wais_rivdis = []
for i in list(range(len(pi_wais["time"]))):
    pi_wais_riv = pi_wais['F_rivdis'][i, :, :].sum()
    pi_wais_rivdis.append(float(pi_wais_riv))

pi_timesteps = picontrol["time"]
pi_control_rivdis = []
for i in list(range(len(picontrol["time"]))):
    pi_control_riv = picontrol['F_rivdis'][i, :, :].sum()
    pi_control_rivdis.append(float(pi_control_riv))

#salinity
lig_wais_salinity = []
for i in list(range(len(lig_wais["time"]))):
    lig_wais_sal = lig_wais['O_sal'][i, :, :, :].sum()
    lig_wais_salinity.append(float(lig_wais_sal))
    
pi_wais_salinity=[]
for i in list(range(len(pi_wais["time"]))):
    pi_wais_sal = pi_wais['O_sal'][i, :, :, :].sum()
    pi_wais_salinity.append(float(pi_wais_sal))

pi_timesteps = picontrol["time"]
pi_control_salinity=[]
for i in list(range(len(picontrol["time"]))):
    pi_control_sal = picontrol['O_sal'][i, :, :, :].sum()
    pi_control_salinity.append(float(pi_control_sal))

fig, ax = plt.subplots()

ax.plot(pi_timesteps, pi_control_rivdis, label="picontrol")
ax.plot(lig_wais_timesteps, lig_wais_rivdis, label="lig_wais_lig")
ax.plot(lig_wais_timesteps, pi_wais_rivdis, label="pi_wais_lig")
ax.set_ylabel("sum rivdis")
ax.legend()

fig, ax = plt.subplots()

ax.plot(pi_timesteps, pi_control_rivdis, label="picontrol")
ax.plot(lig_wais_timesteps, lig_wais_rivdis, label="lig_wais_lig")
ax.plot(lig_wais_timesteps, pi_wais_rivdis, label="pi_wais_lig")
ax.set_ylabel("sum rivdis")
ax.legend()
ax.set_title("Global")

fig, ax2 = plt.subplots()
ax2.plot(pi_timesteps, pi_control_salinity, label="picontrol")
ax2.plot(lig_wais_timesteps, lig_wais_salinity, label="lig_wais_lig")
ax2.plot(lig_wais_timesteps, pi_wais_salinity, label="pi_wais_lig")
ax2.set_ylabel("sum salinity")

ax2.legend()
ax2.set_title("Global")
ax.set_title("Global")

fig, ax2 = plt.subplots()
ax2.plot(pi_timesteps, pi_control_salinity, label="picontrol")
ax2.plot(lig_wais_timesteps, lig_wais_salinity, label="lig_wais_lig")
ax2.plot(lig_wais_timesteps, pi_wais_salinity, label="pi_wais_lig")
ax2.set_ylabel("sum salinity")

ax2.legend()
ax2.set_title("Global")


#%% Calculate river discharge value from Antarctica
lig_wais_timesteps = lig_wais["time"]
lig_wais_salinity = []
lig_wais_rivdis = []
for i in list(range(len(lig_wais["time"]))):
    lig_wais_riv = lig_wais['F_rivdis'][i, :16, :].sum()
    lig_wais_rivdis.append(float(lig_wais_riv))
    
pi_wais_rivdis = []
pi_wais_salinity=[]
for i in list(range(len(pi_wais["time"]))):
    pi_wais_riv = pi_wais['F_rivdis'][i, :16, :].sum()
    pi_wais_rivdis.append(float(pi_wais_riv))

pi_timesteps = picontrol["time"]
pi_control_rivdis = []
pi_control_salinity=[]
for i in list(range(len(picontrol["time"]))):
    pi_control_riv = picontrol['F_rivdis'][i, :16, :].sum()
    pi_control_rivdis.append(float(pi_control_riv))

lig_wais_salinity = []
for i in list(range(len(lig_wais["time"]))):
    lig_wais_sal = lig_wais['O_sal'][i, :, :16, :].sum()
    lig_wais_salinity.append(float(lig_wais_sal))
    
pi_wais_salinity=[]
for i in list(range(len(pi_wais["time"]))):
    pi_wais_sal = pi_wais['O_sal'][i, :, :16, :].sum()
    pi_wais_salinity.append(float(pi_wais_sal))

pi_timesteps = picontrol["time"]
pi_control_salinity=[]
for i in list(range(len(picontrol["time"]))):
    pi_control_sal = picontrol['O_sal'][i, :, :16, :].sum()
    pi_control_salinity.append(float(pi_control_sal))

fig, ax = plt.subplots()

ax.plot(pi_timesteps, pi_control_rivdis, label="picontrol")
ax.plot(lig_wais_timesteps, lig_wais_rivdis, label="lig_wais_lig")
ax.plot(lig_wais_timesteps, pi_wais_rivdis, label="pi_wais_lig")

ax.set_ylabel("sum rivdis")
ax.legend()

#ax2 = ax.twinx()
# ax2.plot(pi_timesteps, pi_control_salinity, linestyle = '--')
# ax2.plot(lig_wais_timesteps, lig_wais_salinity, linestyle = '--')
# ax2.plot(lig_wais_timesteps, pi_wais_salinity, linestyle = '--')

ax.set_title("Antarctica")

#%% zonal sum for last timestep

lats = lig_wais["lat"]
lig_wais_rivdis = []
for i in list(range(len(lig_wais["lat"]))):
    lig_wais_ant = lig_wais['F_rivdis'][-1, i, :].sum()
    lig_wais_rivdis.append(float(lig_wais_ant))
    
pi_wais_rivdis = []
for i in list(range(len(pi_wais["lat"]))):
    pi_wais_ant = pi_wais['F_rivdis'][-1, i, :].sum()
    pi_wais_rivdis.append(float(pi_wais_ant))

pi_rivdis = []
for i in list(range(len(picontrol["lat"]))):
    pi_ant = picontrol['F_rivdis'][-1, i, :].sum()
    pi_rivdis.append(float(pi_ant))

lats = lig_wais["lat"]
lig_wais_salinity = []
for i in list(range(len(lig_wais["lat"]))):
    lig_wais_sal = lig_wais['O_sal'][-1, :, i, :].sum(dim=["depth","lon"])
    lig_wais_salinity.append(float(lig_wais_sal))
    
pi_wais_salinity = []
for i in list(range(len(pi_wais["lat"]))):
    pi_wais_sal = pi_wais['O_sal'][-1, :, i, :].sum(dim=["depth","lon"])
    pi_wais_salinity.append(float(pi_wais_sal))

pi_control_salinity = []
for i in list(range(len(picontrol["lat"]))):
    pi_sal = picontrol['O_sal'][-1, :, i, :].sum(dim=["depth","lon"])
    pi_control_salinity.append(float(pi_sal))

fig, ax = plt.subplots()
#ax2 = ax.twinx()

ax.plot(lats, pi_rivdis, label="picontrol")
ax.plot(lats, lig_wais_rivdis, label="lig_wais_lig")
ax.plot(lats, pi_wais_rivdis, label="pi_wais_lig")
ax.set_ylabel("sum rivdis")

ax.set_title("Zonal total river discharge")
ax.legend()

#%% water balance

