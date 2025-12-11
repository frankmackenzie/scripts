#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 16:41:45 2024

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
lig_wais_file = path_to_data + "lig_wais_lig/280/tavg_cat_water_merged.nc"
lig_wais = xr.open_dataset(lig_wais_file, decode_times=False)

pi_wais_file = path_to_data + "pi_wais_lig/280/tavg_cat_water_merged.nc"
pi_wais = xr.open_dataset(pi_wais_file, decode_times=False)


pi_control_file = path_to_data + "picontrol/280/tavg_cat_water_merged.nc"
pi_control = xr.open_dataset(pi_control_file, decode_times=False)

lon = lig_wais["lon"]
lat = lig_wais["lat"]

#%%

lig_rivdis_time = lig_wais["F_rivdis"].sum(dim = ["lat","lon"])
lig_sal_time = lig_wais["O_sal"].sum(dim = ["depth","lat","lon"])

pi_rivdis_time = pi_wais["F_rivdis"].sum(dim = ["lat","lon"])
pi_sal_time = pi_wais["O_sal"].sum(dim = ["depth","lat","lon"])

pi_control_rivdis_time = pi_control["F_rivdis"].sum(dim = ["lat","lon"])
pi_control_sal_time = pi_control["O_sal"].sum(dim = ["depth","lat","lon"])

fig, ax = plt.subplots()
ax.scatter(lig_rivdis_time, lig_sal_time)
ax.set_title("lig_wais")
ax.set_xlabel("global total rivdis")
ax.set_ylabel("global total sal")


fig, ax = plt.subplots()
ax.scatter(pi_rivdis_time, pi_sal_time)
ax.set_title("pi_wais")
ax.set_xlabel("global total rivdis")
ax.set_ylabel("global total sal")

fig, ax = plt.subplots()
ax.scatter(pi_control_rivdis_time, pi_control_sal_time)
ax.set_title("picontrol")
ax.set_xlabel("global total rivdis")
ax.set_ylabel("global total sal")
