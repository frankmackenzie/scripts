#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:22:00 2024

@author: mackenfr

Look at time slices of high res data containing DCE to generate Figure 13 in my proposal
"""
from matplotlib import animation
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cftime
import numpy as np
import pandas as pd

#%% define functions


    
def plot_polar_stereo(lon, lat, value,cmap,vmin,vmax,title):
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    ax.coastlines() 
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")


    ct=ax.pcolormesh(lon, lat, value,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")

    #colorbar 
    cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)
    cb.set_label(value.attrs["long_name"] + "(" + value.attrs["units"] + ")")

    
    #plt.title(title)
    plt.tight_layout()
    #fig.savefig(title +'.png', transparent=True)
    
#%% load data

file = "/Users/mackenfr/uvic_outputs/AMN-PI/280/seaice5day.nc"#tavg_pi_wais_lig_5day.nc"
data = xr.open_dataset(file, decode_times=False)
icefra = data["O_icefra"]
lon = data["lon"]
lat = data["lat"]

#%% monthly averages
#Generate new time variable that xarray can deal with 
units = 'days since 1800-01-01 00:00'
time_365 = cftime.num2date(np.arange(0, len(data.time)*5, 5), units, '365_day')
time_adj = xr.DataArray(np.arange(time_365.size), coords = [time_365], dims = 'time', name = 'data')
icefra["time"]= time_adj["time"]
#%%
#define month as 30 days (mean of 6 5day snapshots)
icefra['time'] = np.arange(len(icefra['time'])) // 6
grouped_mean = icefra.groupby('time').mean()
grouped_mean_all = icefra.groupby('time').mean(dim=["lat", "lon"])
plt.plot(grouped_mean_all)
print(grouped_mean)

#%%
for x in list(range(310,330)):
    plot_polar_stereo(lon, lat, grouped_mean[x, :,:], "Blues_r", 0, 1, "Ice fraction")
# %%
file = "/Users/mackenfr/uvic_outputs/AMN-LIG/400/sea_ice_5daymean.nc"#tavg_pi_wais_lig_5day.nc"
data = xr.open_dataset(file, decode_times=False, engine="netcdf4")
# %%
for x in list(range(0,73)):
    plot_polar_stereo(lon, lat, data["O_icefra"][x,:,:], "Blues_r", 0, 1, "Ice fraction")

# %%

# %%
