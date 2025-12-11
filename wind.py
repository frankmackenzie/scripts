#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 15:49:09 2025

@author: mackenfr
"""

"""
Created on Mon Jan 27 14:09:46 2025

@author: mackenfr adapted from keller
"""



import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import matplotlib.cm as cm
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import pandas as pd

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

from scipy.interpolate import griddata
import os
#from matplotlib import animation
os.chdir("/Users/mackenfr/scripts")

import pandas as pd
from import_data import import_data
#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
tavgs, ref_df = import_data(path_to_data, ["PRD-PI", "AMN-PI", "TRS-PI","EQB-PI"])
picontrol = tavgs[0]
from batlow import *  
#%%
def plot_winds_global(u,v, title,projection, vmin, vmax, c=None, cmap=None):
    step=3
    fig = plt.figure()

    ax = fig.add_subplot(1, 1, 1, projection=projection)
    #ax = plt.axes(projection=ccrs.SouthPolarStereo())
    if projection == ccrs.PlateCarree():
        ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    else:
        ax.set_extent([-180, 180, -90, 60], ccrs.PlateCarree())
    #windspeed colors background
    #import pdb; pdb.set_trace()
   # if not c is None:
   #     color = plt.pcolormesh(xx, yy, c,cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    NORMZ = colors.Normalize(vmin=0, vmax=1)
   
    qt=ax.quiver(xx[::step,::step], yy[::step,::step], u[::step,::step], v[::step,::step],c[::step,::step], transform=ccrs.PlateCarree(),
                 cmap=cmap,angles="xy", scale=1, scale_units = "dots")
    plt.quiverkey(qt,X=0.75,Y=-0.1,U=12, label='12 m/s', labelpos='E')
    plt.title(title)
    ax.coastlines() 
    #ax.gridlines(linestyle='--')
    #ax.set_boundary(circle, transform=ax.transAxes)

    # Color bar
    if not c is None:
        label="winspeed (m/s)"
        cbar = fig.colorbar(qt, ax=ax, fraction=0.027, pad=0.045, 
                            orientation="horizontal")
        cbar.set_label(label, rotation=0, labelpad=5, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

#%% load data  

# path_to_data = "/Users/mackenfr/uvic_outputs/" 

# configs = ["picontrol","lig_wais_pi", "pi_wais_lig",  "lig_wais_lig"]
# ice = ["PRD-PI", "PRD-LIG", "AMN-PI", "AMN-LIG"]
# co2s = ["280", "300", "350", "400", "450"]#
# co2ccn = ["280", "300", "350", "400", "450"]

# si_maxs = []
# si_mins = []
# tavgs = []
# names = []
# co2_all=[]
# for config in ice:
#     for co2 in co2s:
#         filename_tavg = path_to_data + "/" + config + "/" + co2 + "_wind/tavg.nc"
#         data = xr.open_mfdataset(filename_tavg, decode_times=False, engine = "scipy")
#         names.append(config)
#         tavgs.append(data)
#         co2_all.append(co2)
        
# ref_df = pd.DataFrame({"config": names,
#                        "co2": co2_all})

#%% reference data
path_to_refs = "/Users/mackenfr/UVic_ESCM/data/"
ref_wind = xr.open_mfdataset(path_to_refs + "A_wind.nc")
ref_windsur = xr.open_mfdataset(path_to_refs + "A_windsur.nc")
#%%
lon = tavgs[0]["longitude"]
lat = tavgs[0]["latitude"]
lonV = tavgs[0]["longitude_V"]
latV = tavgs[0]["latitude_V"]
x, y = np.meshgrid(lon,lat)
xx, yy = np.meshgrid(lonV,latV)

#%% input winds, speed

ang = ref_windsur["A_windang"]
mag = ref_windsur["A_windspd"]

windX = mag * np.cos(ang)  
windY = mag * np.sin(ang)

input_windXang = windX.mean(dim="time")
input_windYang = windY.mean(dim="time")

input_windspd = mag.mean(dim="time")

input_windX = input_windXang
input_windY = input_windYang

#%% winds, q & t

windqX = ref_wind["A_windqX"]
windqY = ref_wind["A_windqY"]
input_windqX = windqX.mean(dim="time")
input_windqY = windqY.mean(dim="time")

windtX = ref_wind["A_windtX"]
windtY = ref_wind["A_windtY"]
input_windtX = windtX.mean(dim="time")
input_windtY = windtY.mean(dim="time")

output_windqX = tavgs[0]["A_windqX"][0,:,:]
output_windqY = tavgs[0]["A_windqY"][0,:,:]

output_windtX = tavgs[0]["A_windtX"][0,:,:]
output_windtY = tavgs[0]["A_windtY"][0,:,:]

#%%
output_awindX = tavgs[10]["A_awindX"][0,:,:]
output_awindY = tavgs[10]["A_awindY"][0,:,:]
output_windspd = tavgs[10]["A_windspd"][0,:,:]

total_windX = input_windX.values + output_awindX.values
total_windY = input_windY.values + output_awindY.values
#%%
plot_winds_global(input_windX.values, input_windY.values, "base winds from input", projection=ccrs.PlateCarree(), vmin=0, vmax=20,c=input_windspd, cmap=batlow_map)


#magnitude of awind very small compared to base windfields - not helpful to plot, even for non-control
#%%
#windspeed anomaly
anom = output_windspd.values - input_windspd.values
plot_winds_global(total_windX, total_windY, "windspeed anomaly prd_amn280 - base", projection = ccrs.PlateCarree(),vmin=-2, vmax=2, c=anom, cmap="PiYG")


#%%
step=2
titles = [i + j for i in ref_df["ice"].unique() for j in ref_df["co2"].unique()]
savepath = "/Users/mackenfr/plots_organised/winds/awind_windspd_anoms/"
for title, tavg in zip(titles, tavgs):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    u=tavg["A_awindX"][0,:,:].values - tavgs[0]["A_awindX"][0,:,:].values
    v=tavg["A_awindY"][0,:,:].values - tavgs[0]["A_awindY"][0,:,:].values 
    c=tavg["A_windspd"][0,:,:].values - tavgs[0]["A_windspd"][0,:,:].values
    ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    ax.coastlines() 
    
    #ax.patch.set_alpha(0)   
    ax.set_facecolor('white')
    #ax.gridlines(linestyle='--')
    #color = ax.pcolormesh(xx, yy, c,cmap="PiYG", transform=ccrs.PlateCarree(), vmin=-1, vmax=1)
    NORMZ = colors.Normalize(vmin=-1, vmax=1)
    qt=ax.quiver(xx[::step,::step], yy[::step,::step], u[::step,::step], v[::step,::step], c[::step,::step], cmap="PiYG", norm=NORMZ,transform=ccrs.PlateCarree(),
              angles="xy", scale=0.1, scale_units = "dots")
    plt.title(title)
    #plt.quiverkey(qt,X=0.75,Y=-0.1,U=10, label='10 m/s', labelpos='E')
    cb = fig.colorbar(qt, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label("wind speed m/s", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    #plt.savefig(savepath + title+".png")
  #%%
step = 1
anom = "config"
#for i, data in enumerate(tavgs[5]):
i=9
if anom=="picontrol":
    anom_idx = [0 for i in range(len(ref_df.config))]
elif anom=="config":
    suptitle=title + " anom from PRD"
    anom_idx = [(i % 5) for i in range(len(ref_df.config))]
elif anom=="co2":
    suptitle=title + " anom from 280 ppm"
    anom_idx = [(i // 5) * 5 for i in range(len(ref_df.config))]
else:
    suptitle=title
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
u=tavgs[i]["A_awindX"][0,:,:].values - tavgs[14]["A_awindX"][0,:,:].values
v=tavgs[i]["A_awindY"][0,:,:].values - tavgs[14]["A_awindY"][0,:,:].values 
c=tavgs[i]["A_windspd"][0,:,:].values - tavgs[14]["A_windspd"][0,:,:].values
ax.set_extent([-180, -20, -90, -50], ccrs.PlateCarree())
ax.coastlines() 

#ax.patch.set_alpha(0)   
ax.set_facecolor('white')
#ax.gridlines(linestyle='--')
color = ax.pcolormesh(xx, yy, c,cmap="PiYG", transform=ccrs.PlateCarree(), vmin=-10, vmax=10)

NORMZ = colors.Normalize(vmin=-1, vmax=1)
qt=ax.quiver(xx[::step,::step], yy[::step,::step], u[::step,::step], v[::step,::step], #c[::step,::step], cmap="PiYG", norm=NORMZ,
        transform=ccrs.PlateCarree(),angles="xy", scale=0.1, scale_units = "dots")
plt.title(ref_df["ice"][i] +  ref_df["co2"][i] + "- TRS 450")
plt.quiverkey(qt,X=0.75,Y=-0.1,U=10, label='10 m/s', labelpos='E')
cb = fig.colorbar(color, orientation='horizontal', fraction=0.05, pad=0.04)
cb.set_label("wind speed m/s", fontsize=14)
cb.ax.tick_params(labelsize=14)      

#%%
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

u=tavgs[19]["A_awindX"][0,:,:].values #- tavgs[0]["A_awindX"][0,:,:].values 

ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
ax.coastlines() 

ax.set_facecolor('white')
c = ax.pcolormesh(xx, yy, u,cmap="PiYG", transform=ccrs.PlateCarree(), vmin=-1, vmax=1)

plt.title("amnlig450 awindX")
#plt.quiverkey(qt,X=0.75,Y=-0.1,U=10, label='10 m/s', labelpos='E')
cb = fig.colorbar(c, orientation='horizontal', fraction=0.05, pad=0.04)
cb.set_label("wind speed m/s", fontsize=14)
cb.ax.tick_params(labelsize=14)  
#%%

df_master = pd.DataFrame(columns=["config",
                          "latitudes",
                          "280",
                          "300",
                          "350",
                          "400",
                          "450"])

for name in ref_df["ice"].unique():
    plt.figure(figsize=(10, 4))
    indexes = list(ref_df.query("ice == @name").index)
    co2s_config = ref_df["co2"][indexes[0]:(indexes[-1]+1)]
    tavgs_config = tavgs[indexes[0]:(indexes[-1]+1)]
    latitudes = tavgs[0]["latitude_V"][10:30]
    df_config = pd.DataFrame({"latitudes": latitudes})
    zonal_means=[]
    
    for data, co2 in zip(tavgs_config, co2s_config):
        zonal_mean=data["A_awindX"][0,10:30,:].mean(dim="longitude_V")
        #zonal_means.append(zonal_mean)
        df_config[co2]=zonal_mean
    df_config["config"] = name
    df_master = pd.concat([df_config, df_master])
    
    melted = (df_config.drop(columns="config")
        .melt(id_vars = "latitudes", var_name="co2", value_name="wind_spd")
        )
   
    #import pdb; pdb.set_trace()
    sns.scatterplot(x="latitudes", y="wind_spd", data=melted, hue="co2", palette="viridis_r")
    sns.lineplot(x="latitudes", y="wind_spd", data=melted, hue="co2", legend=False, palette="viridis_r")
    plt.ylim(-1, 2)
    plt.title(name)
    plt.ylabel("Windspeed anomaly m/s")
    plt.xlabel("Latitude")
    
    


#%% latitude of max wind anomaly

max_spd = (df_master[['280', '300', '350', '400', '450', "config"]]
           .groupby("config")
           .max()
           .reset_index()
           .melt(id_vars=["config"], var_name="co2", value_name="max_spd"))

lat_idx = (df_master[['280', '300', '350', '400', '450', "config"]]
           .groupby("config")
           .idxmax()
           .reset_index()
           .melt(id_vars=["config"], var_name="co2", value_name="lat_idx"))

max_spd["latitude"]=latitudes.values[lat_idx["lat_idx"]]
max_spd["co2"]=max_spd["co2"].astype("int")
#%%
plt.figure(figsize=(10, 4))
ax = sns.scatterplot(x="latitude", y="max_spd", data=max_spd, hue="co2", style="config",palette="viridis_r")
ax.collections[0].set_sizes([100])
plt.ylabel("windspeed anomaly m/s")
plt.xlabel("Latitude")
plt.title("Zonal mean maximum windspeed in Westerlies")
   #%%
plt.figure(figsize=(10, 4))
ax = sns.scatterplot(x="latitude", y="co2", data=max_spd, hue="max_spd", style="config",palette="plasma_r", legend="auto")

ax.collections[0].set_sizes([100])
plt.ylabel("CO2 (ppm)")
plt.xlabel("Latitude")
plt.title("Zonal mean maximum windspeed in Westerlies")

norm = plt.Normalize(max_spd['max_spd'].min(), max_spd['max_spd'].max())
sm = plt.cm.ScalarMappable(cmap='plasma_r', norm=norm)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Max. wind speed anomaly') 

handles, labels = ax.get_legend_handles_labels()
new_handles, new_labels = [], []
for h, l in zip(handles, labels):
    if l in max_spd['config'].unique():  # Keep only style legends
        new_handles.append(h)
        new_labels.append(l)

# Add the updated legend (only style)
ax.legend(new_handles, new_labels, title="config")
#%%
plt.figure(figsize=(10, 4))
ax = sns.scatterplot(x="latitude", y="co2", data=max_spd, hue="config", style="config", legend="auto")

ax.collections[0].set_sizes([100])
plt.ylabel("CO2 (ppm)")
plt.xlabel("Latitude")
plt.title("Zonal mean maximum windspeed in Westerlies")
#%%
plt.figure(figsize=(10, 4))
ax = sns.scatterplot(x="latitude", y="max_spd", data=max_spd, hue="config", legend="auto")

ax.collections[0].set_sizes([100])
plt.ylabel("Max speed ")
plt.xlabel("Latitude")
plt.title("Zonal mean maximum windspeed in Westerlies")

#%%
cmap = sns.cubehelix_palette(as_cmap=True)
for config in ref_df["ice"].unique():
    subset = max_spd.query("config== @config")
    plt.figure(figsize=(10, 4))
    plt.ylim(-70, -50)
    plt.title(config)
    sns.lineplot(x="co2", y="latitude", data=subset, color="lightgray")
    ax = sns.scatterplot(x="co2", y="latitude", hue="max_spd", data=subset,legend=False, palette=cmap)
    norm = plt.Normalize(max_spd['max_spd'].min(), max_spd['max_spd'].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Max. wind speed anomaly') 

#%% gradient between max and equator

    