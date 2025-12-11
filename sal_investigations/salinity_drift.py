#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 16:44:37 2024

@author: mackenfr
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd
from matplotlib.cm import ScalarMappable


#%%
path_to_ref = "/Users/mackenfr/uvic_outputs/PRD-PI/280/tavg.nc" 

path_to_lig = "/Users/mackenfr/uvic_outputs/AMN-PI/280_nowind/tavg_cat_sal.nc" 
path_to_control = "/Users/mackenfr/uvic_outputs/PRD-PI/nowind/280_nowind/tavg_cat_sal.nc" 

pi_wais_lig = xr.open_dataset(path_to_lig, decode_times=False)
picontrol = xr.open_dataset(path_to_control, decode_times=False)
ref = xr.open_dataset(path_to_ref, decode_times=False).rename_dims({"latitude": "lat", "longitude": "lon"}).rename_vars({"latitude": "lat", "longitude": "lon"})

#%%
ocean_levels = list(range(19))

cmap = plt.get_cmap('viridis')  # Choose any other colormap
norm = mcolors.Normalize(vmin=0, vmax=len(ocean_levels)-1)


def plot_sal_drift_basins(data, min_lon, max_lon, title=None):   
    
    if min_lon > max_lon:
         lon_range = list(range(0,max_lon+1))+list(range(min_lon, 100))
    else:
         lon_range = list(range(min_lon, max_lon+1))
    maximum = ref["G_areaT"][:,lon_range].max()
    weights = (ref["G_areaT"][:,lon_range]/maximum)
    fig, ax = plt.subplots()
    for i, level in enumerate(ocean_levels): #for ax, level in zip(axs.flat, upper_levels):
        weighted_sal = pi_wais_lig["O_sal"][50:,level, :, lon_range].weighted(weights)
        mean_sal = weighted_sal.mean(dim=["lon", "lat"])
        salinity_drift=(mean_sal)/mean_sal.mean()
        color = cmap(norm(i))  # Get color from colormap
        ax.plot(pi_wais_lig["time"][50:], salinity_drift, color=color)
    ax.set_title(title)
    #ax.set_ylim(-0.0002,0.0006)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # You need to set an empty array to use with ScalarMappable
    
    # Add the colorbar to the figure
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('ocean_levels') 

def plot_sal_drift_lat(data, min_lat, max_lat, title=None):   
    maximum = ref["G_areaT"][min_lat:max_lat,:].max()
    weights = (ref["G_areaT"][min_lat:max_lat,:]/maximum)
    fig, ax = plt.subplots()
    for i, level in enumerate(ocean_levels): #for ax, level in zip(axs.flat, upper_levels):
        weighted = pi_wais_lig["O_sal"][50:,level, min_lat:max_lat, :].weighted(weights)
        mean_sal = weighted.mean(dim=["lon", "lat"])
        salinity_drift=(mean_sal - mean_sal[0])/mean_sal.mean()
        color = cmap(norm(i))  # Get color from colormap
        ax.plot(pi_wais_lig["time"][50:], salinity_drift, color=color)
    ax.set_ylim(-0.0002,0.0012)
    ax.set_title(title)
    
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # You need to set an empty array to use with ScalarMappable
    
    # Add the colorbar to the figure
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('ocean_levels') 


#%%


plot_sal_drift_basins(pi_wais_lig, 0,99, "pi_wais_lig global salinity drift")    
plot_sal_drift_basins(pi_wais_lig, 80,5, "pi_wais_lig atlantic salinity drift")    
plot_sal_drift_basins(pi_wais_lig, 35, 80, "pi_wais_lig pacfific salinity drift")    
plot_sal_drift_basins(pi_wais_lig, 5,35, "pi_wais_lig indian salinity drift")    

#%%

plot_sal_drift_lat(pi_wais_lig, 0,24, "pi_wais_lig 90S-45S salinity drift")    
plot_sal_drift_lat(pi_wais_lig, 24,50, "pi_wais_lig 45S-0 salinity drift")    
plot_sal_drift_lat(pi_wais_lig, 50,74, "pi_wais_lig 0-45N salinity drift")    
plot_sal_drift_lat(pi_wais_lig,74, 99,"pi_wais_lig 45N-90N salinity drift")    

#%%
maximum = ref["G_areaT"][:,:].max()
weights = (ref["G_areaT"][:,:]/maximum)
fig, ax = plt.subplots()
weighted = pi_wais_lig["O_sal"][50:,:, :, :].weighted(weights)
mean_sal = weighted.mean(dim=["lon", "lat", "depth"])
salinity_drift=mean_sal
ax.plot(pi_wais_lig["time"][50:], salinity_drift)
ax.set_title("global mean sal")


# %%

sal_diffs = []
for i in ocean_levels:
    weighted_5k = pi_wais_lig["O_sal"][50,i, :, :].weighted(weights)
    mean_sal_5k = weighted_5k.mean(dim=["lon", "lat"])
    weighted_10k = pi_wais_lig["O_sal"][-1,i, :, :].weighted(weights)
    mean_sal_10k = weighted_10k.mean(dim=["lon", "lat"])
    sal_diff = mean_sal_10k.values - mean_sal_5k.values
    sal_diffs.append(sal_diff)

diff_df = pd.DataFrame({"dS": sal_diffs,
                       "level": ocean_levels,
                       "depth": pi_wais_lig["depth"].values})
#%%
fig, ax = plt.subplots()
sns.scatterplot(x = "depth", y="dS", data=diff_df)
# %%
path_to_amn = "/Users/mackenfr/uvic_outputs/AMN-PI/280/tsi_all.nc" 
path_to_control = "/Users/mackenfr/uvic_outputs/PRD-PI/280/tsi_all.nc" 
amn_tsi = xr.open_dataset(path_to_amn, decode_times=False)
prd_tsi = xr.open_dataset(path_to_control, decode_times=False)

# %%
fig, ax = plt.subplots()
sns.lineplot(x = prd_tsi.time[5000:], y=prd_tsi.O_salsur[5000:])
ax.set_title("prd salsur")

fig, ax = plt.subplots()
sns.lineplot(x = prd_tsi.time[5000:], y=prd_tsi.O_sal[5000:])
ax.set_title("prd sal")

# %%
fig, ax = plt.subplots()
sns.lineplot(x = amn_tsi.time[5000:], y=amn_tsi.O_salsur[5000:])
ax.set_title("amn salsur")

fig, ax = plt.subplots()
sns.lineplot(x = amn_tsi.time[5000:], y=amn_tsi.O_sal[5000:])
ax.set_title("amn sal")

# %%
