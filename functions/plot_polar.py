#import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import colormaps
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np

os.chdir("/Users/mackenfr/scripts")

import pandas as pd
from import_data import import_data
#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
tavgs, ref_df = import_data(path_to_data,  config_names=["PRD-PI", "AMN-PI", "TRS-PI", "EQB-PI"])
picontrol = tavgs[0]
#%%
def get_2d_slice(da, lev=0):
    if "time" in da.dims:
        da = da.isel(time=0)
    if "depth" in da.dims:
        da = da.isel(depth=lev)
    if "depth_W" in da.dims:
        da = da.isel(depth_W=lev)
    return da

def plot_polar_grid(var, cmap, vmin, vmax, lev = 0, proj = ccrs.SouthPolarStereo(), title="", anom=None): 
    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = -45
        figsize=(15,12)
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 45
        north_bound = 90
        figsize=(15,12)
    elif proj == ccrs.PlateCarree():
        # south_bound = -90
        # north_bound = 90
        # figsize=(15,8)
        south_bound = -90
        north_bound = 90
        figsize=(15,15)
        
    if anom=="picontrol":
        suptitle=title + " anom from PRD 280"
        anom_idx = [0 for i in range(len(ref_df.config))]
    elif anom=="config":
        suptitle=title + " anom from PRD"
        anom_idx = [(i % 5) for i in range(len(ref_df.config))]
    elif anom=="co2":
        suptitle=title + " anom from 280 ppm"
        anom_idx = [(i // 5) * 5 for i in range(len(ref_df.config))]
    else:
        suptitle=title

    if anom:
        anomalies = [get_2d_slice(tavgs[i][var], lev) - get_2d_slice(tavgs[anom_idx[i]][var], lev) for i in range(len(tavgs))]
    else:
        anomalies = [get_2d_slice(tavgs[i][var], lev) for i in range(len(tavgs))]

  
    fig, axs = plt.subplots(4, 5, figsize=figsize, subplot_kw={'projection': proj})
    fig.suptitle(suptitle, fontsize="xx-large", y=0.93)

    lons = tavgs[0]["longitude"]
    lats = tavgs[0]["latitude"]

    for i, ax in enumerate(axs.flat):
        data = anomalies[i]

            
        ax.set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(),
                        vmin=vmin, vmax=vmax, 
                        cmap=cmap, shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
    
    for n, col in enumerate(ref_df.co2.unique()):
        axs[0, n].set_title(col)

    for n, row in enumerate(ref_df.config.unique()):
        axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.8 - (n * 0.19), row, va='center', ha='right', fontsize=12, rotation=90)

    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label(f"{tavgs[0][var].attrs.get('long_name', var)} ({tavgs[0][var].attrs.get('units', '')})", fontsize=14)
    #cb.ax.tick_params(labelsize=14)

#     #return fig
# #%%
# # for i in list(range(19)):
# #     plot_polar_grid("O_velZ", "RdYlBu", -0.00001, 0.00001, lev=i, title="", anom="co2")

# # %%
# plot_polar_grid("O_velZ", "RdYlBu",  -0.0000005, 0.0000005, title="velz", anom="config")
# # %%
# plot_polar_grid("O_ventdep", "RdYlBu_r",  -100, 100, title="ventdepth", anom="config")
# plot_polar_grid("O_ventdep", "RdYlBu_r",  -100, 100, title="ventdepth", anom="co2")
# # %%
# plot_polar_grid("F_heat", "jet",  -50, 50, title="heat flux", anom="config")
# plot_polar_grid("F_heat", "jet",  -50, 50, title="heat flux", anom="co2")


# # %%
# plot_polar_grid("O_temp", "coolwarm",  -1, 1, title="O_temp", anom="config")
# plot_polar_grid("O_temp", "coolwarm",  -4, 4, title="O_temp", anom="co2")
# # %%
# plot_polar_grid("O_sal", "PiYG", -1, 1, title="salinity", anom="config")
# plot_polar_grid("O_sal", "PiYG",  -1, 1, title="salinity", anom="co2")
# # %%
# plot_polar_grid("F_precip", "PuOr", -0.00001, 0.00001, title="precip", anom="config")
# plot_polar_grid("F_precip", "PuOr",  -0.00001, 0.00001, title="precip", anom="co2")
# # %%
# #plot_polar_grid("F_precip", "PuOr",  0, 0.0001, proj=ccrs.PlateCarree(), title="precip")
# plot_polar_grid("O_temp", "coolwarm",  -4, 4, proj=ccrs.PlateCarree(), title="otemp",anom="co2")
# plot_polar_grid("O_temp", "coolwarm",  -1, 1, proj=ccrs.PlateCarree(), title="otemp",anom="config")

# %%

# %%
