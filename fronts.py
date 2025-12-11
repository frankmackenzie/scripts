import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import colormaps
from matplotlib.lines import Line2D

import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
import gsw
os.chdir("/Users/mackenfr/scripts")

import pandas as pd
path_to_data = "/Users/mackenfr/uvic_outputs/" 
config_names = ["PRD-PI", "AMN-PI", "TRS-PI", "EQB-PI"] # "PRD-LIG","AMN-LIG"
co2 = "280"
tavgs = []
sea_ice = []
names = []
ices = []
orbs = []
for config in config_names:
    if config == "PRD-PI" or config == "AMN-PI":
            filename_si = path_to_data + "/" + config + "/" + co2 + "/sea_ice_annual.nc"
    if config == "EQB-PI" or config == "TRS-PI":
            filename_si = path_to_data + "/" + config + "/" + co2 + "/sea_ice_mean.nc"
    si_data = xr.open_dataset(filename_si, decode_times=False)
    filename_tavg = path_to_data + "/" + config + "/"+co2+"/tavg.nc"
    data = xr.open_dataset(filename_tavg, decode_times=False, engine = "netcdf4")
    ice, orb = config.split("-")
    sea_ice.append(si_data)
    names.append(config)
    tavgs.append(data)
    ices.append(ice)
    orbs.append(orb)

ref_df = pd.DataFrame({"config": names,
                    "ice": ices,
                    "orb": orbs
                    })
                           
picontrol = tavgs[0]
#%%
#%%
sal = picontrol["O_sal"]
temp = picontrol["O_temp"]

test = picontrol["O_temp"].where(((sal[0, :4, :, :]>34)&((temp[0, 4, :, :]>5))))
mask = ~test.isnull()
new = xr.where(mask[0,0,:,:], 0, 1)
new_masked = new.where(picontrol["G_mskt"]==1, np.nan)



#data = temp[0,0,1]
#%%
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.SouthPolarStereo()})

lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]

ax.set_extent([-180, 180,-90, -40], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
cm = ax.pcolormesh(lons, lats,picontrol["O_temp"][0,0,:,:], transform=ccrs.PlateCarree(),
                vmin=0, vmax=20, 
                cmap="Spectral_r", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
cols = ["green", "cyan", "blue", "red"]
legend_lines = []
for i, data in enumerate(tavgs):
    sal = data["O_sal"]
    temp = data["O_temp"]
    test = data["O_temp"].where(((sal[0, :4, :, :]>34)&((temp[0, 4, :, :]>5))))
    mask = ~test.isnull()
    new = xr.where(mask[0,0,:,:], 0, 1)
    new_masked = new.where(picontrol["G_mskt"]==1, np.nan)
    c1 = ax.contour(lons, lats,new_masked, colors=cols[i], levels = [0,1], transform=ccrs.PlateCarree())
    legend_lines.append(Line2D([0], [0], color=cols[i], lw=2))
cb = fig.colorbar(cm, ax=ax, orientation='horizontal', fraction=0.05, pad=0.04)

ax.legend(legend_lines, list(ref_df.ice), loc='lower left', title='Experiments')
plt.title("SAF")
# %%
fig, axs = plt.subplots(2,2,subplot_kw={'projection': ccrs.SouthPolarStereo()})
legend_lines = []
lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]
for i, ax in enumerate(axs.flat):
    ax.set_extent([-180, 180,-90, -40], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    cm = ax.pcolormesh(lons, lats,picontrol["O_temp"][0,0,:,:], transform=ccrs.PlateCarree(),
                vmin=0, vmax=20, 
                cmap="Spectral_r", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))

    ax.set_title(ref_df.ice[i])

    c1 = ax.contour(lons, lats,tavgs[i]["O_temp"][0, 2, :, :], colors="black", levels = [2], transform=ccrs.PlateCarree())
    c1 = ax.contour(lons, lats,tavgs[i]["O_temp"][0, 6, :, :], colors="black", levels = [2], transform=ccrs.PlateCarree())
cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)

fig.suptitle("PF")

# %%

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.SouthPolarStereo()})

lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]

ax.set_extent([-180, 180,-90, -40], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
cm = ax.pcolormesh(lons, lats,picontrol["G_mskt"], transform=ccrs.PlateCarree(),
                vmin=0, vmax=1, 
                cmap="Greys_r", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
cols = ["green", "cyan", "blue", "red"]
legend_lines = []
for i, data in enumerate(tavgs):
    sal = data["O_sal"]
    temp = data["O_temp"]
    test = data["O_temp"].where(((sal[0, :4, :, :]>34)&((temp[0, 4, :, :]>5))))
    mask = ~test.isnull()
    new = xr.where(mask[0,0,:,:], 0, 1)
    new_masked = new.where(picontrol["G_mskt"]==1, np.nan)
    c1 = ax.contour(lons, lats,new_masked, colors=cols[i], levels = [0,1], transform=ccrs.PlateCarree())
    legend_lines.append(Line2D([0], [0], color=cols[i], lw=2))

ax.legend(legend_lines, list(ref_df.ice), loc='lower left')
plt.title("SAF")
# %%
fig, axs = plt.subplots(2,2,subplot_kw={'projection': ccrs.SouthPolarStereo()})
legend_lines = []
lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]
for i, ax in enumerate(axs.flat):
    ax.set_extent([-180, 180,-90, -40], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    cm = ax.pcolormesh(lons, lats,picontrol["G_mskt"], transform=ccrs.PlateCarree(),
                vmin=0, vmax=1, 
                cmap="Greys_r", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))

    ax.set_title(ref_df.ice[i])

    #c1 = ax.contour(lons, lats,tavgs[i]["O_temp"][0, 2, :, :], colors=cols[i], levels = [2], transform=ccrs.PlateCarree())
    c1 = ax.contour(lons, lats,tavgs[i]["O_temp"][0, 6, :, :], colors=cols[i], levels = [2], transform=ccrs.PlateCarree())

fig.suptitle("PF")
# %%
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.SouthPolarStereo()})

lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]

ax.set_extent([-180, 180,-90, -40], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
cm = ax.pcolormesh(lons, lats,picontrol["G_mskt"], transform=ccrs.PlateCarree(),
                vmin=0, vmax=1, 
                cmap="Greys_r", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
cols = ["green", "cyan", "blue", "red"]
legend_lines = []
for i, data in enumerate(tavgs):
    c1 = ax.contour(lons, lats,tavgs[i]["O_temp"][0, 6, :, :], colors=cols[i], levels = [2], transform=ccrs.PlateCarree())
    legend_lines.append(Line2D([0], [0], color=cols[i], lw=2))

ax.legend(legend_lines, list(ref_df.ice), loc='lower left')
plt.title("PF")

# %%
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.SouthPolarStereo()})

lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]

ax.set_extent([-180, 180,-90, -40], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
cm = ax.pcolormesh(lons, lats,picontrol["G_mskt"], transform=ccrs.PlateCarree(),
                vmin=0, vmax=1, 
                cmap="Greys_r", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
cols = ["green", "cyan", "blue", "red"]
legend_lines = []
for i, data in enumerate(tavgs):
    c1 = ax.contour(lons, lats,tavgs[i]["O_velX"][0, 0, :, :], colors=cols[i], levels = [0.05], transform=ccrs.PlateCarree())
    legend_lines.append(Line2D([0], [0], color=cols[i], lw=2))

ax.legend(legend_lines, list(ref_df.ice), loc='lower left')
plt.title("ACC")
# %%
