#%%
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import colormaps
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
import gsw
os.chdir("/Users/mackenfr/scripts")

import pandas as pd

from import_data import import_data
from calc_ovt import calc_ovt
from common_func import calculate_mean_density, calculate_density, calc_mean, get_2d_slice
#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
config_names = ["PRD-PI", "AMN-PI", "EQB-PI", "PRD-LIG", "AMN-LIG", "EQB-LIG"] # "PRD-LIG","AMN-LIG"
co2 = ["280"]
tavgs, ref_df = import_data(path_to_data, co2, config_names)
#%%
subset = ref_df.loc[ref_df["orb"]=="LIG"].reset_index()
control_subset = ref_df.loc[ref_df["orb"]=="PI"].reset_index()
datasets = [tavgs[i] for i in subset["index"]]
control_datasets = [tavgs[i] for i in control_subset["index"]] 
#%%
def surface_ver2(var, cmap, vmin, vmax, lev = 0, proj = ccrs.SouthPolarStereo(), nbound = -45, contour = None, levels = None,  title=""): 
    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = nbound
        figsize=(10,3)
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 35
        north_bound = 90
        figsize=(7,4)
    elif proj == ccrs.PlateCarree():
        # south_bound = -90
        # north_bound = 90
        # figsize=(15,8)
        south_bound = -90
        north_bound = 90
        figsize=(10,2)

    datasets = [get_2d_slice(tavgs[i][var], lev) for i in subset["index"]]
    control_datasets = [get_2d_slice(tavgs[i][var], lev) for i in control_subset["index"]]  

    fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': proj})
    fig.suptitle(title, fontsize=14, y=1)


    lons = tavgs[0]["longitude"]
    lats = tavgs[0]["latitude"]

    for i, ax in enumerate(axs.flat):
        data = datasets[i] - control_datasets[i]
            
        ax.set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(),
                        vmin=vmin, vmax=vmax, 
                        cmap=cmap, shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
        if contour:
            c1 = ax.contour(lons, lats,contour[i], colors=["black"], levels = levels, transform=ccrs.PlateCarree())
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
        ax.set_title(subset["ice"][i] + ", " + subset["orb"][i] + " - " + control_subset["orb"][i], fontsize=10, pad=2)
        fig.suptitle(title, fontsize=14, va="top")
    cb = fig.colorbar(axs[0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.05)
    cb.set_label(f"{tavgs[0][var].attrs.get('long_name', var)} ({tavgs[0][var].attrs.get('units', '')})", fontsize=14)

#%%
def transect_ver2(var, basin, cmap, vmin,vmax,contour = None, extent="Global",  title="", datasets = datasets, control_datasets=control_datasets):
    if extent == "SO":
        xlim=(-90, -50)
        xextent = 24
        figsize = (10,6)

    else:
        xlim=(-90, 90)
        xextent = 100
        figsize = (12,4)
 
    
    values = []
    density = []
    for data, ctrl in zip(datasets, control_datasets):        
        mean = calc_mean(data, var, basin)
        control = calc_mean(ctrl, var, basin)
        value = mean[0,:,:]-control[0,:,:]
        values.append(value)
        if contour == "density":
            dens = calculate_mean_density(data, basin)[0,:,:]
            density.append(dens)

    if cmap=="jet":
        contour_color = "white"
    else:
        contour_color = "black"

    if contour == "density":
        datasets = zip(values, density)
    else:
        datasets = values #just awful find a better way to do this

    fig, axs = plt.subplots(1, 3, figsize=figsize, sharey = True, sharex =True)
    
    lat=tavgs[0]["latitude"][:xextent]
    depth=tavgs[0]["depth"]
    x, y = np.meshgrid(lat, depth)
    axs[0].invert_yaxis()
    i = 0
    for ax, dataset_pair in zip(axs.flat, datasets):
        ax.set_xlim(xlim)
        value = dataset_pair[0]
        if contour == "density":
            levs = np.linspace(27.4, 28, 13)#np.concatenate([np.linspace(22, 27, 6),np.linspace(27.7, 28, 7)])
            c = ax.pcolormesh(x, y,value[:,:xextent], cmap=cmap, vmin=vmin, vmax=vmax) # i don't know ehy the number of indices is different but not looking into it now
            c1 = ax.contour(x,y,dataset_pair[1][:,:xextent],levels=levs, colors=[contour_color])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)

        elif contour:
            c = ax.pcolormesh(x, y,value[:,:xextent], cmap=cmap, vmin=vmin, vmax=vmax)
            c1 = ax.contour(x,y,value[:,:xextent],levels=np.linspace(vmin, vmax, 21), colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
        else:
            c = ax.pcolormesh(x, y,value[:,:xextent], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(subset["ice"][i] + ", " + subset["orb"][i] + " - " + control_subset["orb"][i])
        i+=1 

    fig.suptitle(title, fontsize=18, y = 1)

    fig.text(0.5, 0.167, 'Latitude', ha='center', fontsize= 14)
    fig.text(0.05, 0.55, 'Depth (m)', va='center', rotation='vertical', fontsize= 14)  
    
    cb = fig.colorbar(axs[0].collections[0], ax=axs, orientation='horizontal', fraction=0.04, pad=0.15)
    cb.ax.tick_params(labelsize=12)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize=12)

# %%
transect_ver2("O_temp", "Global", "coolwarm", -1,1,contour = "density", extent ="Global",  title="")
surface_ver2("O_temp", "coolwarm", -1, 1,  nbound = -45, contour = None, levels = None,  title="SST")
surface_ver2("O_temp", "coolwarm", -1, 1, proj = ccrs.PlateCarree(), nbound = 90, contour = None, levels = None,  title="SST")
surface_ver2("O_temp", "coolwarm", -1, 1, proj = ccrs.NorthPolarStereo())

# %%
transect_ver2("O_sal", "Global", "PiYG", -0.1,0.1,contour = "density", extent ="Global",  title="")
surface_ver2("O_sal", "PiYG", -0.5,0.5,  nbound = -45, contour = None, levels = None,  title="")
surface_ver2("O_sal", "PiYG", -0.5,0.5, proj = ccrs.PlateCarree(), nbound = 90, contour = None, levels = None,  title="")

# %%
transect_ver2("O_temp", "Global", "coolwarm", -1,1,contour = "density", extent ="Global",  title="")
transect_ver2("O_temp", "Atlantic", "coolwarm", -1,1,contour = "density", extent ="SO",  title="Atlantic")
transect_ver2("O_temp", "Pacific", "coolwarm", -1,1,contour = "density", extent ="SO",  title="Pacific")
transect_ver2("O_temp", "Indian", "coolwarm", -1,1,contour = "density", extent ="SO",  title="Indian")
# %%
transect_ver2("O_sal", "Global", "PiYG", -0.1,0.1, contour ="density", extent ="Global",  title="")
transect_ver2("O_sal", "Atlantic", "PiYG", -0.1,0.1, contour ="density", extent ="SO",  title="Atlantic")
transect_ver2("O_sal", "Pacific", "PiYG", -0.1,0.1, contour ="density", extent ="SO",  title="Pacific")
transect_ver2("O_sal", "Indian", "PiYG", -0.1,0.1, contour ="density", extent ="SO",  title="Indian")
# %%
surface_ver2("A_sat", "coolwarm", -1, 1,proj = ccrs.PlateCarree(), title="")
surface_ver2("A_sat", "coolwarm", -1, 1,proj = ccrs.SouthPolarStereo(), title="")

# %%
surface_ver2("A_awindY", "managua", -0.1, 0.1,proj = ccrs.SouthPolarStereo(), title="")
surface_ver2("A_awindX", "managua", -0.1, 0.1,proj = ccrs.SouthPolarStereo(), title="")
surface_ver2("O_tauX", "managua", -0.005, 0.005,proj = ccrs.SouthPolarStereo(), title="")

#%%
def plot_ovt(basin, title, save=False, path = None, datasets = datasets, control_datasets=control_datasets):
    levs=[-15, -10, -5, 0, 5, 10, 15]         
    vmin = -0.8
    vmax = 0.8
    levels = [i * 0.5 for i in levs]
    fig, axs = plt.subplots(4,1, figsize=(8,8), sharex =True, sharey =True, layout='constrained')
    fig.suptitle(title, fontsize=16)

    latitude=tavgs[0]["latitude"]
    depth=tavgs[1]["depth_edges"]#depth_edges
    axs[0].invert_yaxis()
    x, y = np.meshgrid(latitude, depth)
    for i, ax in enumerate(axs.flat):
        ax.invert_yaxis()
        ax.set_xlim((-90,90))
        picontrol = calc_ovt(tavgs[0], region=basin)
        if i==0:
            ax.pcolormesh(x, y, picontrol[0,:,:],
                        vmin=-50, vmax=50, 
                        cmap="seismic", shading="auto")
            c1 = ax.contour(x,y, picontrol[0,:,:],levels=levs, colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
            ax.set_title(ref_df["ice"][i])
            cb = fig.colorbar(ax.collections[0], ax=ax, orientation='vertical', shrink=1)
            cb.set_label("Meridional overturning (Sv)", fontsize=8)
            ax.set_title("PRD-PI")
            ax.set_ylabel("Depth")
        else:
            ovt = calc_ovt(datasets[i-1],region=basin)
            control = calc_ovt(control_datasets[i-1],region=basin)
            value = ovt.values - control.values
            ax.pcolormesh(x, y, value[0,:,:],
                            vmin=vmin, vmax=vmax, 
                            cmap="RdYlGn_r", shading="auto")
            c1 = ax.contour(x,y, value[0,:,:],levels=levels, colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
            ax.set_title(subset["ice"][i-1] + ", " + subset["orb"][i-1] + " - " + control_subset["orb"][i-1], fontsize=10)
            ax.set_ylabel("Depth")
    
    cb = fig.colorbar(axs[1].collections[0], ax=axs, orientation='horizontal', fraction=0.02, pad=0.07)
    cb.set_label("Meridional overturning (Sv)", fontsize=14)
    fig.text(0.5, 0.1, 'Latitude', ha='center', fontsize= 14)

    if save==True:
        fig.savefig(path + title.replace(" ", "_") + ".png", transparent=True)

plot_ovt("Atlantic", "AMOC", save=False, path = None)
# %%
