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
config_names =["PRD-PI", "PRD-LIG", "AMN-PI","AMN-LIG", "EQB-PI",  "EQB-LIG"] 
co2 = ["280"]
tavgs, ref_df = import_data(path_to_data, co2, config_names)

#%%
def surface_ver1(var, cmap, vmin, vmax, lev = 0, proj = ccrs.SouthPolarStereo(), nbound = -45, contour = None, levels = None,  title=""): 
    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = nbound
        figsize=(5,8)
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 40
        north_bound = 90
        figsize=(5,8)
    elif proj == ccrs.PlateCarree():
        south_bound = -90
        north_bound = 90
        figsize=(10,8)
    
    datasets = [get_2d_slice(tavgs[i][var], lev) for i in range(len(tavgs))]
    suptitle = title 
    #

    fig, axs = plt.subplots(3, 2, figsize=figsize, subplot_kw={'projection': proj})
    fig.suptitle(suptitle, fontsize=18, y=0.95)

    lons = tavgs[0]["longitude"]
    lats = tavgs[0]["latitude"]

    for i, ax in enumerate(axs.flat):
        data = datasets[i]
        ax.grid(visible=False)
        ax.set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(),
                        vmin=vmin, vmax=vmax, 
                        cmap=cmap, shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
        if contour:
            c1 = ax.contour(lons, lats,contour[i], colors=["black"], levels = levels, transform=ccrs.PlateCarree())
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)

    for n, col in enumerate(ref_df["orb"].unique()):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(ref_df["ice"].unique()):
        fig.text(0.12, 0.78 - (n * 0.25), row, va='center', ha='right', fontsize=12, rotation=90)

    cb = fig.colorbar(axs[0,0].collections[1], ax=axs, orientation='horizontal', fraction=0.02, pad=0.04)
    cb.set_label(f"{tavgs[0][var].attrs.get('long_name', var)} ({tavgs[0][var].attrs.get('units', '')})", fontsize=10)

surface_ver1("O_temp", "Spectral_r", -2, 20, lev = 0, proj = ccrs.PlateCarree(), nbound = 90, contour = None, levels = None,  title="SST")
surface_ver1("O_temp", "Spectral_r", -2, 20, lev = 0, proj = ccrs.SouthPolarStereo(), nbound =  -45, contour = None, levels = None,  title="SST")
surface_ver1("O_temp", "Spectral_r", -2, 20, lev = 0, proj = ccrs.NorthPolarStereo(), nbound =  90, contour = None, levels = None,  title="SST")

#%%
densities = []
for tavg in tavgs:
    density = calculate_density(tavg)
    surf_density = density[0,0,:,:]
    densities.append(surf_density)
#%%
surface_ver1("F_heat", "jet", -100, 100, proj = ccrs.SouthPolarStereo(), nbound =-60, contour = densities, levels=np.linspace(27.4,  28, 7),  title="heat flux+density")
surface_ver1("O_tauX", "managua", -0.1, 0.1,proj = ccrs.SouthPolarStereo(), title="")

#%%
def transect_ver1(var, basin, cmap, vmin,vmax,contour = None, extent="Global",  title=""):
    if extent == "SO":
        xlim=(-90, -50)
        xextent = 24
        figsize=(8, 12)
    else:
        xlim=(-90, 90)
        xextent = 100
        figsize=(10, 12)
    values = []
    density = []
    for tavg in tavgs:        
        mean = calc_mean(tavg,var, basin)
        #mean = tavg[var][0,:,:,:]
        value=mean.values 
        values.append(value)
        if contour == "density":
            dens = calculate_mean_density(tavg, basin)[0,:,:]
            #dens = calculate_density(tavg)[0,:,:,:]
            density.append(dens)

    if cmap=="jet":
        contour_color = "white"
    else:
        contour_color = "black"

    if contour == "density":
        datasets = zip(values, density)
    else:
        datasets = values #just awful find a better way to do this

    fig, axs = plt.subplots(3, 2, figsize=figsize, sharey = True, sharex =True)
    
    lat=tavgs[0]["latitude"][:xextent]
    depth=tavgs[0]["depth"]
    x, y = np.meshgrid(lat, depth)
    axs[0,0].invert_yaxis()
    i = 0
    for ax, dataset_pair in zip(axs.flat, datasets):
        ax.set_xlim(xlim)
        
        value = dataset_pair[0]
        if contour == "density":
            levs = np.linspace(27.4, 28, 7)#np.concatenate([np.linspace(22, 27, 6),np.linspace(27.7, 28, 7)])
            c = ax.pcolormesh(x, y,value[0, :,:xextent], cmap=cmap, vmin=vmin, vmax=vmax) # i don't know ehy the number of indices is different but not looking into it now
            c1 = ax.contour(x,y,dataset_pair[1][:,:xextent],levels=levs, colors=[contour_color])

            ax.clabel(c1, c1.levels, inline=True, fontsize=10)

        elif contour:
            c = ax.pcolormesh(x, y,value[:,:xextent], cmap=cmap, vmin=vmin, vmax=vmax)
            c1 = ax.contour(x,y,value[:,:xextent],levels=np.linspace(vmin, vmax, 21), colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
        else:
            c = ax.pcolormesh(x, y,value[:,:xextent], cmap=cmap, vmin=vmin, vmax=vmax)
        i+=1 
        

    fig.suptitle(title, fontsize=18, y = 0.93)

    for n, col in enumerate(ref_df["orb"].unique()):
       axs[0, n].set_title(col, fontsize=14)
       
    for n, row in enumerate(ref_df["ice"].unique()):
        fig.text(0.04, 0.78 - (n * 0.25), row, va='center', ha='right', fontsize=14, rotation=90)

    fig.text(0.5, 0.13, 'Latitude', ha='center', fontsize= 12)
    fig.text(0.05, 0.55, 'Depth (m)', va='center', rotation='vertical', fontsize= 12)  
    
    cb = fig.colorbar(axs[0,0].collections[0], ax=axs, orientation='horizontal', fraction=0.02, pad=0.06)
    cb.ax.tick_params(labelsize=12)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize=12)

#%%
transect_ver1("O_temp", "Global", "Spectral_r", -2,20,contour = "density", extent ="Global",  title="Global zonal mean")
transect_ver1("O_temp", "Global", "Spectral_r", -2,20,contour = "density", extent ="SO",  title="SO zonal mean")

transect_ver1("O_temp", "Atlantic", "Spectral_r", -2,10,contour = "density", extent ="SO",  title="Atlantic mean")
transect_ver1("O_temp", "Pacific", "Spectral_r", -2,10,contour = "density", extent ="SO",  title="Pacific mean")
transect_ver1("O_temp", "Indian", "Spectral_r", -2,10,contour = "density", extent ="SO",  title="Indian mean")

# %%
def plot_ovt(basin, title, save=False, path = None):
    levels=[-15, -10, -5, 0, 5, 10, 15]         
    cmap = "seismic"
    vmin = -50
    vmax = 50
  
    fig, axs = plt.subplots(3,2, figsize=(12,8), sharex =True, sharey =True)#(4, 5, figsize=(20,12)
    fig.suptitle(title, fontsize=16, y=0.95)

    latitude=tavgs[0]["latitude"]
    depth=tavgs[1]["depth_edges"]#depth_edges
    axs[0,0].invert_yaxis()
    x, y = np.meshgrid(latitude, depth)
    for i, ax in enumerate(axs.flat):

        ax.set_xlim((-90,90))
        ovt = calc_ovt(tavgs[i],region=basin)
        ax.pcolormesh(x, y, ovt[0,:,:],
                        vmin=vmin, vmax=vmax, 
                        cmap=cmap, shading="auto")
        c1 = ax.contour(x,y, ovt[0,:,:],levels=levels, colors=['black'])
        ax.clabel(c1, c1.levels, inline=True, fontsize=10)
   
    for n, col in enumerate(ref_df["orb"].unique()):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(ref_df["ice"].unique()):
        fig.text(0.05, 0.78 - (n * 0.25), row, va='center', ha='right', fontsize=12, rotation=90)    
    cb = fig.colorbar(axs[0, 0].collections[0], ax=axs, orientation='horizontal', fraction=0.02, pad=0.07)

    cb.set_label("Meridional overturning (Sv)", fontsize=14)
    
    fig.text(0.5, 0.13, 'Latitude', ha='center', fontsize= 14)
    fig.text(0.06, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize= 14) 

    if save==True:
        fig.savefig(path + title.replace(" ", "_") + ".png", transparent=True)

plot_ovt("Atlantic", "AMOC", save=False, path = None)
