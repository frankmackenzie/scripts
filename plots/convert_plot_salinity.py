import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import colormaps
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
import gsw
os.chdir("/Users/mackenfr/scripts")
from import_data import import_data
#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
tavgs, ref_df = import_data(path_to_data, co2s=["280"], config_names=["PRD-PI", "AMN-PI", "EQB-PI"])
import pandas as pd
#%%
pressure = tavgs[0]["depth"]
lons = tavgs[0]["longitude"]
lats = tavgs[0]["latitude"]
#%%
def plot_polar_grid(proj = ccrs.SouthPolarStereo(), nbound = -45, contour = None, levels = None,  title=""): 
    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = nbound
        figsize=(8,8)
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 45
        north_bound = 90
        figsize=(8,8)
    elif proj == ccrs.PlateCarree():
        # south_bound = -90
        # north_bound = 90
        # figsize=(15,8)
        south_bound = -90
        north_bound = 90
        figsize=(10,6)
        
    suptitle = title 
    #

    fig, axs = plt.subplots(2, 2, figsize=figsize, subplot_kw={'projection': proj})
    fig.suptitle(suptitle, fontsize=18, y=0.93)



    for i, ax in enumerate(axs.flat):
        sal = tavgs[i]["O_sal"][0, 0, :, :]
        SA = gsw.conversions.SA_from_SP(sal,0, lons, lats)
        ax.set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ax.pcolormesh(lons, lats, SA, transform=ccrs.PlateCarree(),
                        vmin=32, vmax=37, 
                        cmap="gist_ncar", shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
        if contour:
            c1 = ax.contour(lons, lats,contour[i+1], colors=["black"], levels = levels, transform=ccrs.PlateCarree())
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
        ax.set_title(ref_df["ice"][i])
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label("Salinity (g/kg))", fontsize=14)


#%%
def plot_polar_grid_anom(cmap, vmin, vmax, nbound = -45,  lev = 0, proj = ccrs.SouthPolarStereo(), contour = None, levels = None,  title=""): 
    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = nbound
        figsize=(8,3)
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 45
        north_bound = 90
        figsize=(8,3)
    elif proj == ccrs.PlateCarree():
        # south_bound = -90
        # north_bound = 90
        # figsize=(15,8)
        south_bound = -90
        north_bound = 90
        figsize=(10,2)
        
    suptitle = title + " anom from PRD"
    #

    fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': proj})
    fig.suptitle(suptitle, fontsize=18, y=0.99)

    lons = tavgs[0]["longitude"]
    lats = tavgs[0]["latitude"]
    refsal = tavgs[0]["O_sal"][0, 0, :, :]
    refSA = gsw.conversions.SA_from_SP(refsal,0, lons, lats)          
    for i, ax in enumerate(axs.flat):
        sal = tavgs[i+1]["O_sal"][0, 0, :, :]
        SA = gsw.conversions.SA_from_SP(sal,0, lons, lats)    
        data=SA-refSA        
        ax.set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(),
                        vmin=vmin, vmax=vmax, 
                        cmap=cmap, shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
        if contour:
            c1 = ax.contour(lons, lats,contour[i], colors=["black"], levels = levels, transform=ccrs.PlateCarree())
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
        ax.set_title(ref_df["ice"][i+1])
    cb = fig.colorbar(axs[0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
#%%
def calc_mean(data, basin):
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    G_mskhr = tavgs[0]["G_mskhr"]
    if basin=="Global":
        mean = data[:,:,:].mean(dim=["longitude"])
    elif basin in basins:
        mask = data.where(G_mskhr==basins[basin])
        mean = mask[:,:,:].mean(dim=["longitude"])
    return(mean)

def plot_transect_grid(basin, cmap, vmin,vmax,contour = None, xlim=(-90, 90),  title=""):
    fig, axs = plt.subplots(1, 3, figsize=(10, 6), sharey = True, sharex =True)
    axs[0].invert_yaxis()
    lat=tavgs[0]["latitude"][:24]
    depth=tavgs[0]["depth"]
    x, y = np.meshgrid(lat, depth)
    suptitle=title

    values = []
    density = []

    for tavg in tavgs:        
        sal = tavg["O_sal"][0, :, :, :]
        SA = gsw.conversions.SA_from_SP(sal, pressure, lons, lats)
        #mean=calc_mean(SA, basin)
        transect = SA[:,:,72]
        values.append(transect)

    if cmap=="jet":
        contour_color = "white"
    else:
        contour_color = "black"

    if contour == "density":
        datasets = zip(values, density)
    else:
        datasets = values #just awful find a better way to do this
    i = 0
    for ax, dataset_pair in zip(axs.flat, datasets):
        #ax.invert_yaxis()
        ax.set_xlim(xlim)
        ax.set_ylim(2000)
        value = dataset_pair
        if contour == "density":
            levs = np.linspace(27.4, 28, 7)
            c = ax.pcolormesh(x, y,value[:,:24], cmap=cmap, vmin=vmin, vmax=vmax) # i don't know ehy the number of indices is different but not looking into it now
            c1 = ax.contour(x,y,dataset_pair[1][0,:,:24],levels=levs, colors=[contour_color])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)

        elif contour:
            #return(value)
            c = ax.pcolormesh(x, y,value[:, :24], cmap=cmap, vmin=vmin, vmax=vmax)
            c1 = ax.contour(x,y,value[:,:24],levels=np.linspace(vmin, vmax, 21), colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
        else:
            c = ax.pcolormesh(x, y,value[:,:24], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(ref_df["ice"][i])
        i+=1 
        

    fig.suptitle(suptitle, fontsize=18, y = 0.98)

    fig.text(0.5, 0.16, 'Latitude', ha='center', fontsize= 14)
    fig.text(0.05, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize= 14)  
    
    cb = fig.colorbar(axs[0].collections[0], ax=axs, orientation='horizontal', fraction=0.04, pad=0.12)
    cb.ax.tick_params(labelsize=12)
    cb.set_label("Salintiy g/Kg", fontsize=12)
    #fig.savefig("/Users/mackenfr/Desktop/eric_meeting/" + title.replace(" ", "_") + ".png", transparent=True)
    
plot_transect_grid("Global", "gist_rainbow", 33,35,contour = True, xlim=(-90, -50),  title="")
#input_amn_sal_path = "/Users/mackenfr/scripts/pism_to_uvic/outputs/O_salsur_LIG.nc"
#input_amn_sal = xr.open_dataset(input_amn_sal_path)
#tavgs[0]=input_amn_sal
#let it be known that the fresh surface in the AMN is not due to initialisation
# %%
