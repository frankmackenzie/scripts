#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 13:51:51 2024

@author: mackenfr
"""

import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
os.chdir("/Users/mackenfr/scripts")

from calc_ovt import calc_ovt
#%% paths
path_to_data = "/Users/mackenfr/uvic_outputs/" 


#%%

def plot_polar(index, var, cmap,vmin,vmax, title, anom=True):
    tavg = tavgs[index][var]
    if anom:
        value=tavg[0:,:].values - tavgs[0][var][0,:,:].values
    else:
        value=tavg[0,:,:].values
    fig = plt.figure()
    
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
     # Limit the map to -45 degrees latitude and below.
    ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
    ax.coastlines() 
    ax.gridlines(linestyle='--', draw_labels=True, color="gray")
   
    ct=ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"], value,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")
   
    #colorbar 
    cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)
    cb.set_label(tavg.attrs["long_name"] + " (" + tavg.attrs["units"] + ")")
   
    
    plt.title(title)
    plt.tight_layout()
    return(fig)
    
    
def plot_polar_grid(var, cmap,vmin,vmax, title="", anom=None): 
    values = []
    control=0
    index_co2=0
    index_config=1
    for tavg in tavgs:
        if anom == "picontrol":
            suptitle=title + " anom from PRD-PI280"
            value=tavg[var].values - tavgs[0][var].values
        elif anom == "co2":
            suptitle=title + " anom from 280ppm"
            value=tavg[var].values - tavgs[control][var].values
            index_co2+=1
            if index_co2!=0 and index_co2%5==0:
                control+=5
        elif anom == "config":
            suptitle = title + " anom from PRD-PI"
            if index_config%5 ==0:
                control = 4
            else:
                control=(index_config%5)-1
            #import pdb; pdb.set_trace()
            value=tavg[var].values - tavgs[control][var].values
            index_config+=1
        else:
            value=tavg[var].values 
        values.append(value)        
    fig, axs = plt.subplots(4, 5, figsize=(15, 12), subplot_kw={'projection': ccrs.SouthPolarStereo()})
    fig.suptitle(suptitle, fontsize="xx-large", y = 0.93)
    for ax, dataset in zip(axs.flat, values):
        #import pdb; pdb.set_trace()
        if len(dataset.shape) == 4:
            data = dataset[0,0,:,:]
        if len(dataset.shape) == 3:
            data = dataset[0,:,:]
        ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
        ax.coastlines() 
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")

        ct = ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"], data,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")
    
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(configs):
        axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.8 - (n * 0.18), row, va='center', ha='right', fontsize=12, rotation=0)
        
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    return(fig)

def plot_polar_grid_depth(var, cmap,vmin,vmax, title="", anom=None): 
    values = []
    control=0
    index_co2=0
    index_config=1
    for tavg in tavgs:
        if anom == "picontrol":
            suptitle=title + " anom from PRD-PI280"
            value=tavg[var].values - tavgs[0][var].values
        elif anom == "co2":
            suptitle=title + " anom from 280ppm"
            value=tavg[var].values - tavgs[control][var].values
            index_co2+=1
            if index_co2!=0 and index_co2%5==0:
                control+=5
        elif anom == "config":
            suptitle = title + " anom from PRD-PI"
            if index_config%5 ==0:
                control = 4
            else:
                control=(index_config%5)-1
            value=tavg[var].values - tavgs[control][var].values
            index_config+=1
        else:
            value=tavg[var].values
        values.append(value)      
    fig, axs = plt.subplots(4, 5, figsize=(15, 12), subplot_kw={'projection': ccrs.SouthPolarStereo()})
    fig.suptitle(title, fontsize="xx-large", y = 0.93)
    for ax, dataset in zip(axs.flat, values):
        if len(dataset.shape) == 4:
            data = dataset[0,depth,:,:]
        if len(dataset.shape) == 3:
            data = dataset[0,:,:]
        ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
        ax.coastlines() 
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        #import pdb; pdb.set_trace()

        ct = ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"], data,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")
    
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(configs):
        axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.8 - (n * 0.18), row, va='center', ha='right', fontsize=12, rotation=0)

        
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    return(fig)

def plot_global_grid(var, cmap,vmin,vmax, title="", anom=None):    
    values = []
    control=0
    index_co2=0
    index_config=1
    for tavg in tavgs:
        if anom == "picontrol":
            suptitle=title + " anom from PRD-PI280"
            value=tavg[var].values - tavgs[0][var].values
        elif anom == "co2":
            suptitle=title + " anom from 280ppm"
            value=tavg[var].values - tavgs[control][var].values
            index_co2+=1
            if index_co2!=0 and index_co2%5==0:
                control+=5
        elif anom == "config":
            suptitle = title + " anom from PRD-PI"
            if index_config%5 ==0:
                control = 4
            else:
                control=(index_config%5)-1
            value=tavg[var].values - tavgs[control][var].values
            index_config+=1
        else:
            value=tavg[var].values[0,:,:]
            suptitle=title
        values.append(value) 
    #import pdb;pdb.set_trace()
    fig, axs = plt.subplots(4, 5, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.suptitle(suptitle, fontsize="xx-large", y = 0.93)
    for ax, dataset in zip(axs.flat, values):
           
        if len(dataset.shape) == 4:
            data = dataset[0,0,:,:]
        if len(dataset.shape) == 3:
            data = dataset[0,:,:]
        else:
            data = dataset[:,:]

        ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        ax.coastlines() 
        #ax.patch.set_alpha(0)   
        ax.set_facecolor('white')
        #ax.gridlines(linestyle='--')
    
        ct=ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"],data,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap)#, shading="auto")
    
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(configs):
        axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.8 - (n * 0.18), row, va='center', ha='right', fontsize=12, rotation=0)

        
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    return(fig)

def plot_NA_grid(var, cmap,vmin,vmax, title="", anom=None):    
    values = []
    control=0
    index_co2=0
    index_config=1
    for tavg in tavgs:
        if anom == "picontrol":
            suptitle=title + " anom from PRD-PI280"
            value=tavg[var].values - tavgs[0][var].values
        elif anom == "co2":
            suptitle=title + " anom from 280ppm"
            value=tavg[var].values - tavgs[control][var].values
            index_co2+=1
            if index_co2!=0 and index_co2%5==0:
                control+=5
        elif anom == "config":
            suptitle = title + " anom from PRD-PI"
            if index_config%5 ==0:
                control = 4
            else:
                control=(index_config%5)-1
            value=tavg[var].values - tavgs[control][var].values
            index_config+=1
        else:
            value=tavg[var].values
        values.append(value) 
    fig, axs = plt.subplots(4, 5, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.suptitle(suptitle, fontsize="xx-large", y = 0.93)
    for ax, dataset in zip(axs.flat, values):
           
        if len(dataset.shape) == 4:
            data = dataset[0,0,:,:]
        if len(dataset.shape) == 3:
            data = dataset[0,:,:]

        ax.set_extent([-70, 25, 35, 85], ccrs.PlateCarree())
        ax.coastlines() 
        #ax.patch.set_alpha(0)   
        ax.set_facecolor('white')
        #ax.gridlines(linestyle='--')
    
        ct=ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"],data,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap)#, shading="auto")
    
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(configs):
        axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.8 - (n * 0.18), row, va='center', ha='right', fontsize=12, rotation=0)

        
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    return(fig)

def plot_ovt_grid(data, basin, title="", anom=None):
    fig, axs = plt.subplots(4, 5, figsize=(30, 15), sharey = True, sharex =True)
    
    #axs[0,0].invert_yaxis()
    levels=[-30, -15, 0, 5, 15, 30]#[-15, -10, -5, 0, 5, 10, 15]
    cmap="PuOr_r"
    vmin = -20
    vmax = 20
    lat=data[0]["latitude_V"][:]#[9:24] for SO
    depth=data[0]["depth_edges"]
    x, y = np.meshgrid(lat, depth)
    
    values = []
    control=0
    index_co2=0
    index_config=1
    for tavg in tavgs:
        ovt = calc_ovt(tavg, basin)
        picontrol = calc_ovt(tavgs[0], basin)
        if anom == "picontrol":
            suptitle=title + " anom from PRD-PI280"
            value=ovt.values - picontrol.values
        elif anom == "co2":
            control_ovt = calc_ovt(tavgs[control], basin)
            suptitle=title + " anom from 280ppm"
            value= ovt.values - control_ovt.values
            index_co2+=1
            if index_co2!=0 and index_co2%5==0:
                control+=5
        elif anom == "config":
            control_ovt = calc_ovt(tavgs[control], basin)
            suptitle = title + " anom from PRD-PI"
            if index_config%5 ==0:
                control = 4
            else:
                control=(index_config%5)-1
            value=ovt.values - control_ovt.values
            index_config+=1
        else:
            value=ovt
        values.append(value) 
        #import pdb; pdb.set_trace()
    fig.suptitle(basin + " " + suptitle, fontsize=24, y = 0.93)
    if anom=="picontrol":
        axs[0,0].pcolormesh(x, y,picontrol[0,:,:], cmap="bwr", vmin=-30, vmax=30)
        c1 = axs[0,0].contour(x,y,picontrol[0,:,:],levels=levels, colors=['black'])
        axs[0,0].clabel(c1, c1.levels, inline=True, fontsize=10)
        fig.colorbar(axs[0, 0].collections[0], ax=axs[0,0], orientation='vertical',fraction=0.03, pad=0.025)
        
        for ax, dataset in zip(axs.flat[1:], values[1:]):
            ax.invert_yaxis()
            c = ax.pcolormesh(x, y,dataset[0,:,:], cmap=cmap, vmin=vmin, vmax=vmax)
            c1 = ax.contour(x,y,picontrol[0,:,:],levels=levels, colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)
#adapt this to other anom options later
       
    else:
        for ax, dataset in zip(axs.flat[1:], values[1:]):
            ax.invert_yaxis()
            c = ax.pcolormesh(x, y,dataset[0,:,:], cmap=cmap, vmin=vmin, vmax=vmax)
            c1 = ax.contour(x,y,dataset[0,:,:],levels=levels, colors=['black'])
            ax.clabel(c1, c1.levels, inline=True, fontsize=10)

        
    # fig.text(0.5, 0.17, 'Latitude', ha='center', fontsize= 18)
    # fig.text(0.09, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize= 18)  
    fig.text(0.5, 0.16, 'Latitude', ha='center', fontsize= 18)
    fig.text(0.06, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize= 18)  
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col, fontsize=18)
       
    for n, row in enumerate(configs):
        axs[0, n].set_ylabel(row)
        fig.text(0.08, 0.8 - (n * 0.18), row, va='center', ha='right', fontsize=18, rotation=90)
        
    cb = fig.colorbar(axs[0, 4].collections[0], ax=axs, orientation='horizontal', fraction=0.01, pad=0.08)
    cb.set_label("Overturning (Sv)", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    return(fig)
       

#%% load data
configs = ["PRD-PI","AMN-PI"]# "PRD-LIG", "AMN-PI", "AMN-LIG", "TRS-LIG", "TRS-PI", "EQB-PI", "EQB-LIG"]
co2ccn = ["280", "300", "350", "400", "450"]
tavgs = []
si_maxs = []
si_mins = []
names = []
co2_all = []
for config in configs:
    for co2 in co2ccn:
        filename_tavg = path_to_data + "/" + config + "/" + co2 + "/tavg.nc"

        name = config 
        data = xr.open_dataset(filename_tavg, decode_times=False, engine = "scipy")
        names.append(name)
        tavgs.append(data)
        co2_all.append(co2)

        
#%% 
plot_global_grid("O_temp", "coolwarm",-5,5, "SST anomaly from PI")

plot_polar_grid( "O_temp", "coolwarm",-5,5, anom="picontrol,", "SST anomaly from PI")

plot_polar_grid("A_sat", "RdYlBu_r",-10,10, "SAT anomaly from PI ")

plot_ovt_grid(tavgs,  "Indian", "Indian ocean overturning anomaly")

plot_ovt_grid(tavgs,  "Atlantic", "Atlantic ocean overturning anomaly")


plot_polar_grid("F_heat", "Spectral", -70,70, "Heat flux")

plot_polar_grid("O_ventdep", "PiYG", -100,100, "Ventdepth anomaly")

plot_polar_grid("G_kmt", "autumn", 0,19, " ")

plot_polar_grid("O_velX", "PRGn", -0.1,0.1, " ")

plot_polar_grid("O_velX", "PRGn", -0.05,0.05, " Ocean vel anom")

plot_polar_grid("O_velY", "PiYG", -0.02,0.02, " Ocean velY anom")


plot_polar_grid("A_windqX", "PuOr", -0.5,0.5, " awind moisture anom")

plot_global_grid_anom("A_windqX", "PuOr", -0.5,0.5, " awind moisture anom")

#%%

#%%
for depth in list(range(0,10)):
    plot_polar_grid_depth(depth, "O_velX", "PRGn", -0.05,0.05, " Ocean vel anom")


#%%
plot_global_grid("A_windspd", "PiYG",-2,2,anom="bonus")
#%%
#heat flux
savepath = "/Users/mackenfr/plots_organised/"
#%%
path = savepath+"sst/"
names = ["control_anom_global", "co2_anom_global", "config_anom_global"]#"PRD-PI_abs", "control_anom_polar", "co2_anom_polar",

plots=[
#(plot_polar(0,"O_temp", "coolwarm", -10,10, "SST PRD-PI", anom=False)),
#(plot_polar_grid("O_temp", "coolwarm", -4,4, "SST", "picontrol")),
#(plot_polar_grid("O_temp", "coolwarm", -4,4, "SST", "co2")),
#(plot_polar_grid("O_temp", "coolwarm", -4,4, "SST", "config")),
(plot_global_grid("O_temp", "coolwarm", -4,4, "SST", anom="picontrol")),
(plot_global_grid("O_temp", "coolwarm", -4,4, "SST", anom="co2")),
(plot_global_grid("O_temp", "coolwarm", -4,4, "SST", anom="config"))

]
#%%
for plot, name in zip(plots, names):
    plot.savefig(path + name+".png")
     
#%%
path = savepath+"heat_flux/"
names = ["PRD-PI_abs", "control_anom_polar", "co2_anom_polar","control_anom_global", "co2_anom_global", "config_anom_global"]#


#%%       
       
(plot_polar(0,"F_heat", "Spectral", -70,70, "Heat flux PRD-PI", anom=False)),
(plot_polar_grid("F_heat", "Spectral", -70, 70, "Heat flux", "picontrol")),
(plot_polar_grid("F_heat", "Spectral", -70, 70, "Heat flux", "co2")),
(plot_polar_grid("F_heat", "Spectral", -70, 70, "Heat flux", "config")),
(plot_global_grid("F_heat", "Spectral", -70, 70, "Heat flux", anom="picontrol")),
(plot_global_grid("F_heat", "Spectral", -70, 70, "Heat flux", anom="co2")),
(plot_global_grid("F_heat", "Spectral", -70, 70, "Heat flux", anom="config"))
]

for plot, name in zip(plots, names):
    plot.savefig(path + name+".png")
    
#%%
plot_NA_grid("O_temp", "coolwarm", -4,4, "SST", anom="picontrol")
plot_NA_grid("O_temp", "coolwarm", -4,4, "SST", anom="co2")
plot_NA_grid("O_temp", "coolwarm", -4,4, "SST", anom="config")

#%%
plot_global_grid("A_albatm", "cubehelix",-0.005,0.005, anom ="co2")
plot_global_grid("A_albatm", "cubehelix",-0.005,0.005, anom ="config")

plot_global_grid("A_albsur", "cubehelix",-0.005,0.005, anom ="co2")
plot_global_grid("A_albsur", "cubehelix",-0.005,0.005, anom ="config")

plot_global_grid("A_albplt", "cubehelix",-0.005,0.005, anom ="picontrol")

plot_global_grid("A_albplt", "cubehelix",-0.005,0.005, anom ="co2")
plot_global_grid("A_albplt", "cubehelix",-0.005,0.005, anom ="config")
plot_global_grid("A_albplt", "gnuplot",0,1, anom =None)
plot_global_grid("A_albsur", "gnuplot",0,1, anom =None)

plot_global_grid("A_albplt", "gnuplot",0,1, anom =None)


