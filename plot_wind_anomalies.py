#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:12:22 2024

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

savepath = "/Users/mackenfr/Desktop/"
#%%  
def plot_polar_grid_winds(var,cmap,vmin,vmax, title):
    anoms = []
    i=0
    for tavg in tavgs:

        if i%2 == 0:
            anom=tavgs[i+1][var].values - tavgs[i][var].values #wind minus default
            anoms.append(anom)
        i +=1
        
    fig, axs = plt.subplots(4, 5, figsize=(15, 12), subplot_kw={'projection': ccrs.SouthPolarStereo()})
    fig.suptitle(title, fontsize="xx-large", y = 0.93)
    for ax, dataset in zip(axs.flat, anoms):
        #import pdb; pdb.set_trace()

        if len(dataset.shape) == 4:
            data = dataset[0,0,:,:]
        if len(dataset.shape) == 3:
            data = dataset[0,:,:]
        #import pdb; pdb.set_trace()
        ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
        ax.coastlines() 
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ct = ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"], data,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")
    
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(configs):
        #axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.78 - (n * 0.19), row, va='center', ha='right', fontsize=12, rotation=0)
        
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='vertical', fraction=0.05, pad=0.04)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize="large")


def plot_global_grid_winds(var,cmap,vmin,vmax, title):    
    anoms = []
    i=0
    for tavg in tavgs:
        if i%2 == 0:
            anom=tavgs[i+1][var].values - tavgs[i][var].values
            anoms.append(anom)
        i +=1
    fig, axs = plt.subplots(4, 5, figsize=(20, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.suptitle(title, fontsize="xx-large", y = 0.93)
    for ax, dataset in zip(axs.flat, anoms):
           
        if len(dataset.shape) == 4:
            data = dataset[0,0,:,:]
        if len(dataset.shape) == 3:
            data = dataset[0,:,:]

        ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        ax.coastlines() 
        #ax.patch.set_alpha(0)   
        ax.set_facecolor('white')
        #ax.gridlines(linestyle='--')
    
        ct=ax.pcolormesh(tavgs[0]["longitude"], tavgs[0]["latitude"],data,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap)#, shading="auto")
    
        plt.title(title)
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col)
       
    for n, row in enumerate(configs):
        #axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.78 - (n * 0.19), row, va='center', ha='right', fontsize=12, rotation=0)
        
    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='vertical', fraction=0.05, pad=0.04)
    cb.set_label(tavgs[0][var].attrs["long_name"] + " (" + tavgs[0][var].attrs["units"] + ")", fontsize="large")

def plot_ovt_grid_wind(data, basin, title):
    fig, axs = plt.subplots(4, 5, figsize=(25, 13), sharey = True, sharex =True)
    fig.suptitle(title, fontsize=24, y = 0.93)
    axs[0,0].invert_yaxis()
    levels=[-15, -10, -5, 0, 5, 10, 15]
    cmap="PuOr"
    vmin = -30
    vmax = 30
    lat=data[0]["latitude_V"]
    depth=data[0]["depth_edges"]
    x, y = np.meshgrid(lat, depth)
    anoms = []
    i = 0
    for tavg in tavgs:
        if i%2 == 0:
            ovt = calc_ovt(tavg, basin)
            ovt1 = calc_ovt(tavgs[i+1], basin)
            anom =  ovt1.values - ovt.values #wind - defaul
            anoms.append(anom)
            #import pdb; pdb.set_trace()
        i +=1
    for ax, dataset in zip(axs.flat, anoms):
        ax.invert_yaxis()
               
        c = ax.pcolormesh(x, y,dataset[0,:,:], cmap=cmap, vmin=vmin, vmax=vmax)
        c1 = ax.contour(x,y,dataset[0,:,:],levels=levels, colors=['black'])
        ax.clabel(c1, c1.levels, inline=True, fontsize=10)
    fig.text(0.5, 0.09, 'Latitude', ha='center', fontsize= 14)
    fig.text(0.1, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize= 14)  
    for n, col in enumerate(co2ccn):
       axs[0, n].set_title(col, fontsize=18)
       
    for n, row in enumerate(configs):
        #axs[0, n].set_ylabel(row)
        fig.text(0.085, 0.78 - (n * 0.19), row, va='center', ha='right', fontsize=18, rotation=90)
        
    cb = fig.colorbar(axs[0, 4].collections[0], ax=axs, orientation='vertical', fraction=0.025, pad=0.01)
    cb.set_label("Overturning (Sv)", fontsize=18)
    cb.ax.tick_params(labelsize=20)