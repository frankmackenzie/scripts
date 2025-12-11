#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:22:47 2024

@author: mackenfr
"""
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import numpy as np

#%% from xoverturning on github
def calc_ovt(data, region="Global"):
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    if region=="Global":
        #convert velocity to volume
        transport = data["O_velY"] * data["G_dxu"] * data["G_dzt"]
    elif region in basins:
        G_mskhr_aligned = data["G_mskhr"].sel(latitude=data["latitude_V"], longitude=data["longitude_V"], method="nearest")
        mask = data[["O_velY","G_dxu", "G_dzt"]].where(G_mskhr_aligned==basins[region])
        transport = (mask["O_velY"] * mask["G_dxu"] * mask["G_dzt"])

    zonalsum = transport.sum(dim="longitude_V")
    # integrate from bottom
    integ_layers_from_bottom = zonalsum.cumsum(dim="depth") - zonalsum.sum(dim="depth")
    # the result of the integration over layers is evaluated at the interfaces
    # with psi = 0 as the bottom boundary condition for the integration
    bottom_condition = xr.zeros_like(integ_layers_from_bottom.isel({"depth": -1}))
    psi_raw = xr.concat([integ_layers_from_bottom, bottom_condition], dim="depth")
    psi_raw = psi_raw.chunk({"depth": len(psi_raw["depth"])})  # need to rechunk to new size
    # convert kg.s-1 to Sv (1e6 m3.s-1)
    psi_Sv = psi_raw / 1.0e6
    return(psi_Sv)

#%% SO - swap lon lat        
def calc_ovt_so(data, region="Global"):
    mask = data[["O_velX","G_dyu", "G_dzt"]].where(data["latitude_V"]<-50)
    transport = (mask["O_velX"] * mask["G_dyu"] * mask["G_dzt"])

    meridsum = transport.sum(dim="latitude_V")
    # integrate from bottom
    integ_layers_from_bottom = meridsum.cumsum(dim="depth") - meridsum.sum(dim="depth")
    # the result of the integration over layers is evaluated at the interfaces
    # with psi = 0 as the bottom boundary condition for the integration
    bottom_condition = xr.zeros_like(integ_layers_from_bottom.isel({"depth": -1}))
    psi_raw = xr.concat([integ_layers_from_bottom, bottom_condition], dim="depth")
    #psi_raw = psi_raw.chunk({"depth": len(psi_raw["depth"])})  # need to rechunk to new size
    # convert kg.s-1 to Sv (1e6 m3.s-1)
    psi_Sv = psi_raw / 1.0e6
    return(psi_Sv)

def plot_ovt_SO(data, ovt, anom, title, save=False, path=None):
    levels=[-15, -10, -5, 0, 5, 10, 15]
    if anom == True:
        cmap="PiYG"
        vmin = -5
        vmax = 5
        levels = [i * 0.5 for i in levels]
    else:
        cmap="PuOr"
        vmin = -50
        vmax = 50
    lon=data["longitude_V"]
    depth=data["depth_edges"]

    x, y = np.meshgrid(lon, depth)

    fig, ax = plt.subplots()
    ax.set_xlabel("longitude")
    ax.set_ylabel("depth (m)")
    ax.invert_yaxis()
    
    ax.set_title(title)

    c = ax.pcolormesh(x, y,ovt[0,:,:], cmap=cmap, vmin=vmin, vmax=vmax)
    c1 = ax.contour(x,y,ovt[0,:,:],levels=levels, colors=['black'])
    ax.clabel(c1, c1.levels, inline=True, fontsize=10)

    cbar = plt.colorbar(c, label = "Zonal overturning Southern Ocean (Sv)")
    
    if save==True:
        fig.savefig(path + title.replace(" ", "_") + ".png", transparent=True)
    