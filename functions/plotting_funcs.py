#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 19:44:26 2023


@author: mackenfr
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import sys
import xarray as xr
import os

# In[]


def plot_polar_stereo(lon, lat, value,cmap,vmin,vmax, title, save=False, path=None):
    
    fig = plt.figure()

    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    # Limit the map to -45 degrees latitude and below.
    ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
    ax.coastlines() 
    ax.gridlines(linestyle='--', draw_labels=True, color="gray")

    ct=ax.pcolormesh(lon, lat, value,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")

    #colorbar 
    cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)
    cb.set_label(value.attrs["long_name"] + " (" + value.attrs["units"] + ")")

    
    plt.title(title)
    plt.tight_layout()
    if save==True:
        fig.savefig(path + "plots/" + title +'.png', transparent=True, dpi=300)
        
def plot_global(lon, lat, value,cmap,vmin,vmax, title, save=False, path=None):    
    
    fig = plt.figure()

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    ax.coastlines() 
    #ax.patch.set_alpha(0)   
    ax.set_facecolor('white')
    #ax.gridlines(linestyle='--')

    ct=ax.pcolormesh(lon,lat,value,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap)#, shading="auto")

    #colorbar 
    cb = plt.colorbar(ct,orientation="vertical",fraction=0.023, pad=0.04)
    cb.set_label(value.attrs["long_name"] + " (" + value.attrs["units"] + ")", size = "small")


    #lines to show location of transects
    #plt.plot([180, 180], [-180, 180],transform=ccrs.PlateCarree(), linestyle="--", color="black")
    #plt.plot([-30, -30], [-180, 180],transform=ccrs.PlateCarree(), linestyle="--", color="black")

    plt.title(title)
    plt.tight_layout()
    if save==True:
        fig.savefig(path + "plots/" + title +'.png', transparent=True, dpi=300)
    
def plot_depth_transect(data, cmap, lon, title):
    lat=data["latitude"]
    depth=data["depth"]
    
    x, y = np.meshgrid(lat, depth)
    
    fig, ax = plt.subplots()
    ax.set_xlabel("latitude")
    #ax.set_ylabel(data.attrs["long_name"] + "(" + data.attrs["units"] + ")")
    values = data[0,:,:,lon]
    ax.contourf(x,y,values, cmap=cmap)
    ax.invert_yaxis()
    ax.legend()
    
    c = ax.contourf(x, y, values, cmap=cmap)
    cbar = plt.colorbar(c)

    
    plt.title(title)
    plt.show()
    
    

# =============================================================================
# #%%
# 
# files = os.listdir()
# suffix=files[0].split("_")[-1]
# 
# for file in files:
#     if file.endswith('.nc'):
#         #import pdb; pdb.set_trace()
#         data = xr.open_dataset(file, decode_times=False)
#         lon = data["longitude"]
#         lat = data["latitude"]
#         if file == "L_ice_LIG.nc":
#             plot_polar_stereo(lon, lat, data["L_icethk"][0], "Blues", -10, 10, "L_icethk")
#             plot_polar_stereo(lon, lat, data["L_icefra"][0], "BuGn", 0, 1, "L_icefra")
#         else:
#             var=file.rstrip("_"+suffix)
#             if var.startswith("O") and not var.endswith("sur"):#only looking at surface of O vars
#                 value=data[var][0]
#             elif var.startswith("O") and var.endswith("sur"):
#                 value=data[var].mean("time") #surface values provided as monthly means; convert to annual for this purpose
#             else:
#                 value=data[var]
#             vmin=float(value.min())
#             vmax=float(value.max())
#             title=var
#             cmap="PiYG"
#             plot_polar_stereo(lon, lat, value,cmap,vmin,vmax, title, save=False, path=None)
#     
#     
# #%%
# #generated data (lig)
# g_kmt_lig = xr.open_dataset("G_kmt_LIG.nc", decode_times=False)
# g_mskt_lig = xr.open_dataset("G_mskt_LIG.nc", decode_times=False)
# g_mskhreg_lig = xr.open_dataset("G_mskhreg_LIG.nc", decode_times=False)
# 
# #uvic reference data (pd)
# 
# g_kmt_pd = xr.open_dataset("/Users/mackenfr/scripts/pism_to_uvic/inputs/uvic/G_kmt.nc", decode_times=False)
# g_mskt_pd = xr.open_dataset("/Users/mackenfr/scripts/pism_to_uvic/inputs/uvic/G_mskt.nc", decode_times=False)
# 
# lon = data["longitude"]
# lat = data["latitude"]
# 
# kmt_bin_lig=(g_kmt_lig.where(g_kmt_lig["G_kmt"]==0, 1))
# kmt_bin_pd=(g_kmt_pd.where(g_kmt_pd["G_kmt"]==0, 1))
# 
# mskt_minus_kmt_lig=g_mskt_lig["G_mskt"]-kmt_bin_lig
# mskt_minus_kmt_pd=g_mskt_pd["G_mskt"]-kmt_bin_pd
# 
# cmap="PiYG"
# vmin=-1
# vmax=1
# path="path"
# title="kmt"
# 
# plot_polar_stereo(lon, lat, mskt_minus_kmt_lig["G_kmt"], cmap, vmin, vmax, "G_mskt_LIG-G_kmt_LIG")
# plot_polar_stereo(lon, lat, mskt_minus_kmt_pd["G_kmt"], cmap, vmin, vmax, "G_mskt_PD-G_kmt_PD")
# 
# #%%
# 
# 
# plot_polar_stereo(lon, lat, g_kmt_pd["G_kmt"], cmap, 0, 19, "G_kmt_PD")
# plot_polar_stereo(lon, lat, g_mskt_pd["G_mskt"], cmap, 0, 1, "G_mskt_PD")
# 
# 
# kmt_pd_minus_lig=g_kmt_pd["G_kmt"]-g_kmt_lig["G_kmt"]
# mskt_pd_minus_lig=g_mskt_pd["G_mskt"]-g_mskt_lig["G_mskt"]
# 
# plot_polar_stereo(lon, lat, kmt_pd_minus_lig, cmap, 0, 19, "G_kmt_PD-G_kmt_LIG")
# 
# plot_polar_stereo(lon, lat, mskt_pd_minus_lig, cmap, 0, 1, "G_mskt_PD-G_mskt_LIG")
# 
# 
# #%%
# filename = "/Users/mackenfr/scripts/pism_to_uvic/outputs/LIG_RIS/G_mskt_LIG.nc"
# data = xr.open_dataset(filename, decode_times=False)
# lon = data["longitude"]
# lat = data["latitude"]
# value=data["G_mskt"]
# cmap="coolwarm"
# vmin=0
# vmax=4
# path="path"
# title="mskt"
# 
# plot_polar_stereo(lon, lat, value, cmap, vmin, vmax, path, title)
# 
# filename = "/Users/mackenfr/scripts/pism_to_uvic/outputs/LIG_RIS/L_elev_LIG.nc"
# data = xr.open_dataset(filename, decode_times=False)
# lon = data["longitude"]
# lat = data["latitude"]
# value=data["L_elev"]
# cmap="coolwarm"
# vmin=-1000
# vmax=4000
# path="path"
# title="L_elev"
# 
# plot_polar_stereo(lon, lat, value, cmap, vmin, vmax, path, title)
# 
# =============================================================================
