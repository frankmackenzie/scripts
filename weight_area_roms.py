#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d



#%%

roms_path = "/Users/mackenfr/Downloads/ROMS/"
pd_roms = xr.open_dataset(roms_path + "023_grd_ACOSM9km24_SeaWAIS_PD_e02.nc")
#%%
lat_roms= pd_roms["lat_rho"]
lon_roms= pd_roms["lon_rho"]
h = pd_roms["h"]
landmask = pd_roms["mask_rho"]
pm = pd_roms["pm"]
pn = pd_roms["pn"]
# %%


masked_cont = xr.where(((h<1000) & (landmask==1)), 1, 0)

fig = plt.figure()
ax = plt.axes(projection=ccrs.SouthPolarStereo())
# Limit the map to -45 degrees latitude and below.
ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(lon_roms, lat_roms, masked_cont,transform=ccrs.PlateCarree(),cmap="Spectral", shading="auto")

#colorbar 
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)

#%%
area_roms = (1/pm)*(1/pn)
shelf_area = masked_cont*area_roms
# %%
uvic_path = "/Users/mackenfr/uvic_outputs/PRD-PI/280/tavg.nc"
pd_uvic = xr.open_dataset(uvic_path, decode_times=False)
lons_uvic = pd_uvic["longitude"]
lats_uvic = pd_uvic["latitude"]

area_uvic = pd_uvic["G_areaT"]
kmt_uvic = pd_uvic["G_kmt"]
lon_edges_uvic = pd_uvic["longitude_edges"]
lat_edges_uvic = pd_uvic["latitude_edges"]

#match roms domain

# %
#%%
lon_roms_flat = lon_roms.values.ravel()
lat_roms_flat = lat_roms.values.ravel()
shelf_area_flat = shelf_area.values.ravel()
#%%

test_area_binned, xx, yy, _ = binned_statistic_2d(
    lon_roms_flat, lat_roms_flat, area_roms.values.ravel(),
    statistic='sum', 
    bins=[lon_edges_uvic, lat_edges_uvic]
)
shelf_area_binned, xx, yy, _ = binned_statistic_2d(
    lon_roms_flat, lat_roms_flat, shelf_area_flat,
    statistic='sum', 
    bins=[lon_edges_uvic, lat_edges_uvic]
)


# %%
fig = plt.figure()
ax = plt.axes(projection=ccrs.SouthPolarStereo())
# Limit the map to -45 degrees latitude and below.
ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(xx, yy, test_area_binned.T,transform=ccrs.PlateCarree(),cmap="YlGnBu", shading="auto")
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)

fig = plt.figure()
ax = plt.axes(projection=ccrs.SouthPolarStereo())
# Limit the map to -45 degrees latitude and below.
ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(xx, yy, shelf_area_binned.T,transform=ccrs.PlateCarree(),cmap="YlGnBu", shading="auto")
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)


# %%
#fraction of cell that is shelf
#first mask out land
masked_cont = xr.where(((h<1000) & (landmask==1)), 1, 0)

shelf_area_uvic = xr.where(kmt_uvic>0, shelf_area_binned.T, 0)

fig = plt.figure()
ax = plt.axes(projection=ccrs.SouthPolarStereo())
# Limit the map to -45 degrees latitude and below.
ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(xx, yy, shelf_area_uvic,transform=ccrs.PlateCarree(),cmap="YlGnBu", shading="auto")
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)

# %%
shelf_area_fraction = shelf_area_uvic/area_uvic
# %%

fig = plt.figure()
ax = plt.axes(projection=ccrs.SouthPolarStereo())
# Limit the map to -45 degrees latitude and below.
ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(xx, yy, shelf_area_fraction,transform=ccrs.PlateCarree(),cmap="YlGnBu", shading="auto")
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)

# %%
