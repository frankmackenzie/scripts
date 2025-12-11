#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:19:04 2023

@author: mackenfr
"""
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import os

os.chdir("/Users/mackenfr/scripts")

from plotting_funcs import plot_polar_stereo, plot_global 

#%% paths
path_to_data = "/Users/mackenfr/uvic_outputs/" 

savepath = "/Users/mackenfr/Desktop/"
#%%
lig_file = path_to_data + "lig_wais_lig/280/tavg.10000.01.01.nc"
lig = xr.open_dataset(lig_file, decode_times=False)
#if more than one timestep

pi_file = path_to_data + "pi_wais_lig/280/tavg.10000.01.01.nc"
pi_wais_lig = xr.open_dataset(pi_file, decode_times=False)

picontrol_file = path_to_data + "picontrol/280/tavg.10000.01.01.nc"#"/Users/mackenfr/data/uvic_outputs/pi_wais_lig/pi_wais_lig_tavg.09000.01.01.nc"
picontrol = xr.open_dataset(picontrol_file, decode_times=False)

#if average from dif no of timesteps
lig, picontrol = xr.align(lig, picontrol, join='override')

lon = lig["longitude"]
lat = lig["latitude"]

#%% sst
vmin = -1
vmax = 1
 
var = "O_temp"
anom_lig = lig[var]-picontrol[var]
anom_lig.attrs = lig[var].attrs

plot_global(lon, lat, anom_lig[0, 0, :, :], "coolwarm",vmin,vmax,"SST lig_wais_lig - picontrol", save=True, path=savepath)
plot_polar_stereo(lon, lat, anom_lig[0, 0, :, :], "coolwarm",vmin,vmax, "SST lig_wais_lig - picontrol", save=True, path=savepath)
 

anom_pi = pi_wais_lig[var]-picontrol[var]
anom_pi.attrs = pi_wais_lig[var].attrs

plot_global(lon, lat, anom_pi[0, 0, :, :], "coolwarm",vmin,vmax, "SST pi_wais_lig - picontrol")
plot_polar_stereo(lon, lat, anom_pi[0, 0, :, :], "coolwarm",vmin,vmax, "SST pi_wais_lig - picontrol")

vmin = -10
vmax = 30
plot_global(lon, lat, lig[var][0, 0, :, :], "coolwarm",vmin,vmax,"SST lig_wais_lig", save=True, path=savepath)
plot_polar_stereo(lon, lat, lig[var][0, 0, :, :], "coolwarm",vmin,vmax, "SST lig_wais_lig", save=True, path=savepath)

plot_global(lon, lat, pi_wais_lig[var][0, 0, :, :], "coolwarm",vmin,vmax,"SST pi_wais_lig", save=True, path=savepath)
plot_polar_stereo(lon, lat, pi_wais_lig[var][0, 0, :, :], "coolwarm",vmin,vmax, "SST pi_wais_lig", save=True, path=savepath)

#%% salinity
vmin=-1.5
vmax=1.5
var = "O_sal"
anom_lig = lig[var]-picontrol[var]
anom_lig.attrs = lig[var].attrs

plot_global(lon, lat, anom_lig[0, 0, :, :], "PiYG",vmin,vmax, "sal_global_anom", save=True, path=savepath)#"pi_wais_lig - picontrol"
plot_polar_stereo(lon, lat, anom_lig[0, 0, :, :], "PiYG",vmin,vmax, "sal_polar_anom",save=True, path=savepath)#"pi_wais_lig - picontrol"

vmin=30
vmax=40

plot_global(lon, lat, lig[var][0, 0, :, :], "PiYG",vmin,vmax, "sal_global_lig", save=True, path=savepath)#"pi_wais_lig - picontrol"
plot_polar_stereo(lon, lat, lig[var][0, 0, :, :], "PiYG",vmin,vmax, "sal_polar_lig",save=True, path=savepath)#"pi_wais_lig - picontrol"

# anom_pi = pi_wais_lig[var]-picontrol[var]
# anom_pi.attrs = pi_wais_lig[var].attrs

#%% salinity
#for i in list(range(19)):
   # plot_polar_stereo(lon, lat, lig[var][0, i, :, :], "PiYG",30,40, path_to_data,"lig_absolute_sal_polar")
#plot_polar_stereo(lon, lat, anom_pi[0, 0, :, :], "PiYG",vmin,vmax, "pi_wais_lig - picontrol")
#%% temperature

var = "A_sat"
vmin=-10
vmax=10
anom_lig = lig[var]-picontrol[var]
anom_lig.attrs = lig[var].attrs
anom_pi = pi_wais_lig[var]-picontrol[var]
anom_pi.attrs = pi_wais_lig[var].attrs

plot_global(lon, lat, anom_lig[0, :, :], "RdYlBu_r",-2,2, "SAT lig_wais_lig - picontrol", save=True, path=savepath)#"pi_wais_lig - picontrol"
plot_polar_stereo(lon, lat, anom_lig[0, :, :], "RdYlBu_r",vmin,vmax,"SAT lig_wais_lig - picontrol", save=True, path=savepath)#"pi_wais_lig - picontrol"

plot_global(lon, lat, anom_pi[0, :, :], "RdYlBu_r",-2,2,"SAT pi_wais_lig - picontrol", save=True, path=savepath)#"pi_wais_pi - picontrol"
plot_polar_stereo(lon, lat, anom_pi[0, :, :], "RdYlBu_r",vmin,vmax,"SAT pi_wais_lig - picontrol", save=True, path=savepath)#"pi_wais_pi - picontrol"
#%%
vmin=-30
vmax=30
plot_global(lon, lat, lig[var][0, :, :], "RdYlBu_r",vmin,vmax,"SAT lig_wais_lig", save=True, path=savepath)
plot_polar_stereo(lon, lat, lig[var][0, :, :], "RdYlBu_r",vmin,vmax, "SAT lig_wais_lig", save=True, path=savepath)

plot_global(lon, lat, pi_wais_lig[var][0, :, :], "RdYlBu_r",vmin,vmax,"SAT pi_wais_lig", save=True, path=savepath)
plot_polar_stereo(lon, lat, pi_wais_lig[var][0, :, :], "RdYlBu_r",vmin,vmax, "SAT pi_wais_lig", save=True, path=savepath)

#plot_polar_stereo(lon, lat, anom_pi[0, :, :], "RdYlBu_r",vmin,vmax, "pi_wais_lig - picontrol")#%% temperature

#%% precipitation
vmin=0
vmax=0.00005
var = "F_precip"
anom_lig = lig[var]-picontrol[var]
anom_lig.attrs = lig[var].attrs

plot_global(lon, lat, lig[var][0, :, :], "BrBG",vmin,vmax, path_to_data,"precip_global")#"pi_wais_lig - picontrol"
plot_global(lon, lat, picontrol[var][0, :, :], "BrBG",vmin,vmax, path_to_data,"precip_global")#"pi_wais_lig - picontrol"

#plot_polar_stereo(lon, lat, anom_lig[0, :, :], "BrBG",vmin,vmax,path_to_data,"precip_global")#"pi_wais_lig - picontrol"

# anom_pi = pi_wais_lig[var]-picontrol[var]
# anom_pi.attrs = pi_wais_lig[var].attrs

# plot_global(lon, lat, anom_pi[0, :, :], "RdYlBu_r",vmin,vmax, "pi_wais_lig - picontrol")
# plot_polar_stereo(lon, lat, anom_pi[0, :, :], "RdYlBu_r",vmin,vmax, "pi_wais_lig - picontrol")



#%% sea ice edge

var = "O_icefra"
for i in [lig, pi_wais_lig, picontrol]:
    ice_mask = i[var][0, :, :].where(i["O_icefra"] > 0.15)
    ice_edge = ice_mask.notnull() & ice_mask.shift(latitude = -1).isnull()
    plot_polar_stereo(lon, lat, ice_edge[:, :, 0],"BuGn" ,vmin,vmax, i.attrs["run_stamp"][94:])#"pi_wais_lig - picontrol"


#%%
vmin=0
vmax=1
plot_polar_stereo(lon, lat, lig[var][0, :, :], "PRGn",vmin,vmax,"siconc_lig_polar")#"pi_wais_lig - picontrol"
plot_polar_stereo(lon, lat, pi_wais_lig[var][0, :, :], "PRGn",vmin,vmax,"siconc_pi_polar")
plot_polar_stereo(lon, lat, picontrol[var][0, :, :], "PRGn",vmin,vmax, "siconc_picontrol_polar")#"pi_wais_lig - picontrol"


#%% river discharge
var = "F_rivdis"
anom = lig[var]-picontrol[var]
anom.attrs = lig[var].attrs

plot_polar_stereo(lon, lat, lig["F_rivdis"][0, :, :], "blues",0,0.0001, path_to_data,"siconc_anom_polar")#"pi_wais_lig - picontrol"
plot_polar_stereo(lon, lat, lig["F_rivdis"][0, :, :], "blues",0,0.0001, path_to_data,"siconc_anom_polar")#"pi_wais_lig - picontrol"

#plot_global(lon, lat, lig["F_rivdis"][0, :, :], "blues",0,0.0005, path_to_data,"siconc_anom_polar")#"pi_wais_lig - picontrol"

fig = plt.figure()

ax = plt.axes(projection=ccrs.NorthPolarStereo())
# Limit the map to -60 degrees latitude and below.
ax.set_extent([-180, 180, 90, 40], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(lon, lat, lig["F_rivdis"][0, :, :],transform=ccrs.PlateCarree(),vmin=0,vmax=0.0001,cmap="blues", shading="auto")

#colorbar 
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)


#plt.title(title)
plt.tight_layout()


fig = plt.figure()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
#ax.patch.set_alpha(0)   
ax.set_facecolor('white')
 #ax.gridlines(linestyle='--')

ct=ax.pcolormesh(lon,lat,lig["F_rivdis"][0, :, :],transform=ccrs.PlateCarree(),vmin=0,vmax=0.0001,cmap="blues")#, shading="auto")

 #colorbar 
cb = plt.colorbar(ct,orientation="vertical",fraction=0.023, pad=0.04)





#%% icefra arctic
anom = lig["O_icefra"] - picontrol["O_icefra"]

fig = plt.figure()

ax = plt.axes(projection=ccrs.NorthPolarStereo())
# Limit the map to -60 degrees latitude and below.
ax.set_extent([-180, 180, 90, 50], ccrs.PlateCarree())
ax.coastlines() 
ax.gridlines(linestyle='--', draw_labels=True, color="gray")

ct=ax.pcolormesh(lon, lat, anom[0, :, :],transform=ccrs.PlateCarree(),vmin=-1,vmax=1,cmap="PRGn", shading="auto")

#colorbar 
cb = plt.colorbar(ct,orientation="vertical",fraction=0.046, pad=0.04)



#%% checking masks

# =============================================================================
# 
# vmax=2
# vmin=0
# plot_polar_stereo(lon, lat, lig["G_mskt"], "GnBu",vmin,vmax, "msk_lig")
# plot_polar_stereo(lon, lat, picontrol["G_mskt"], "GnBu",vmin,vmax, "msk_pi")
# ocean_lvls = np.where(lig["G_kmt"].values >= 3, 1, 0)
# lig["G_kmt"].values = ocean_lvls
# 
# ocean_lvls = np.where(picontrol["G_kmt"].values >= 3, 1, 0)
# picontrol["G_kmt"].values = ocean_lvls
# 
# vmax=1
# vmin=0
# var = "G_kmt"
# plot_polar_stereo(lon, lat, lig[var][:, :], "GnBu",vmin,vmax, path_to_data,"msk_lig")#"pi_wais_lig - picontrol"
# plot_polar_stereo(lon, lat, picontrol[var][:, :], "GnBu",vmin,vmax, path_to_data,"msk_pi")#"pi_wais_lig - picontrol"
# 
# n=1063
# =============================================================================
#plot_polar_stereo(lon, lat, data["O_icefra"][n,:, :], "Blues_r",vmin,vmax, path_to_data,"polynya")#"pi_wais_lig - picontrol"
