#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 08:30:49 2024

@author: mackenfr
"""
#%%
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
import pandas as pd
import seaborn as sns
os.chdir("/Users/mackenfr/scripts")
from calc_ovt import calc_ovt
#%% paths
path_to_data = "/Users/mackenfr/uvic_outputs/" 
config_names = ["PRD-PI", "AMN-PI", "TRS-PI", "EQB-PI"]
def import_data(path_to_data):
    config_names = ["PRD-PI", "AMN-PI", "TRS-PI", "EQB-PI"]# "PRD-LIG","AMN-LIG"
    co2s = ["280", "300", "350", "400", "450"]#
    tavgs = []
    names = []
    co2_all=[]
    ices = []
    orbs = []
    for config in config_names:
        for co2 in co2s:
            filename_tavg = path_to_data + "/" + config + "/" + co2 + "/tavg.nc"
            data = xr.open_dataset(filename_tavg, decode_times=False, engine = "netcdf4")
            ice, orb = config.split("-")
            names.append(config)
            tavgs.append(data)
            ices.append(ice)
            orbs.append(orb)
            co2_all.append(co2)
        
    ref_df = pd.DataFrame({"config": names,
                       "co2": co2_all,
                       "ice": ices,
                       "orb": orbs
                       })
                           
    return tavgs,ref_df

tavgs, ref_df = import_data(path_to_data)

picontrol = tavgs[0]
#%%
#names = ["PRD-PI", "AMN-PI"]
#["picontrol","pi_wais_lig","pi_wais_part" , "pi_wais_no"]

min_lon=0#80
max_lon =99#5

co2s = ["280", "300", "350", "400", "450"]#
for co2 in co2s:
    idx = ref_df[ref_df.co2==co2].index
    fig, ax = plt.subplots()
    for i in idx:
        if min_lon > max_lon:
            lon_range = list(range(0,max_lon+1))+list(range(min_lon, 100))
        else:
            lon_range = list(range(min_lon, max_lon+1))
        #zonal_mean_pi = tavgs[0]["F_heat"][0,:,lon_range].mean(dim="longitude")
        #import pdb; pdb.set_trace()
        
        zonal_mean = tavgs[i]["F_heat"][0,:,lon_range].mean(dim="longitude")
        zonal_mean_anom = zonal_mean #- zonal_mean_pi
        ax.plot(tavgs[i]["latitude"], zonal_mean_anom[:], label=ref_df.at[i, "config"])
    ax.set_ylim(-60,60)
    ax.legend()
    ax.set_title("Heat flux Global " +co2)
# %%
min_lon=0#80
max_lon =99#5
config_names = ["PRD-PI", "AMN-PI",  "TRS-PI", "EQB-PI",]
for config in config_names:
    idx = ref_df[ref_df.config==config].index
    fig, ax = plt.subplots()
    zonal_means = pd.DataFrame()
    for i in idx:
        co2 = ref_df.co2[i] 
        zonal_mean = tavgs[i]["F_heat"][0,:,lon_range].mean(dim="longitude")
        zonal_mean_anom = zonal_mean #- zonal_mean_pi
        zonal_means[co2] = zonal_mean_anom.values
    zonal_means["latitude"] = tavgs[0]["latitude"].values
    zonal_means = zonal_means.melt(id_vars="latitude", var_name="co2", value_name = "mean_flux")
    sns.lineplot(x="latitude", y="mean_flux", hue="co2", data = zonal_means, palette="plasma_r")
    ax.set_ylim(-60,60)
    ax.legend()
    ax.set_title("Heat flux Global " +config)
# %% by sector'
# Weddell: 60W - 0
# DML: 0-60E
# Indian: 60E - 160E
# Ross: 160E - 150W
# Amundsen: 150W - 60W
ocean_mask = picontrol["G_areaT"].where(picontrol["G_mskt"] == 1)
sectors = {"Weddell": list(range(82,100)),
           "DML": list(range(0,17)),
           "Indian": list(range(17,45)),
           "Ross": list(range(45,58)),
           "Amundsen":list(range(58,82))}
area = picontrol["G_areaT"]
for config in config_names:
    idx = ref_df[ref_df.config==config].index
    fig, ax = plt.subplots()
    sector_means = pd.DataFrame({"co2":[float(i) for i in co2s]})
    for sector in sectors:
        sector_mean = []
        for i in idx:
            ocean_area = area[:16, sectors[sector]].where(picontrol["G_mskt"][:16, sectors[sector]] == 1)
            max_area = ocean_area.max()
            area_weights = area[:16, sectors[sector]]/max_area
            heat_flux = (tavgs[i]["F_heat"][0,:16,sectors[sector]]*area_weights).mean()
            sector_mean.append(float(heat_flux))
        sector_means[sector] = sector_mean
    melted_means = sector_means.melt(id_vars="co2", var_name = "sector", value_name="heat flux")
    sns.lineplot(x="co2", y="heat flux", hue="sector", data=melted_means)
    sns.scatterplot(x="co2", y="heat flux", hue="sector", data=melted_means, legend=False)  
    ax.set_title(config +" <60S")
    ax.set_ylabel("downward heat flux")
    

# %% same but anom

area = picontrol["G_areaT"]
for config in config_names:
    idx = ref_df[ref_df.config==config].index
    fig, ax = plt.subplots()
    sector_means = pd.DataFrame({"co2":[float(i) for i in co2s]})
    for sector in sectors:
        sector_mean = []
        for i in idx:
            ocean_area = area[:24, sectors[sector]].where(picontrol["G_mskt"][:24, sectors[sector]] == 1)
            max_area = ocean_area.max()
            area_weights = area[:24, sectors[sector]]/max_area
            heat_flux_prd = (tavgs[i%5]["F_heat"][0,:24,sectors[sector]]*area_weights).mean()
            heat_flux = (tavgs[i]["F_heat"][0,:24,sectors[sector]]*area_weights).mean()
            anom = heat_flux - heat_flux_prd
            sector_mean.append(float(anom))
        sector_means[sector] = sector_mean
    melted_means = sector_means.melt(id_vars="co2", var_name = "sector", value_name="heat flux")
    sns.lineplot(x="co2", y="heat flux", hue="sector", data=melted_means)
    sns.scatterplot(x="co2", y="heat flux", hue="sector", data=melted_means, legend=False)  
    ax.set_title(config+" anomaly <45")
    ax.set_ylabel("downward heat flux")
        

# %%
basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}

area = picontrol["G_areaT"]
for config in config_names:
    idx = ref_df[ref_df.config==config].index
    fig, ax = plt.subplots()
    sector_means = pd.DataFrame({"co2":[float(i) for i in co2s]})
    for sector in basins:
        sector_mean = []
        for i in idx:
            j = basins[sector]
            ocean_area = area[:24, :].where(picontrol["G_mskhr"][:24,:] == j)
            max_area = ocean_area.max()
            area_weights = area[:24,:]/max_area
            heat_flux_prd = (tavgs[i%5]["F_heat"][0,:24,:].where(picontrol["G_mskhr"][:24,:] == j)*area_weights).mean()
            heat_flux = (tavgs[i]["F_heat"][0,:24,:].where(picontrol["G_mskhr"][:24,:] == j)*area_weights).mean()
            anom = heat_flux - heat_flux_prd
            sector_mean.append(float(anom))
        sector_means[sector] = sector_mean
    melted_means = sector_means.melt(id_vars="co2", var_name = "sector", value_name="heat flux")
    sns.lineplot(x="co2", y="heat flux", hue="sector", data=melted_means)
    sns.scatterplot(x="co2", y="heat flux", hue="sector", data=melted_means, legend=False)  
    ax.set_title(config+" anomaly <45")
    ax.set_ylabel("downward heat flux")
        

# %%
area = picontrol["G_areaT"]
for config in config_names:
    idx = ref_df[ref_df.config==config].index
    fig, ax = plt.subplots()
    sector_means = pd.DataFrame({"co2":[float(i) for i in co2s]})
    for sector in basins:
        sector_mean = []
        for i in idx:
            j = basins[sector]
            ocean_area = area[:24, :].where(picontrol["G_mskhr"][:24,:] == j)
            max_area = ocean_area.max()
            area_weights = area[:24,:]/max_area
            heat_flux_prd = (tavgs[i%5]["F_heat"][0,:24,:].where(picontrol["G_mskhr"][:24,:] == j)*area_weights).mean()
            heat_flux = (tavgs[i]["F_heat"][0,:24,:].where(picontrol["G_mskhr"][:24,:] == j)*area_weights).mean()
            anom = heat_flux# - heat_flux_prd
            sector_mean.append(float(anom))
        sector_means[sector] = sector_mean
    melted_means = sector_means.melt(id_vars="co2", var_name = "sector", value_name="heat flux")
    sns.lineplot(x="co2", y="heat flux", hue="sector", data=melted_means)
    sns.scatterplot(x="co2", y="heat flux", hue="sector", data=melted_means, legend=False)  
    ax.set_title(config+" abs <45")
    ax.set_ylabel("downward heat flux")
        
# %%
