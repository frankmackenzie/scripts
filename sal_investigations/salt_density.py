#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 13:30:22 2025

@author: mackenfr
"""
#%%
 import pandas as pd
import seaborn as sns
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import gsw
#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 

def import_data(path_to_data):
    config_names = ["PRD-PI", "AMN-PI"]
    co2s = ["280", "300", "350", "400", "450"]#
    tavgs = []
    names = []
    co2_all=[]
    ices = []
    orbs = []
    for config in config_names:
        for co2 in co2s:
            filename_tavg = path_to_data + "/" + config + "/" + co2 + "/tavg.nc"
            data = xr.open_dataset(filename_tavg, decode_times=False, engine = "scipy")
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
#%%
def calculate_density(data):
    lon = data["longitude"]
    lat = data["latitude"]
    SP = data["O_sal"] #practical salinity
    pt = data["O_temp"] #potential temp
    p = data["depth"] #sea pressure (approximate depth in m ~= pressure in dbar)
    SA = gsw.conversions.SA_from_SP(SP, p, lon, lat) #convert practical salinity into absolute salinity
    CT = gsw.conversions.CT_from_pt(SA, pt) #potential temperature to conservative temperature
    sigma0 = gsw.density.sigma0(SA, CT) #calculate potential density anomaly (anomaly from surface)
    return(sigma0)
#%%
lon_range = list(range(0,6))+list(range(80, 100))
max_area = tavgs[0]["G_areaT"][:, :].max()
weights = tavgs[0]["G_areaT"][:, :]/max_area

max_area_nh = tavgs[0]["G_areaT"][:, 90:].max()
weights_nh = tavgs[0]["G_areaT"][:, 90:]/max_area_nh

max_area_weddell = tavgs[0]["G_areaT"][:, :16].max()
weights_weddell  = tavgs[0]["G_areaT"][:, :16]/max_area_weddell


NH_salt = []
weddell_salt = []
for i in tavgs:
    nh = i["F_salt"][0, :, 90:]*weights_nh
    nh_mean= nh.mean(dim=["latitude", "longitude"])
    NH_salt.append(float(nh_mean.values))
    weddell = i["F_salt"][0, :, :16]*weights_weddell
    weddell_mean= weddell.mean(dim=["latitude", "longitude"])
    weddell_salt.append(float(weddell_mean.values))

# %%
salt_df = ref_df.copy()
salt_df["weddell"] = weddell_salt
salt_df["nh"] = NH_salt
plt.figure(figsize=(10, 4))
sns.scatterplot(x="co2", y="weddell", hue = "orb", style="ice", data = salt_df)
plt.title("SH")


plt.figure(figsize=(10, 4))
sns.scatterplot(x="co2", y="nh", hue = "orb", style="ice", data = salt_df)
plt.title("NH")


# %%
NA_density = []
NA_sal = []
NA_temp = []
weddell_density = []
weddell_sal = []
weddell_temp = []
for i in tavgs:
    density = calculate_density(i)
    weddell_dens = density[0,:,10:14,86:90].mean(dim=["longitude", "latitude"])
    weddell_osal = i["O_sal"][0,:,10:14,86:90].mean(dim=["longitude", "latitude"])
    weddell_otemp = i["O_temp"][0,:,10:14,86:90].mean(dim=["longitude", "latitude"])
    weddell_sal.append(weddell_osal.values)
    weddell_temp.append(weddell_otemp.values)
    weddell_density.append(weddell_dens.values)

    NA_dens = density[0,:,79:83,91:94].mean(dim=["longitude", "latitude"])
    NA_osal = i["O_sal"][0,:,79:83,91:94].mean(dim=["longitude", "latitude"])
    NA_otemp = i["O_temp"][0,:,79:83,91:94].mean(dim=["longitude", "latitude"])
    NA_sal.append(NA_osal.values)
    NA_temp.append(NA_otemp.values)
    NA_density.append(NA_dens.values)
# %%
ref_df["full_name"]=ref_df["config"]+ref_df["co2"]
#idxs = ref_df[ref_df.co2=="280"].index
columns = [ref_df["full_name"][i] for i in ref_df.index]
dens_values = [weddell_density[i] for i in ref_df.index]
sal_values = [weddell_sal[i] for i in ref_df.index]
temp_values = [weddell_temp[i] for i in ref_df.index]

ws_density_df = pd.DataFrame({k:v for (k,v) in zip(columns, dens_values)})
ws_temp_df = pd.DataFrame({k:v for (k,v) in zip(columns, temp_values)})
ws_sal_df = pd.DataFrame({k:v for (k,v) in zip(columns, sal_values)})

ws_density_df["depth"] = tavgs[0]["depth"]
ws_temp_df["depth"] = tavgs[0]["depth"]
ws_sal_df["depth"] = tavgs[0]["depth"]

# %%
fig, axs = plt.subplots(2, 5, figsize=(12, 9), sharey = True, sharex =True)
axs[0,0].invert_yaxis()
for ax, config in zip(axs.flat, ws_density_df):
    ax1 = ax.twiny()
    ax2 = ax.twiny()

    plot = sns.lineplot(x=config, y="depth",  data = ws_density_df, ax=ax,  orient='y')
    plot1 = sns.lineplot(x=config, y="depth",  data = ws_sal_df, ax=ax1, color="green", orient='y')
    plot2 = sns.lineplot(x=config, y="depth",  data = ws_temp_df, ax=ax2,  color="orange", orient='y')

    plot.set(xlabel='', ylabel='depth')
    plot1.set(xlabel='', ylabel='')
    plot2.set(xlabel='', ylabel='')

for n, col in enumerate(ref_df["co2"].unique()):
    axs[0, n].set_title(col, fontsize=15)
for n, row in enumerate(ref_df["config"].unique()):
    fig.text(0.02, 0.71 - (n * 0.42),row, va='center', ha='right', fontsize=12, rotation=90)
fig.suptitle("Weddell sea", y=0.99, size=20)
fig.supxlabel("density", y=0.05)

# %%
columns = [ref_df["full_name"][i] for i in ref_df.index]
dens_values = [NA_density[i] for i in ref_df.index]
sal_values = [NA_sal[i] for i in ref_df.index]
temp_values = [NA_temp[i] for i in ref_df.index]

na_density_df = pd.DataFrame({k:v for (k,v) in zip(columns, dens_values)})
na_temp_df = pd.DataFrame({k:v for (k,v) in zip(columns, temp_values)})
na_sal_df = pd.DataFrame({k:v for (k,v) in zip(columns, sal_values)})

na_density_df["depth"] = tavgs[0]["depth"]
na_temp_df["depth"] = tavgs[0]["depth"]
na_sal_df["depth"] = tavgs[0]["depth"]
# %% na water definitions from liu and tanhua - these are neutral densities not potential sigma0
#saiw = [27.7, 27.85] #>27.7
#nadw = [27.85, 28.1]
#ventilated pycnocline - lumpkin and speer 2003 - <27.5 

#%%
fig, axs = plt.subplots(2, 5, figsize=(12, 9), sharey = True, sharex =True)
axs[0,0].invert_yaxis()
for ax, config in zip(axs.flat, ws_density_df):
    ax1 = ax.twiny()
    ax2 = ax.twiny()

    plot = sns.lineplot(x=config, y="depth",  data = na_density_df, ax=ax,  orient='y')
    plot1 = sns.lineplot(x=config, y="depth",  data = na_sal_df, ax=ax1, color="green", orient='y')
    plot2 = sns.lineplot(x=config, y="depth",  data = na_temp_df, ax=ax2,  color="orange", orient='y')

    plot.set(xlabel='', ylabel='depth')
    plot1.set(xlabel='', ylabel='')
    plot2.set(xlabel='', ylabel='')

for n, col in enumerate(ref_df["co2"].unique()):
    axs[0, n].set_title(col, fontsize=15)
for n, row in enumerate(ref_df["config"].unique()):
    fig.text(0.02, 0.71 - (n * 0.42),row, va='center', ha='right', fontsize=12, rotation=90)
fig.suptitle("North Atlantic", y=0.99, size=20)
fig.supxlabel("density", y=0.05)
#%% simpler vers dens only
fig, axs = plt.subplots(2, 5, figsize=(6, 6), sharey = True, sharex =True)
axs[0,0].invert_yaxis()
for ax, config in zip(axs.flat, na_density_df):
    #ax.axvspan(xmin=saiw[0], xmax=saiw[1], facecolor="green", alpha=0.3)
    #ax.axvspan(xmin=nadw[0], xmax=nadw[1], facecolor="blue", alpha=0.3)
    ax.axvspan(xmin=27, xmax=27.5, facecolor="red", alpha=0.3)
    hi = sns.scatterplot(x=config, y="depth",  data = na_density_df, ax=ax)
    hi.set(xlabel='', ylabel='depth')



for n, col in enumerate(ref_df["co2"].unique()):
    axs[0, n].set_title(col, fontsize=15)
for n, row in enumerate(ref_df["config"].unique()):
    fig.text(0.02, 0.71 - (n * 0.42),row, va='center', ha='right', fontsize=12, rotation=90)
fig.suptitle("North Atlantic", y=0.96, size=20)
fig.supxlabel("density", y=0.05)

# %% eyeballing the pycnocline depth and density in weddell
prd_depth = [27.7, 27.7,27.7,27.6,27.5]
prd_dens = [3,3, 3, 3, 3]
amn_depth = [27.6, 27.6, 27.6, 27.6, 27.5]
amn_dens = [2, 3, 3, 3, 3]
#ok so from this it seems like the pycnocline tends to be around 27.5-27.7, 
# but I wouldn;t say the depth is changing in the wedell . Only thing I'd say 
# is that the surface starts denser in amn and gets less dense with co2 (SI loss)

#in NA - much shallower pycnocline: in the first depth level, then int water (?)
#perists to c. 1000m
prd_depth = []
prd_dens = []
amn_depth = []
amn_dens = []
#%% mean density and density gradient   
na_density_diff = pd.DataFrame(columns = na_density_df.columns)
for col in na_density_diff.columns:
    na_density_diff["shifted"] = na_density_df[col].shift()
    na_density_diff[col] = na_density_diff["shifted"]-na_density_df[col]
na_density_diff["depth"] = na_density_df["depth"]

#%%
fig, axs = plt.subplots(2, 5, figsize=(6, 6), sharey = True, sharex =True)
axs[0,0].invert_yaxis()
for ax, config in zip(axs.flat, na_density_diff):
    hi = sns.scatterplot(x=config, y="depth",  data = na_density_diff, ax=ax)
    hi.set(xlabel='', ylabel='depth')



for n, col in enumerate(ref_df["co2"].unique()):
    axs[0, n].set_title(col, fontsize=15)
for n, row in enumerate(ref_df["config"].unique()):
    fig.text(0.02, 0.71 - (n * 0.42),row, va='center', ha='right', fontsize=12, rotation=90)
fig.suptitle("na", y=0.96, size=20)
fig.supxlabel("chnage in density", y=0.05)

# %%
