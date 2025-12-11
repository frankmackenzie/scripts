#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 15:51:01 2024

@author: mackenfr
"""
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import os
import numpy as np
os.chdir("/Users/mackenfr/scripts")

from calc_ovt import calc_ovt

path_to_data = "/Users/mackenfr/uvic_outputs/" 

configs = ["picontrol"] #, "lig_wais_pi", "pi_wais_lig",  "lig_wais_lig"]
co2s = ["280_wind",  "300_wind",  "350_wind", "400_wind", "450_wind"]
co2ccn = ["280", "300", "350", "400", "450"]
tsis = []
names = []
co2_all = []

for config in configs:
    for co2 in co2s:
        filename_tsi = path_to_data + "/" + config + "/" + co2 + "/tsi_all.nc"
        data = xr.open_dataset(filename_tsi, decode_times=False, engine = "scipy")
        tsis.append(data)
        names.append(config)
        co2_all.append(co2.strip("_wind"))
#%%
        
O_temp = []
for i in tsis:
    otempsur = i["O_tempsur"][:100].mean()
    O_temp.append(float(otempsur))
    
table = pd.DataFrame({"Experiment": names,"CO2": pd.to_numeric(co2_all), "Mean SST": O_temp})

plt.figure(figsize=(10, 6))

# Filter out NaN values
filtered_df = table[~table["Mean SST"].isna()]
colors = cm.get_cmap('viridis', len(filtered_df["Experiment"].unique()))

for i, experiment in enumerate(filtered_df["Experiment"].unique()):
     subset = filtered_df[filtered_df["Experiment"] == experiment]
     subset.sort_values(by=["CO2"], inplace=True)
     plt.scatter(subset["CO2"], subset["Mean SST"], label=experiment)

plt.xlabel('CO$_2$ (ppm)')

plt.ylabel("global mean SST ($^\circ$C)")
plt.legend()
#plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(scientific_format))
