
#%%
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import colormaps
import seaborn as sns
import xarray as xr
import cartopy.crs as ccrs
import os
import numpy as np
import gsw
os.chdir("/Users/mackenfr/scripts")

import pandas as pd

from import_data import import_data
from calc_ovt import calc_ovt
from common_func import calculate_mean_density, calculate_density, calc_mean, get_2d_slice

#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
config_names = ["PRD-PI", "PRD-LIG", "AMN-PI","AMN-LIG", "EQB-PI",  "EQB-LIG"] 
co2 = ["280"]
tavgs, ref_df = import_data(path_to_data, co2, config_names)
#%%
subset = ref_df.loc[ref_df["orb"]=="LIG"].reset_index()
control_subset = ref_df.loc[ref_df["orb"]=="PI"].reset_index()
datasets = [tavgs[i] for i in subset["index"]]
control_datasets = [tavgs[i] for i in control_subset["index"]] 
#%%
anom_df = pd.DataFrame({"Latitude": tavgs[0]["latitude"].values})
lig_df = pd.DataFrame({"Latitude": tavgs[0]["latitude"].values})
pi_df = pd.DataFrame({"Latitude": tavgs[0]["latitude"].values})

for i in list(range(3)):
    lig_ind = subset["index"][i]
    lig = tavgs[lig_ind]
    pi_ind = control_subset["index"][i]
    pi = tavgs[pi_ind]
    zmean_lig = lig["A_sat"].mean(dim="longitude")
    zmean_pi = pi["A_sat"].mean(dim="longitude")
    lig_df[subset["ice"][i]] = zmean_lig[0,:].values
    pi_df[subset["ice"][i]] = zmean_pi[0,:].values
    anom = zmean_lig[0,:] - zmean_pi[0,:]
    anom_df[subset["ice"][i]] = anom.values


sat_pi = pi_df.melt(id_vars= "Latitude", var_name="Ice sheet", value_name="SAT")
sat_lig = lig_df.melt(id_vars= "Latitude", var_name="Ice sheet", value_name="SAT")
sat_anom = anom_df.melt(id_vars= "Latitude", var_name="Ice sheet", value_name="SAT anom")

fig,ax = plt.subplots()
sns.lineplot(x="Latitude", y="SAT", data = sat_pi, hue = "Ice sheet", legend=True)

fig,ax = plt.subplots()
sns.lineplot(x="Latitude", y="SAT", data = sat_lig, hue = "Ice sheet", legend=True)

fig,ax = plt.subplots()
sns.lineplot(x="Latitude", y="SAT anom", data = sat_anom, hue = "Ice sheet", legend=True)


# %%
prdpi_anom_df = pd.DataFrame({"Latitude": tavgs[0]["latitude"].values})

for i in list(range(3)):
    lig_ind = subset["index"][i]
    lig = tavgs[lig_ind]
    pi = tavgs[0]
    zmean_lig = lig["A_sat"].mean(dim="longitude")
    zmean_pi = pi["A_sat"].mean(dim="longitude")
    lig_df[subset["ice"][i]] = zmean_lig[0,:].values
    pi_df[subset["ice"][i]] = zmean_pi[0,:].values
    anom = zmean_lig[0,:] - zmean_pi[0,:]
    prdpi_anom_df[subset["config"][i]] = anom.values
prdpi_sat_anom = anom_df.melt(id_vars= "Latitude", var_name="Ice sheet", value_name="SAT anom")

fig,ax = plt.subplots()
sns.lineplot(x="Latitude", y="SAT anom", data = prdpi_sat_anom, hue = "Ice sheet", legend=True)
ax.set_title("from PRDPI")
ax.set_xlim(-50, 90)
ax.set_ylim(-1, 1)

# %%
