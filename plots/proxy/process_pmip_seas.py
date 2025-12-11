#%%
import xarray as xr
import cftime
import datetime
import os
import numpy as np

models = ['CESM2',
 'NorESM1',
 'FGOALS-f3',
 'FGOALS-g3',
 'NESM3',
 'NorESM2',
 'HadGEM3',
 'ACCESS',
 'GISS']

pi_path = "/Users/mackenfr/proxy_compare_2025/pmip_data/pcal_picontrol/"
lig_path = "/Users/mackenfr/proxy_compare_2025/pmip_data/pcal_127k/"
save_path = "/Users/mackenfr/proxy_compare_2025/pmip_data/anoms_processed/"
seas = []
ann=[]

for model in models:
    data_pi = xr.open_dataset(pi_path + model + ".nc", chunks={"time":100})
    data_lig = xr.open_dataset(lig_path + model + ".nc",chunks={"time":100})

    new_monlen_pi = data_pi.rmonlen
    new_monlen_lig = data_lig.rmonlen

    weights_pi = (new_monlen_pi.groupby("time.season") / new_monlen_pi.groupby("time.season").sum())
    weights_lig = (new_monlen_lig.groupby("time.season") / new_monlen_lig.groupby("time.season").sum())

    seas_pi = (data_pi["tos"]*weights_pi).groupby("time.season").sum(dim="time")
    seas_lig = (data_lig["tos"]*weights_lig).groupby("time.season").sum(dim="time")

    summer_anom = ((seas_lig[0,:,:] - seas_pi[0,:,:])
                   .where(~np.isnan(data_pi["tos"][0,:,:]))
                   .to_dataset(name="tos"))

    summer_anom.to_netcdf(save_path + model + "summer_anom.nc")

    weights_ann_pi = (new_monlen_pi / new_monlen_pi.sum())
    pi_ann = (data_pi["tos"]*weights_ann_pi).sum(dim="time")

    weights_ann_lig = (new_monlen_lig / new_monlen_lig.sum())
    lig_ann = (data_lig["tos"]*weights_ann_lig).sum(dim="time")

    ann_anom = ((lig_ann - pi_ann)
                .where(~np.isnan(data_pi["tos"][0,:,:])
                .to_dataset(name="tos"))
                )
    
    ann_anom.to_netcdf(save_path + model + "ann_anom.nc")
    ann.append(ann)
    seas.append(summer_anom)
# %%

# %%