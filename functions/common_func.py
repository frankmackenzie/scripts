
import xarray as xr
import os
import numpy as np
import gsw
os.chdir("/Users/mackenfr/scripts")

import pandas as pd

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

def calculate_mean_density(data, basin):
    sigma0 = calculate_density(data)
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    G_mskhr = data["G_mskhr"]
    if basin=="Global":
        mean = sigma0.mean(dim=["longitude"])
    elif basin in basins:
        mask = sigma0.where(G_mskhr==basins[basin])
        mean = mask[:,:,:,:].mean(dim=["longitude"])
    return(mean)

def calc_mean(data, var, basin):
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    G_mskhr = data["G_mskhr"]
    if basin=="Global":
        mean = data[var][:,:,:,:].mean(dim=["longitude"])
    elif basin in basins:
        mask = data[var].where(G_mskhr==basins[basin])
        mean = mask[:,:,:,:].mean(dim=["longitude"])
    return(mean)

def get_2d_slice(da, lev=0):
    if "time" in da.dims:
        da = da.isel(time=0)
    if "depth" in da.dims:
        da = da.isel(depth=lev)
    if "depth_W" in da.dims:
        da = da.isel(depth_W=lev)
    return da