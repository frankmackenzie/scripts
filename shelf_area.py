#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 13:57:59 2024

@author: mackenfr
"""
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os
os.chdir("/Users/mackenfr/scripts")
from plotting_funcs import plot_polar_stereo
path_to_data = "/Users/mackenfr/uvic_outputs/" 

def shelf_area(mask):
    path = path_to_data + mask + "/280/tavg.10000.01.01.nc"
    data = xr.open_dataset(path, decode_times=False)
    SO = data["G_kmt"][:20, :]
    shelf_mask = SO.where((SO>0) & (SO<7))#closest cell to 1000m has bottom edge at 980m
    area_mask = data["G_areaT"].where(np.isnan(shelf_mask) == False)
    area_sum = area_mask.sum()
    return(area_sum)
    
area_pi = float(shelf_area("picontrol"))
area_lig = float(shelf_area("pi_wais_lig"))

print(area_pi)
print(area_lig)
print(area_lig/area_pi)