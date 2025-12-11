#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 12:14:22 2023

load two nc files, extract the G_MSKRG from one dataset and assign it to another, save the new dataset as an nc file

@author: mackenfr
"""

import xarray as xr

file1 = "/Users/mackenfr/UVic_ESCM/outputs/pi_wais_lig/tavg_100y_avg_timereset.nc"
file2 = "/Users/mackenfr/UVic_ESCM/outputs/pi_wais_lig/tavg_pi_wais_lig_5day.nc"
data1 = xr.open_dataset(file1, decode_cf=True, decode_times=False)
data2 = xr.open_dataset(file2, decode_cf=True, decode_times=False)

data1["G_mskhr"]=data2["G_mskhr"]

data1.to_netcdf("tavg_100y_avg_timereset_adj.nc")