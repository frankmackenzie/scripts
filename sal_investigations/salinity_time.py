#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 11:37:06 2023

@author: mackenfr
"""
import matplotlib.pyplot as plt
import xarray as xr
import cftime
import numpy as np

lig_path = "/Users/mackenfr/data/uvic_outputs/lig_wais_lig/"
picontrol_path = "/Users/mackenfr/data/uvic_outputs/picontrol/"

lig_tsi_file =lig_path+ "tsi_lig_wais_lig.nc"
#pi_tsi = get this from gns

lig_tavg_file = lig_path+"tavg.10000.01.01.nc"


lig_tavg= xr.open_dataset(lig_tavg_file, decode_times=False)

lig_tsi=xr.open_dataset(lig_tsi_file, decode_times=False)
