#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 14:52:19 2023

takes the outputs from Nick Golledge et al.'s PISM runs from their 2021 paper, and maps 
them to the relevant input files to be used with UVic ESCM. Due to some issues 
that arise from mapping data from a fine to a coarse mesh,  a few manual 
adjustments need to be made for certain masks (lig_wais in particular) to be 
'accepted' by the UVic program.


@author: mackenfr
"""
from netCDF4 import Dataset, MFDataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import xarray as xr
import os

os.chdir("/Users/mackenfr/scripts")

from plotting_funcs import plot_global

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

from scipy.interpolate import griddata

#%% define functions

def plot_polar_stereo(lon,lat,value,cmap,vmin,vmax,title):
    
    fig = plt.figure()

    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    ax.coastlines() 
    ax.gridlines(linestyle='--')
    #ax.set_boundary(circle, transform=ax.transAxes)

    ct=ax.pcolormesh(lon,lat,value,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading="auto")

    #colorbar 
    cb = plt.colorbar(ct,orientation="vertical")
    cb.set_label('m')
    
    plt.title(title)
        
def bin_ocean_depths_UVic(depth_m, lower_bounds):
    #Takes depth of ocean in metres (bathymetry) and bins it for UVic G_kmt.nc input file 
#   UVic kmt input requires depth to be classified into 19 levels

#   inputs: 
#   depth_m = the depth in metres to be classified
#   lower_bounds = the bottom (in meters) of each level; must be an array with 19 elements
#   level 0 = not ocean / very shallow

#   outputs:
#   bin_level = number 0-19 indicating ocean depth level for UVic model

#   translated from Matlab function bin_ocean_depths_UVic.m  03/06/2022

    if(depth_m >= 0):
        bin_level = 0
    elif(depth_m >= -lower_bounds[1]):
        bin_level = 1
    elif(depth_m >= -lower_bounds[2]):
        bin_level = 2
    elif(depth_m >= -lower_bounds[3]):
        bin_level = 3
    elif(depth_m >= -lower_bounds[4]):
        bin_level = 4
    elif(depth_m >= -lower_bounds[5]):
        bin_level = 5
    elif(depth_m >= -lower_bounds[6]):
        bin_level = 6
    elif(depth_m >= -lower_bounds[7]):
        bin_level = 7
    elif(depth_m >= -lower_bounds[8]):
        bin_level = 8
    elif(depth_m >= -lower_bounds[9]):
        bin_level = 9
    elif(depth_m >= -lower_bounds[10]):
        bin_level = 10
    elif(depth_m >= -lower_bounds[11]):
        bin_level = 11
    elif(depth_m >= -lower_bounds[12]):
        bin_level = 12
    elif(depth_m >= -lower_bounds[13]): 
        bin_level = 13
    elif(depth_m >= -lower_bounds[14]):
        bin_level = 14
    elif(depth_m >= -lower_bounds[15]):
        bin_level = 15
    elif(depth_m >= -lower_bounds[16]):
        bin_level = 16    
    elif(depth_m >= -lower_bounds[17]):
        bin_level = 17
    elif(depth_m >= -lower_bounds[18]):
        bin_level = 18    
    else:
        bin_level = 19    
        
    return bin_level

def write_netcdf_depth(filename,latdepth,new_lat,new_lon,new_var,varname,units):
    
    ds = Dataset(filename, 'w', format='NETCDF4')

    #time = ds.createDimension('time', None)
    lat = ds.createDimension('latitude', 100)
    lon = ds.createDimension('longitude', 100)
    depth = ds.createDimension('depth', 19)

    #times = ds.createVariable('time', 'f4', ('time',))
    lats = ds.createVariable('latitude', 'f4', ('latitude',))
    lons = ds.createVariable('longitude', 'f4', ('longitude',))
    depths = ds.createVariable('depth','f4', ('depth',))
    O_var_new = ds.createVariable(varname, 'f4', ('depth','latitude', 'longitude',))
    O_var_new.units = units

    lats[:] = new_lat
    lons[:] = new_lon
    depths[:] = latdepth

    O_var_new[:,:,:] = new_var

    ds.close()
    
    
def write_netcdf_Osur(filename,new_time,new_lat,new_lon,new_var,varname,units):
    
    ds = Dataset(filename, 'w', format='NETCDF4')

    time = ds.createDimension('time', 12)
    lat = ds.createDimension('latitude', 100)
    lon = ds.createDimension('longitude', 100)
    #depth = ds.createDimension('depth', 19)

    times = ds.createVariable('time', 'f8', ('time',))
    lats = ds.createVariable('latitude', 'f4', ('latitude',))
    lons = ds.createVariable('longitude', 'f4', ('longitude',))
    #depths = ds.createVariable('depth','f4', ('depth',))
    O_var_new = ds.createVariable(varname, 'f4', ('time','latitude', 'longitude',))
    O_var_new.units = units

    lats[:] = new_lat
    lons[:] = new_lon
    times[:] = new_time

    O_var_new[:,:,:] = new_var

    ds.close()    
    
#%% path to save files
uvic_filedir_lig = '/Users/mackenfr/scripts/pism_to_uvic/outputs/'

    


#%% UVic orignal files

uvic_filedir = '/Users/mackenfr/scripts/pism_to_uvic/inputs/uvic/'

ffilename1 = uvic_filedir + 'L_elev.nc'
ffilename2 = uvic_filedir + 'G_mskt.nc'
ffilename3 = uvic_filedir + 'G_mskhreg.nc'
ffilename4 = uvic_filedir + 'G_kmt.nc'
ffilename5 = uvic_filedir + 'L_rivers.nc'
ffilename6 = uvic_filedir + 'L_ice.nc'
ffilename7 = uvic_filedir + 'G_grid.nc'

uvic_elev_data = xr.open_dataset(ffilename1)
uvic_mskt_data = xr.open_dataset(ffilename2)
uvic_mskhreg_data = xr.open_dataset(ffilename3)
uvic_kmt_data = xr.open_dataset(ffilename4)
uvic_rivers_data = xr.open_dataset(ffilename5)
uvic_ice_data = xr.open_dataset(ffilename6)
uvic_grid_data = xr.open_dataset(ffilename7)

uvic_elev = uvic_elev_data.variables['L_elev'][:]
uvic_mskt = uvic_mskt_data.variables['G_mskt'][:]
uvic_mskhreg = uvic_mskhreg_data.variables['G_mskhreg'][:]
uvic_kmt = uvic_kmt_data.variables['G_kmt'][:]
uvic_depth = uvic_grid_data.variables['depth'][:] #depth of t grid (upper bounds of each level)
uvic_depthW = uvic_grid_data.variables['depth_W'][:] #depth of W grid = mid point of depth_t(i) - depth_t(i-1)
uvic_rivers = uvic_rivers_data.variables['L_rivdis'][:]
uvic_icefra = uvic_ice_data.variables['L_icefra'][:]
uvic_icethk = uvic_ice_data.variables['L_icethk'][:]
uvic_time = uvic_ice_data.variables['time'][:]

uvic_lat = uvic_elev_data.variables['latitude'][:]
uvic_lon = uvic_elev_data.variables['longitude'][:]


xx, yy = np.meshgrid(uvic_lon, uvic_lat)

#%%


uvic_kmt[40:50,60:70]=10

plot_global(xx, yy, uvic_kmt, "viridis", 0, 19, "title")
#%%
filename = uvic_filedir_lig + 'G_kmt_table.nc'
ds = Dataset(filename, 'w', format='NETCDF4')

#time = ds.createDimension('time', None)
lat = ds.createDimension('latitude', 100)
lon = ds.createDimension('longitude', 100)

#times = ds.createVariable('time', 'f4', ('time',))
lats = ds.createVariable('latitude', 'f4', ('latitude',))
lons = ds.createVariable('longitude', 'f4', ('longitude',))
G_kmt_new = ds.createVariable('G_kmt', 'f4', ('latitude', 'longitude',))

lats[:] = uvic_lat
lons[:] = uvic_lon

G_kmt_new[ :, :] = uvic_kmt

ds.close()

#%% fill in values for all ocean input files where kmt differs
cmap = "cool"
filelist = glob.glob('/Users/mackenfr/scripts/pism_to_uvic/inputs/uvic/O_*')
filelist.sort()

#ignore these files, they don't need to be changed for these model runs
# or they are monthly surface files and need to be dealt with separately
vars2ignore = ['O_cfc','O_diffac','O_tau','O_tidenrg','O_tempsur','O_salsur']
depth_diff = uvic_kmt_lig - uvic_kmt
indices_depth_diff = np.where((depth_diff>0) & (lig_mskt_diff==0)) #where depth of ocean has changed, but was already ocean
indices_ocean = np.where((lig_mskt_diff>=0.99)) #indices where land changed to ocean
indices_land = np.where((lig_mskt_diff<=-0.99)) #ocean to land


#load in all files from original data directory

#varlist=[]
for f in filelist:
    
    t1 = f.split('/')[-1]
    var = t1.split('.')[0]
    #varlist.append(var)
    
    if (var in vars2ignore):
        #ignore this file
        continue

    else:
        O_data = Dataset(f)
        uvic_var = O_data.variables[var][:]
        fillval = O_data.variables[var].fill_value
        units = O_data.variables[var].units
        
        #find differences in land mask and kmt and fill in new ocean cells 
        #with mean value
        
        #loop through indices (these are the cells that have changed from land to ocean)
        # find depth change
        for i in range(len(indices_ocean[0])):
            
            yi = indices_ocean[0][i] #lat
            xi = indices_ocean[1][i] #lon
            
            level_change = depth_diff[yi,xi]
            #if change is negative or zero, do nothing
            #if change is positve, add additional values at those depths
            if level_change<=0:
                continue
            else:
                for d in range(0,level_change):
                    current_dlevel = uvic_kmt[yi,xi] + d
                    
                    #take mean of variable south of -60 at depth
                    tmp_mean = np.nanmean(uvic_var[current_dlevel,0:16,:])
                    #print(tmp_mean)                    
                    uvic_var[current_dlevel,yi,xi] = tmp_mean
        
        import pdb; pdb.set_trace()
            
        #extend files where depth has changed (but cells were already ocean)  
# =============================================================================
#         for i in range(len(indices_depth_diff[0])):
#             
#             yi = indices_depth_diff[0][i] #lat
#             xi = indices_depth_diff[1][i] #lon    
#             
#             level_change = depth_diff[yi,xi]
#             
#             for d in range(0,level_change):
#                 current_dlevel = uvic_kmt[yi,xi] + d
#                 
#                 #take mean of variable south of -60 at depth
#                 tmp_mean = np.nanmean(uvic_var[current_dlevel,0:16,:])
#                 #print(tmp_mean)                    
#                 uvic_var[current_dlevel,yi,xi] = tmp_mean
#                 
# =============================================================================
                
        #find cells that have changed from ocean to land and fill in
        # !!! this needs to be done last, otherwise fillval is incorporated in tmp_mean
        #for i in range(len(indices_land[0])):
                    
            #yi = indices_land[0][i] #lat
            #xi = indices_land[1][i] #lon    
                    
            #uvic_var[:,yi,xi] = fillval
       # plot_polar_stereo(xx, yy, uvic_var[0,:,:], cmap, vmin, vmax, 'title')            
  
        #write new file 
        new_filename = uvic_filedir_lig + var + '_LIG.nc'
        write_netcdf_depth(new_filename,uvic_depth,uvic_lat,uvic_lon,uvic_var,var,units)
                 
#%% need to do surface quantities separately

filelist = glob.glob('/Users/mackenfr/scripts/pism_to_uvic/inputs/uvic/O_*sur.nc')
filelist.sort()

indices_ocean = np.where((lig_mskt_diff>0.99)) #indices where land changed to ocean 
indices_land = np.where((lig_mskt_diff<-0.99)) #ocean changed to land

vars2ignore = []

for f in filelist:
    
    t1 = f.split('/')[-1]
    var = t1.split('.')[0]
    
    if (var in vars2ignore):
        #ignore this file
        continue

    else:
        O_data = Dataset(f)
        uvic_var = O_data.variables[var][:]
        
        fillval = O_data.variables[var].fill_value
        times = O_data.variables['time'][:]
        units = O_data.variables[var].units
    
        #find differences in land mask and kmt and fill in new ocean cells 
        #with mean value
        
        #loop through indices (these are the cells that have changed land/ocean)
        for i in range(len(indices_ocean[0])):
            
            yi = indices_ocean[0][i] #lat
            xi = indices_ocean[1][i] #lon
            
            for t in range(0,12):
                #preserve monthly mean
                #take mean of variable south of -60 at surface
                tmp_mean = np.nanmean(uvic_var[t,0:16,:])
                #print(tmp_mean)                    
                uvic_var[t,yi,xi] = tmp_mean
                    
        #for i in range(len(indices_land[0])):
            
            #yi = indices_land[0][i] #lat
            #xi = indices_land[1][i] #lon
            
            #put in fill value everywhere
            #uvic_var[:,yi,xi] = fillval
        #plot_polar_stereo(xx, yy, uvic_var[0,:,:], cmap, vmin, vmax, 'title')            
                        
        #write new file 
        new_filename = uvic_filedir_lig + var + '_LIG.nc'
        write_netcdf_Osur(new_filename,times,uvic_lat,uvic_lon,uvic_var,var,units)
        


#%% plot pism ice regrid

cmap = cm.get_cmap("PiYG")
vmin=0
vmax=1000

plot_polar_stereo(xx,yy,uvic_icefra[21,:,:],cmap,vmin,vmax,'PD ice uvic')

plot_polar_stereo(xx,yy,pd_icethk_uvicgrid,cmap,vmin,vmax,'PD ice regrid')
plot_polar_stereo(xx,yy,lig_icethk_uvicgrid,cmap,vmin,vmax,'LIG ice regrid')


#%% plot usurf anomaly

cmap = cm.get_cmap("PiYG")
vmin=0
vmax=5

plot_polar_stereo(xx,yy,lig_anom_icethk,cmap,vmin,vmax,'icethk anom')

#%% plot new elev
cmap = cm.get_cmap("PiYG")
vmin=1
vmax=2000
plot_polar_stereo(xx,yy,uvic_elev_lig,cmap,vmin,vmax,'L_elev LIG')

lig_bin = (uvic_elev_lig/uvic_elev_lig)
elev_diff = uvic_elev_lig + (1-lig_mask_uvicgrid_glob )

#plot_polar_stereo(xx,yy,elev_diff,cmap,vmin,vmax,'L_elev diff from mskt')



#%% plot land mask
cmap = cm.get_cmap("PiYG")
vmin=0
vmax=4

#plot_polar_stereo(lig_lon,lig_lat,lig_mask[0,:,:],cmap,vmin,vmax,'LIG land mask')
vmin=0
vmax=1

plot_polar_stereo(xx,yy,lig_mask_uvicgrid,cmap,vmin,vmax,'LIG mask regrid')
#%% difference between G_mskt
cmap = cm.get_cmap("PiYG")

vmin=0
vmax=4
plot_polar_stereo(xx,yy,uvic_mskhreg_lig,cmap,vmin,vmax,'G_mskhreg')

vmin=-1
vmax=1
plot_polar_stereo(xx,yy,lig_mskt_diff,cmap,vmin,vmax,'G_mskt diff LIG')

vmin=0
vmax=19
plot_polar_stereo(xx,yy,uvic_kmt_lig,cmap,vmin,vmax,'G_kmt')

#%% temp - delete or move up

def plot_global_grid(dataset, x, y, cmap,vmin,vmax, title):    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    fig.suptitle(title, fontsize="xx-large", y = 0.93)           
    ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    ax.coastlines() 
    #ax.patch.set_alpha(0)   
    ax.set_facecolor('white')
    #ax.gridlines(linestyle='--')

    ct=ax.pcolormesh(x, y,dataset[0,:,:],transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap)#, shading="auto")

    plt.title(title)
  
