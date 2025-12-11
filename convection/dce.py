#%%
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from calc_ovt import calc_ovt
import gsw
import cartopy.crs as ccrs

#%%
path_to_data = "/Users/mackenfr/uvic_outputs/AMN-PI/400/" 
tavgs = xr.open_dataset(path_to_data + "dce.nc", decode_times=False)
ref_data = xr.open_dataset(path_to_data + "tavg.nc", decode_times=False)
new_dims = (ref_data[["G_dxu", "G_dzt", "G_mskhr"]]
            .expand_dims({"time": tavgs["time"]})
            .rename({"longitude": "lon", "latitude": "lat"}))
tavgs_expanded = tavgs.assign(new_dims)
#%%
def calc_ovt(data, region="Global"):
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    if region=="Global":
        #convert velocity to volume
        transport = data["O_velY"] * data["G_dxu"] * data["G_dzt"]
    elif region in basins:
        G_mskhr_aligned = data["G_mskhr"].sel(lat=data["latitude_V"], lon=data["longitude_V"], method="nearest")
        mask = data[["O_velY","G_dxu", "G_dzt"]].where(G_mskhr_aligned==basins[region])
        transport = (mask["O_velY"] * mask["G_dxu"] * mask["G_dzt"])

    zonalsum = transport.sum(dim="longitude_V")
    # integrate from bottom
    integ_layers_from_bottom = zonalsum.cumsum(dim="depth") - zonalsum.sum(dim="depth")
    # the result of the integration over layers is evaluated at the interfaces
    # with psi = 0 as the bottom boundary condition for the integration
    bottom_condition = xr.zeros_like(integ_layers_from_bottom.isel({"depth": -1}))
    psi_raw = xr.concat([integ_layers_from_bottom, bottom_condition], dim="depth")
    psi_raw = psi_raw.chunk({"depth": len(psi_raw["depth"])})  # need to rechunk to new size
    # convert kg.s-1 to Sv (1e6 m3.s-1)
    psi_Sv = psi_raw / 1.0e6
    return(psi_Sv)
# %%
amoc = calc_ovt(tavgs_expanded, region= "Atlantic")
amoc_max = amoc[:,4:,50:].max(dim=["latitude_V", "depth"])#below 380m (upper bound of 4th depth level) and above 0N - from  Weiffenbach 2024 (note it's the same if taken as whole atl)

mean_t = tavgs["O_temp"].mean(dim=["lat", "lon","depth"])
mean_sst = tavgs["O_temp"][:,0,:,:].mean(dim=["lat", "lon"])

mean_sal = tavgs["O_sal"].mean(dim=["lat", "lon","depth"])

df = pd.DataFrame()
df["time"]=tavgs["time"]
df["amoc"] = amoc_max.values
df["temp"] = mean_t.values
df["sst"] = mean_sst.values
df["sal"] = mean_sal.values

#%%
fig, ax = plt.subplots(figsize=(10,3))
fig.subplots_adjust(right=0.6)
ax.set_xlim(3000,7000)

ax1= ax.twinx()
ax2= ax.twinx()
ax3= ax.twinx()

sns.lineplot(x="time", y="amoc", data = df, ax=ax)
sns.lineplot(x="time", y="temp", data = df, ax= ax1, color="red")
sns.lineplot(x="time", y="sst", data = df, ax= ax2, color="green")
sns.lineplot(x="time", y="sal", data = df, ax= ax3, color="orange")

ax2.spines.right.set_position(("axes", 1.15))
ax3.spines.right.set_position(("axes", 1.3))


ax.yaxis.label.set_color("blue")
ax1.yaxis.label.set_color("red")
ax2.yaxis.label.set_color("green")
ax3.yaxis.label.set_color("orange")

# %%

def calc_mean(data, var, ts, basin):
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    G_mskhr = data["G_mskhr"][0,:,:]
    if basin=="Global":
        mean = data[var][ts,:,:,:].mean(dim=["lon"])
    elif basin in basins:
        mask = data[var].where(G_mskhr==basins[basin])
        mean = mask[ts,:,:,:].mean(dim=["lon"])
    return(mean)
#%%
def calculate_density(data, ts, basin, return_mean=True):
    lon = data["lon"]
    lat = data["lat"]
    SP = data["O_sal"][ts,:,:,:] #practical salinity
    pt = data["O_temp"][ts,:,:,:] #potential temp
    p = data["depth"] #sea pressure (approximate depth in m ~= pressure in dbar)

    SA = gsw.conversions.SA_from_SP(SP, p, lon, lat) #convert practical salinity into absolute salinity
    CT = gsw.conversions.CT_from_pt(SA, pt) #potential temperature to conservative temperature
    sigma0 = gsw.density.sigma0(SA, CT) #calculate potential density anomaly (anomaly from surface)
    basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
    G_mskhr = data["G_mskhr"][0,:,:]
    if basin=="Global":
        mean = sigma0.mean(dim=["lon"])
    elif basin in basins:
        mask = sigma0.where(G_mskhr==basins[basin])
        mean = mask[:,:,:].mean(dim=["lon"])
    if return_mean:
        return(mean)
    else:
        return(sigma0)


#%%

#%%
steps = {54: "amoc 'on'", 57: "amoc collapse?", 62: "amoc 'off'", 67: "end of amoc 'off'", 41: "amoc reestablishing?"}
for i in steps:
    atl_dens = calculate_density(tavgs_expanded, i, "Atlantic")
    atl_temp = calc_mean(tavgs_expanded, "O_temp", i, "Atlantic")
    atl_sal = calc_mean(tavgs_expanded, "O_sal", i, "Atlantic")
    atl_ovt = amoc[i,:,:]

    fig, ax = plt.subplots(3,1)
    lat=tavgs["lat"]
    depth=tavgs["depth"]
    depth1= atl_ovt["depth"]
    x, y = np.meshgrid(lat, depth)
    xx,yy= np.meshgrid(lat, depth1)
    levs = np.concatenate([np.linspace(26, 27, 2),np.linspace(27.7, 28, 6)])#np.linspace(26, 28, 21)#np.concatenate([np.linspace(26, 27.8, 9),np.linspace(27.85, 29, 12)])
    ax[0].invert_yaxis()
    cm1 = ax[0].pcolormesh(x, y,atl_temp, cmap="Spectral_r")
    c1 = ax[0].contour(x,y,atl_dens,levels=levs, colors="black")
    ax[0].clabel(c1, c1.levels, inline=True, fontsize=10)
    
    cb1 = fig.colorbar(cm1, orientation='vertical')
    cb1.set_label("temp")


    ax[1].invert_yaxis()
    cm2 = ax[1].pcolormesh(x, y,atl_sal, cmap="gist_ncar", vmin = 34, vmax=36)
    c1 = ax[1].contour(x,y,atl_dens,levels=levs, colors="black")
    ax[1].clabel(c1, c1.levels, inline=True, fontsize=10)
    
    cb2 = fig.colorbar(cm2, orientation='vertical')
    cb2.set_label("sal")


    ax[2].invert_yaxis()
    cm2 = ax[2].pcolormesh( xx,yy, atl_ovt, cmap="RdBu", vmin=-20, vmax=20)
    c1 = ax[2].contour( xx,yy, atl_ovt, colors="black")
    ax[2].clabel(c1, c1.levels, inline=True, fontsize=10)
    
    cb2 = fig.colorbar(cm2, orientation='vertical')
    cb2.set_label("sal")

    fig.suptitle(steps[i])
    
#return(fig)


# %%
for i in steps:
    fheat = tavgs_expanded["F_heat"][i, :, :]
    ventdep = tavgs_expanded["O_ventdep"][i, :, :]
    dens = calculate_density(tavgs_expanded, i, "Global", return_mean=False)

    fig, ax = plt.subplots(1,2,figsize = (10,4), subplot_kw={'projection': ccrs.PlateCarree()})
    lat=tavgs["lat"]
    lon=tavgs["lon"]

    ax[0].set_extent([-75, 30, -90, 90], ccrs.PlateCarree())
    ax[0].coastlines() 
    ax[0].gridlines(linestyle='--', draw_labels=False, color="gray")
    cm1 = ax[0].pcolormesh(lon, lat, fheat, cmap="jet", vmin=-200, vmax=200)

    ax[1].set_extent([-75, 30, 40, 90], ccrs.PlateCarree())
    ax[1].coastlines() 
    ax[1].gridlines(linestyle='--', draw_labels=False, color="gray")
    cm2 = ax[1].pcolormesh(lon, lat, ventdep, cmap="cubehelix", vmin=200, vmax=500)
    
    levs = np.concatenate([np.linspace(27, 27.7, 5)])#
    c1 = ax[0].contour(lon,lat, dens[0,:,:], levels=levs, colors="black")
    ax[0].clabel(c1, c1.levels, inline=True, fontsize=10)
    
    #c1 = ax[1].contour(lon,lat, ventdep, colors="white")
    #ax[1].clabel(c1, c1.levels, inline=True, fontsize=10)

    cb1 = fig.colorbar(cm1, orientation='horizontal', fraction=0.04)
    cb1.set_label("temp")

   

    cb2 = fig.colorbar(cm2, orientation='horizontal',fraction=0.04)
    cb2.set_label("sal")

    fig.suptitle(steps[i])

# %%
window = list(range(3000,7000))
