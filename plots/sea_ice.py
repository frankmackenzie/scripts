#%%

import matplotlib as mpl
import matplotlib.pyplot as plt
import cftime
import xarray as xr
import cartopy.crs as ccrs
import os
import seaborn as sns
import numpy as np
import matplotlib.colors as cm
from matplotlib.lines import Line2D
import pandas as pd

from plot_polar import plot_polar_grid
os.chdir("/Users/mackenfr/scripts")

#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
#%%
def convert_time(data):
    newdata = data.copy()
    newdata["time"] = data["time"] - data["time"][0]
    units = 'common_years since 1970-1-1 00:00:00'
    time_365 = cftime.num2date(newdata.time, units, '365_day')
    newtime = xr.DataArray(time_365, coords = [time_365], dims = 'time', name = 'data')
    newdata["time"] = newtime
    newdata["time"] = newdata.indexes["time"].to_datetimeindex()

    return(newdata)

#just a note that the averages for eqb and trs were taken by different methods (ydaymean of 10year runstep) than the amn and prd ones (ensmean of 100year run) - to be amended
def import_data(path_to_data, config_names=["PRD-PI", "AMN-PI","TRS-PI", "EQB-PI"]):# 
    config_names = config_names # "PRD-LIG","AMN-LIG"
    co2s = ["280"]#, "300", "350", "400", "450"]#
    tavgs = []
    names = []
    co2_all=[]
    ices = []
    orbs = []
    for config in config_names:
        for co2 in co2s:
            if config == "PRD-PI" or config == "AMN-PI":
                filename_tavg = path_to_data + "/" + config + "/" + co2 + "/sea_ice_annual.nc"
            if config == "EQB-PI" or config == "TRS-PI":
                filename_tavg = path_to_data + "/" + config + "/" + co2 + "/sea_ice_mean.nc"
            data = xr.open_dataset(filename_tavg, decode_times=False)#, engine = "netcdf4")
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
# %%
tavgs,ref_df = import_data(path_to_data)
# %% mask out all land
prd = xr.open_dataset(path_to_data+"PRD-PI/280/tavg.nc", decode_times=False, engine = "netcdf4")
amn = xr.open_dataset(path_to_data+"AMN-PI/280/tavg.nc", decode_times=False, engine = "netcdf4")
trs = xr.open_dataset(path_to_data+"TRS-PI/280/tavg.nc", decode_times=False, engine = "netcdf4")
eqb = xr.open_dataset(path_to_data+"EQB-PI/280/tavg.nc", decode_times=False, engine = "netcdf4")
prd_mask = prd["G_mskt"].rename({"latitude":"lat", "longitude":"lon"})
amn_mask = amn["G_mskt"].rename({"latitude":"lat", "longitude":"lon"})
trs_mask = trs["G_mskt"].rename({"latitude":"lat", "longitude":"lon"})
eqb_mask = eqb["G_mskt"].rename({"latitude":"lat", "longitude":"lon"})

area = prd["G_areaT"].rename({"latitude":"lat", "longitude":"lon"})
max_area = area.max()
time = tavgs[0].time
weights = (area/max_area).expand_dims(dim={"time":time})
#%%
#temp - fix this when fix averaging method
tavgs[10] = tavgs[10].rename({"latitude":"lat", "longitude":"lon"})
tavgs[15] = tavgs[15].rename({"latitude":"lat", "longitude":"lon"})

#%% rather than sum of icefra, total area of cells where icefra >=0.15
#mask out prd and amn to compare bc then we are only comparing cells that are ocean in all. don't need to mask out trs/int and eqb as there is no extra land in those
#possible to also compare absolute sea ice for each by masking out each ice sheet

icefra = []
i=0
for data in tavgs:
    i+=1
    ice = data["O_icefra"][:,:,:]
    mask = xr.where(((prd_mask[:,:]==1) & (amn_mask[:,:]==1)), ice, 0) 
    ice_edge = xr.where(mask>=0.15, 1, 0) 
    icefra.append(ice_edge)

#%% total sea ice over time
nh_area = []
sh_area = []

for data in icefra:
    icefra_area_sh = data*area[:25,:]
    icefra_area_nh = data*area[75:,:]
    sum_sh = icefra_area_sh.sum(dim=["lat", "lon"])
    sum_nh = icefra_area_nh.sum(dim=["lat", "lon"])
    sh_area.append(sum_sh)
    nh_area.append(sum_nh)
#%%
total_area = ref_df[["config", "co2"]].copy()
total_area["SH"] = sh_area
total_area["NH"] = nh_area

# %% SH over time
cmap = mpl.colormaps['plasma_r']
colors = cmap(np.linspace(0.2, 1, 5))

time_ref = convert_time(tavgs[0])
new_time = pd.to_datetime(time_ref["time"])#.strftime("%d-%m")
#%%
time_dm =new_time.strftime("%d-%m")
fig, ax = plt.subplots()

color=0
for i, label in zip(sh_area[:5], ref_df.co2[:5]):
    plt.plot(time_dm, i, label=label, color=colors[color])
    color+=1
plt.legend()
plt.xticks(time_dm[::3], rotation="vertical")
plt.title("SH ice area PRD")
color=0
fig, ax = plt.subplots()
for i, label in zip(sh_area[5:10], ref_df.co2[5:10]):
    plt.plot(time_dm, i, label=label, color=colors[color])
    color+=1
plt.legend()
plt.xticks(time_dm[::3], rotation="vertical")
plt.title("SH ice area AMN")

fig, ax = plt.subplots()

color=0
for i, label in zip(sh_area[10:15], ref_df.co2[10:15]):
    plt.plot(time_dm, i, label=label, color=colors[color])
    color+=1
plt.legend()
plt.xticks(time_dm[::3], rotation="vertical")
plt.title("SH ice area TRS")
color=0
fig, ax = plt.subplots()
for i, label in zip(sh_area[15:], ref_df.co2[15:]):
    plt.plot(time_dm, i, label=label, color=colors[color])
    color+=1
plt.legend()
plt.xticks(time_dm[::3], rotation="vertical")
plt.title("SH ice area EQB")
#%%

# %% 
total_area['max_SH'] = total_area['SH'].apply(lambda da: float(da.max()))
total_area['max_SH_date'] = total_area['SH'].apply(lambda da: new_time[da.argmax().item()])
total_area['max_SH_index'] = total_area['SH'].apply(lambda da: da.argmax().item())

total_area['max_NH'] = total_area['NH'].apply(lambda da: float(da.max()))
total_area['max_NH_date'] = total_area['NH'].apply(lambda da: new_time[da.argmax().item()])
total_area['max_NH_index'] = total_area['NH'].apply(lambda da: da.argmax().item())

total_area['min_SH'] = total_area['SH'].apply(lambda da: float(da.min()))
total_area['min_SH_date'] = total_area['SH'].apply(lambda da: new_time[da.argmin().item()])
total_area['min_SH_index'] = total_area['SH'].apply(lambda da: da.argmin().item())

total_area['min_NH'] = total_area['NH'].apply(lambda da: float(da.min()))
total_area['min_NH_date'] = total_area['NH'].apply(lambda da: new_time[da.argmin().item()])
total_area['min_NH_index'] = total_area['NH'].apply(lambda da: da.argmin().item())
#%% date of annual max and min SI extent


#%%max and min sea ice extent
fig, ax = plt.subplots()

sns.scatterplot(x="co2", y="max_SH", data = total_area, hue = "config")
sns.lineplot(x="co2", y="max_SH", data = total_area, hue = "config", legend=False)

fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="max_NH", data = total_area, hue = "config")
sns.lineplot(x="co2", y="max_NH", data = total_area, hue = "config", legend=False)

fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="min_SH", data = total_area, hue = "config")
sns.lineplot(x="co2", y="min_SH", data = total_area, hue = "config", legend=False)

fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="min_NH", data = total_area, hue = "config")
sns.lineplot(x="co2", y="min_NH", data = total_area, hue = "config", legend=False)

#%% regional ice area
basins = {"Atlantic":1, "Pacific": 2, "Indian": 3}
regions = prd["G_mskhr"].rename({"latitude":"lat", "longitude":"lon"})
atl_area = []
pac_area = []
ind_area = []
amn_area = []

area_mask_atl = area.where(regions==1)
area_mask_pac = area.where(regions==2)
area_mask_ind = area.where(regions==3)

for i, data in enumerate(icefra):
    index = total_area["max_SH_index"][i]
    data_masked_atl = data[:,:,index].where(regions==1)
    data_masked_pac = data[:,:,index].where(regions==2)
    data_masked_ind = data[:,:,index].where(regions==3)
    icefra_area_atl = data_masked_atl*area_mask_atl
    icefra_area_pac = data_masked_pac*area_mask_pac 
    icefra_area_ind = data_masked_ind*area_mask_ind
    icefra_area_amn = icefra_area_pac[:,61:78]
    sum_atl = icefra_area_atl.sum(dim=["lat", "lon"])
    sum_pac = icefra_area_pac.sum(dim=["lat", "lon"])
    sum_ind = icefra_area_ind.sum(dim=["lat", "lon"])
    sum_amn = icefra_area_amn.sum(dim=["lat", "lon"])

    atl_area.append(float(sum_atl))
    pac_area.append(float(sum_pac))
    ind_area.append(float(sum_ind))
    amn_area.append(float(sum_amn))

regional_area = pd.DataFrame({
                            "co2": ref_df["co2"],
                            "ice sheet": ref_df["ice"],
                            "Atlantic":atl_area,
                            "Pacific": pac_area,
                            "Indian": ind_area,
                            "Amundsen": amn_area
                            })    
#%%
fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="Atlantic", data = regional_area[["co2", "ice sheet", "Atlantic"]], hue = "ice sheet")
sns.lineplot(x="co2", y="Atlantic", data = regional_area[["co2", "ice sheet", "Atlantic"]], hue = "ice sheet", legend=False)
ax.set_ylabel("Max ice area")
ax.set_title("Atlantic")

fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="Pacific", data = regional_area[["co2", "ice sheet", "Pacific"]], hue = "ice sheet")
sns.lineplot(x="co2", y="Pacific", data = regional_area[["co2", "ice sheet", "Pacific"]], hue = "ice sheet", legend=False)
ax.set_ylabel("Max ice area")
ax.set_title("Pacific")

fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="Indian", data = regional_area[["co2", "ice sheet", "Indian"]], hue = "ice sheet")
sns.lineplot(x="co2", y="Indian", data = regional_area[["co2", "ice sheet", "Indian"]], hue = "ice sheet", legend=False)
ax.set_ylabel("Max ice area")
ax.set_title("Indian")
#%%
fig, ax = plt.subplots()
sns.scatterplot(x="co2", y="Amundsen", data = regional_area[["co2", "ice sheet", "Amundsen"]], hue = "ice sheet")
sns.lineplot(x="co2", y="Amundsen", data = regional_area[["co2", "ice sheet", "Amundsen"]], hue = "ice sheet", legend=False)
ax.set_ylabel("Max ice area")
ax.set_title("Amundsen")
# %% south pole icefra

def plot_extent(title, threshhold, minmax,  proj=ccrs.SouthPolarStereo()):

    proj=ccrs.SouthPolarStereo()
    fig, axs = plt.subplots(2,2, figsize=(8,8), subplot_kw={'projection': proj})
    fig.suptitle(title, fontsize="xx-large", y=0.95)

    lons = tavgs[0]["lon"]
    lats = tavgs[0]["lat"]

    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = -50
        refstring = "SH"
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 45
        north_bound = 90
        refstring = "NH"

    string = minmax + "_" + refstring + "_index"

    for i, dataset in enumerate(tavgs): 
        index = total_area[string][i]
        data = dataset["O_icefra"][index, :, :]
        if i <= 4:
            j=(0,0)
        elif i <= 9:
            j=(0,1)
        elif i <= 14:
            j=(1,0)
        else:
            j=(1,1)

        axs[j].set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        axs[j].coastlines()
        axs[j].gridlines(linestyle='--', draw_labels=False, color="gray")
        axs[j].contour(lons, lats, data, transform=ccrs.PlateCarree(),
                        colors=colors[i%5],
                        levels=[threshhold,1] 
                    )
    axs[0,0].set_title("PRD")
    axs[0,1].set_title("AMN")
    axs[1,0].set_title("TRS")
    axs[1,1].set_title("EQB")

    axs[0,0].pcolormesh(lons, lats, tavgs[0]["O_icethk"][total_area[string][0], :, :], transform=ccrs.PlateCarree(), cmap="Blues", vmin=0, vmax=5)
    axs[0,1].pcolormesh(lons, lats, tavgs[5]["O_icethk"][total_area[string][5], :, :], transform=ccrs.PlateCarree(), cmap="Blues", vmin=0, vmax=5)
    axs[1,0].pcolormesh(lons, lats, tavgs[10]["O_icethk"][total_area[string][10], :, :], transform=ccrs.PlateCarree(), cmap="Blues", vmin=0, vmax=5)
    axs[1,1].pcolormesh(lons, lats, tavgs[15]["O_icethk"][total_area[string][15], :, :], transform=ccrs.PlateCarree(), cmap="Blues", vmin=0, vmax=5)

    handles = [Line2D([0], [0], color=colors[i], lw=2) for i in range(5)]
    legend_labels = [str(label) for label in ref_df["co2"].unique()]
    axs[0,0].legend(handles, legend_labels, title="CO$_{2}$", loc="lower left", fontsize="small")
    
    cb = fig.colorbar(axs[0, 0].collections[-1], ax=axs, orientation='vertical', fraction=0.02, pad=0.06)
    cb.set_label("ice thickness at 280 ppm (m)", fontsize=14)
#%%
plot_extent("max sea ice edge", 0.15, "max", ccrs.SouthPolarStereo())
#plot_extent("min sea ice edge", 0.15, "min", ccrs.SouthPolarStereo())
#plot_extent("max sea ice edge", 0.15, "max", ccrs.NorthPolarStereo())
#plot_extent("min sea ice edge", 0.15, "min", ccrs.NorthPolarStereo())

# %%

def plot_polar_grid(var, cmap, vmin, vmax, proj = ccrs.SouthPolarStereo(), title="", anom=None, minmax="max"): 
    if proj == ccrs.SouthPolarStereo():
        south_bound = -90
        north_bound = -60
        figsize=(15,12)
        refstring = "SH"
    elif proj == ccrs.NorthPolarStereo():
        south_bound = 45
        north_bound = 90
        figsize=(15,12)
        refstring = "NH"

    if anom=="picontrol":
        suptitle=title + " anom from PRD 280"
        anom_idx = 0

    elif anom=="config":
        suptitle=title + " anom from PRD"
        anom_idx = [(i % 5) for i in range(len(ref_df.config))]
    elif anom=="co2":
        suptitle=title + " anom from 280 ppm"
        anom_idx = [(i // 5) * 5 for i in range(len(ref_df.config))]
    else:
        suptitle=title

    string = minmax + "_" + refstring + "_index"
    index = [total_area[string][i] for i in range(len(ref_df.config))]
    #mask pi
    if anom:
        anomalies = [xr.where(prd_mask==1, tavgs[i][var][index[i],:, :] - tavgs[anom_idx[i]][var][index[anom_idx[i]],:, :],0)  for i in range(len(tavgs))]
    else:
        anomalies = [tavgs[i][var][index[i],:, :] for i in range(len(tavgs))]
  
    fig, axs = plt.subplots(4, 5, figsize=figsize, subplot_kw={'projection': proj})
    fig.suptitle(suptitle, fontsize="xx-large", y=0.99)

    lons = tavgs[0]["lon"]
    lats = tavgs[0]["lat"]

    for i, ax in enumerate(axs.flat):
        data = anomalies[i]
            
        ax.set_extent([-180, 180,south_bound, north_bound], ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines(linestyle='--', draw_labels=False, color="gray")
        ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(),
                        vmin=vmin, vmax=vmax, 
                        cmap=cmap, shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
    
    for n, col in enumerate(ref_df.co2.unique()):
        axs[0, n].set_title(col)

    for n, row in enumerate(ref_df.config.unique()):
        axs[0, n].set_ylabel(row)
        fig.text(0.12, 0.8 - (n * 0.19), row, va='center', ha='right', fontsize=12, rotation=90)

    cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.02, pad=0.04)
    cb.set_label(f"{tavgs[0][var].attrs.get('long_name', var)} ({tavgs[0][var].attrs.get('units', '')})", fontsize=14)
#%%
plot_polar_grid("O_icethk", "Blues", 0, 1, proj = ccrs.SouthPolarStereo(), title="Annual max ice thickness", anom=None, minmax="max")
plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.SouthPolarStereo(), title="Annual max ice thickness ", anom="co2", minmax="max")
plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.SouthPolarStereo(), title="Annual max ice thickness", anom="config",minmax="min")

plot_polar_grid("O_icethk", "Blues", 0, 1, proj = ccrs.SouthPolarStereo(), title="Annual min ice thickness", anom=None, minmax="min")
plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.SouthPolarStereo(), title="Annual min ice thickness", anom="co2", minmax="min")
plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.SouthPolarStereo(), title="Annual min ice thickness", anom="config", minmax="min")

# plot_polar_grid("O_icethk", "Blues", 0, 1, proj = ccrs.NorthPolarStereo(), title="Annual max ice thickness", anom=None, minmax="max")
# plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.NorthPolarStereo(), title="Annual max ice thickness", anom="co2", minmax="max")
# plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.NorthPolarStereo(), title="Annual max ice thickness ", anom="config", minmax="max")

# plot_polar_grid("O_icethk", "Blues", 0, 1, proj = ccrs.NorthPolarStereo(), title="Annual min ice thickness", anom=None, minmax="min")
# plot_polar_grid("O_icethk", "managua_r", -1, 1, proj = ccrs.NorthPolarStereo(), title="Annual min ice thickness", anom="co2", minmax="min")
# plot_polar_grid("O_icethk", "managua_r", -0.1, 0.1, proj = ccrs.NorthPolarStereo(), title="Annual min ice thickness (note scale)", anom="config",minmax="min")

# %%
plot_polar_grid("O_icethk", "plasma", 0, 12, proj = ccrs.SouthPolarStereo(), title="Annual max ice thickness", anom=None, minmax="max")
plot_polar_grid("O_icethk", "plasma", 0, 12, proj = ccrs.SouthPolarStereo(), title="Annual min ice thickness", anom=None, minmax="min")

#%%
refstring = "SH"
subset_280 = total_area.loc[total_area['co2']=="280"]
string = "max_" + refstring + "_index"
tavgs_subset = [tavgs[i] for i in list(range(0,20,5))]
index = [total_area[string][i] for i in list(range(20))]

values = [tavgs[i]["O_icefra"][index[i],:, :] for i in list(range(0,20,5))]

masks = [prd_mask, amn_mask, trs_mask, eqb_mask]

fig, axs = plt.subplots(2, 2, figsize=(8,8), subplot_kw={'projection': ccrs.SouthPolarStereo()})
fig.suptitle("280ppm", fontsize="xx-large", y=0.99)

lons = tavgs[0]["lon"]
lats = tavgs[0]["lat"]

for i, ax in enumerate(axs.flat):
    
    data = values[i]
    #mask out land
    sea_ice = xr.where(masks[i]==1, data, 0)
    ax.set_extent([-180, 180,-90, -45], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    ax.pcolormesh(lons, lats, sea_ice, transform=ccrs.PlateCarree(),
                    vmin=0, vmax=1, 
                    cmap=cmap, shading="auto")#, norm=colors.SymLogNorm(linthresh=0.000001, linscale=1, vmin=vmin, vmax=vmax))
    ax.set_title(ref_df.config.unique()[i])
cb = fig.colorbar(axs[0, 0].collections[1], ax=axs, orientation='horizontal', fraction=0.02, pad=0.04)
cb.set_label(f"sea ice fraction", fontsize=14)
# %%
