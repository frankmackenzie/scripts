
import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import seaborn as sns
#%%
path_to_data = "/Users/mackenfr/proxy_compare_2025/raw/"
capron_raw = pd.read_csv(path_to_data + "capron_2017.csv", header=12)

def getT(data, longitude, latitude):
    # extract array of lats and lons
    lats=ref_data["latitude"][:]
    lons=ref_data["longitude"][:]   
    lons_adj = np.where(lons>180, lons-360, lons)
    #find values nearest to given coords
    lat_index=(np.abs(lats - latitude)).argmin()
    lon_index=(np.abs(lons_adj - longitude)).argmin()
    T = data[lat_index, lon_index]
    return(float(T))

#%%
sat = (capron_raw[["Station",
                    "Latitude",
                    "Longitude",
                    "Surface air temperature [°C]",
                    "Type",
                    "127 ka Median PIAn [°C]",
                    "127 ka 2s PIAn [°C]"]]
                    .loc[~pd.isnull(capron_raw["Surface air temperature [°C]"])&(capron_raw["Latitude"]<0)]
                    .rename(columns = {"Station": "Site",
                                       "127 ka Median PIAn [°C]": "SAT anom",
                                       "127 ka 2s PIAn [°C]": "SAT anom error"})
)

# %% load uvic data

path_to_uvic = "/Users/mackenfr/proxy_compare_2025/processed_uvic/"
path_to_ref_data = "/Users/mackenfr/uvic_outputs/"

ref_data = xr.open_dataset(path_to_ref_data + "PRD-PI/280/tavg.nc" , decode_times=False)
mask = xr.open_dataset(path_to_ref_data+"masks.nc", decode_times=False, engine = "netcdf4")


icesheets = ["prdpi","prdlig", "amnlig", "eqblig"]
ref_is = []
uvic = []
for ice in icesheets:
    filename = ice + "annual.nc"
    data = xr.open_dataset(path_to_uvic + filename, decode_times=True)
    uvic.append(data)
    ref_is.append(ice.upper())
# %%
for i, ice in enumerate(ref_is):
    sat[ice] = sat.apply(lambda row: getT(uvic[i]["A_sat"], row["Longitude"],row["Latitude"]), axis=1)
    if i>0:
        sat[ice +" anom"] = sat[ice]-sat["PRDPI"]
# %%

norm = plt.Normalize(-10, 10)
lons = ref_data["longitude"]
lats = ref_data["latitude"]
control = uvic[0]["A_sat"]
for i, ice in enumerate(ref_is[1:]):
    data = uvic[i+1]["A_sat"] 
    masked_data = xr.where((mask["G_mskt"][i+1]==0), data, np.nan) 
    anom = masked_data - control
    fig, ax = plt.subplots(subplot_kw={'projection':ccrs.SouthPolarStereo()})
    ax.set_extent([-180, 180,-90, -50], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
    pcm = ax.pcolormesh(lons, lats, anom, cmap = "Spectral_r",transform= ccrs.PlateCarree(),norm=norm)

    sns.scatterplot(data = sat, 
                x = sat["Longitude"],
                y= sat["Latitude"], 
                hue = sat["SAT anom"], 
                marker = "o",
                palette="Spectral_r",
                transform= ccrs.PlateCarree(),
                hue_norm=norm,
                legend=False,
                edgecolors="black",
                s=80)


    cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.05, pad=0.04)
    ax.set_title(ice + " SAT anom")
# %%

outpath = "/Users/mackenfr/proxy_compare_2025/processed_proxy/"
sat.to_csv(outpath + "sat.csv")
# %%
