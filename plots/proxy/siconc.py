import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#%%
path_to_data = "/Users/mackenfr/proxy_compare_2025/raw/"
chadwick_raw = pd.read_csv(path_to_data + "Chadwick_2021.csv")
hadisst_raw = xr.load_dataset(path_to_data + "HadISST_ice.nc")

def get_sic(data, longitude, latitude):
    #extract array of lats and lons
    lats=ref_data["latitude"][:]
    lons=ref_data["longitude"][:]   
    #find values nearest to given coords
    lat_index=(np.abs(lats - latitude)).argmin()
    lon_index=(np.abs(lons - longitude)).argmin()
    siconc = data[lat_index, lon_index]
    return(float(siconc))

hadisst_pi = hadisst_raw["sic"][:360]
septembers = hadisst_pi.loc[{"time": hadisst_pi.time.dt.month == 9}] #select september
hadisst_mean = septembers.mean(dim="time")
#%%
#need to interpolate chadwick to get 127k
#do this for both age models

def add_row(group):
    new_row = group.iloc[[0]].copy()
    new_row["Age [ka BP]"] = 127.0
    new_row["SST"] = np.nan
    group.iloc[[-1]] = new_row
    return group

subset_lr04 = (chadwick_raw[["Event",
                            "Latitude",
                            "Longitude",
                            "Age [ka BP] (Age model, LR04 Lisiecki & Ra...)", 
                            "Sea ice con (9) [%] (Modern analog technique (MAT))"
                            ]]
                            .rename(columns={"Age [ka BP] (Age model, LR04 Lisiecki & Ra...)": "Age [ka BP]",
                                            "SST sum [¬∞C] (Modern analog technique (MAT))": "SST",
                                            "Event": "Site"})
)

subset_edc3 = (chadwick_raw[["Event",
                             "Latitude",
                             "Longitude",
                             "Age [ka BP] (Age model, EDC3 (EPICA Ice Do...)", 
                             "Sea ice con (9) [%] (Modern analog technique (MAT))"
                             ]]
                            .rename(columns={"Age [ka BP] (Age model, EDC3 (EPICA Ice Do...)": "Age [ka BP]",
                                        "SST sum [¬∞C] (Modern analog technique (MAT))": "SST",
                                        "Event": "Site"})
)
chadwick_interp_lr04 = (subset_lr04.groupby("Site", group_keys=False)
 .apply(add_row)
 .sort_values(["Site", "Age [ka BP]"])
 .interpolate()
)
chadwick_interp_lr04["Source"] = "Chadwick et al., 2021 (LR04)"

chadwick_interp_edc3 = (subset_edc3.groupby("Site", group_keys=False)
 .apply(add_row)
 .sort_values(["Site", "Age [ka BP]"])
 .interpolate()
)
chadwick_interp_edc3["Source"] = "Chadwick et al., 2021 (EDC3)"
# add back in lonlats
chadwick_127k = pd.concat([chadwick_interp_edc3.loc[chadwick_interp_edc3["Age [ka BP]"] == 127.0],
                         chadwick_interp_lr04.loc[chadwick_interp_lr04["Age [ka BP]"] == 127.0]])

siconc = (chadwick_127k.drop(columns=["Age [ka BP]"])
          .rename(columns={"Sea ice con (9) [%] (Modern analog technique (MAT))": "siconc%"}))
siconc["siconc"] = siconc["siconc%"]/100
#%%
siconc["HadISST"] = siconc.apply(lambda row: get_sic(hadisst_mean, row["Longitude"],row["Latitude"]), axis=1)
siconc["siconc anom"] = siconc["siconc"] - siconc["HadISST"]
# %% load uvic data

path_to_uvic = "/Users/mackenfr/proxy_compare_2025/processed_uvic/"
path_to_ref_data = "/Users/mackenfr/uvic_outputs/"

ref_data = xr.open_dataset(path_to_ref_data + "PRD-PI/280/tavg.nc" , decode_times=False)
mask = xr.open_dataset(path_to_ref_data+"masks.nc", decode_times=False, engine = "netcdf4")


icesheets = ["prdpi","prdlig", "amnlig", "eqblig"]
ref_is = []
uvic_siconc = []
for ice in icesheets:
    filename = ice + "sep_siconc.nc"
    data = xr.open_dataset(path_to_uvic + filename, decode_times=True)
    uvic_siconc.append(data)
    ref_is.append(ice.upper())
# %%
def get_si_uvic(data, longitude, latitude):
    #extract array of lats and lons
    lats=ref_data["latitude"][:]
    lons=ref_data["longitude"][:]   
    lons_adj = np.where(lons>180, lons-360, lons)
    #find values nearest to given coords
    lat_index=(np.abs(lats - latitude)).argmin()
    lon_index=(np.abs(lons_adj - longitude)).argmin()
    siconc = data[lat_index, lon_index]
    return(float(siconc))

for i, ice in enumerate(ref_is):
    siconc[ice] = siconc.apply(lambda row: get_si_uvic(uvic_siconc[i]["O_icefra"], row["Longitude"],row["Latitude"]), axis=1)
    if i>0:
        siconc[ice +" anom"] = siconc[ice]-siconc["PRDPI"]

# %%
norm = plt.Normalize(-1, 1)
lons = ref_data["longitude"]
lats = ref_data["latitude"]
control = uvic_siconc[0]["O_icefra"]
fig, axs = plt.subplots(1, 3, figsize=(10,4), subplot_kw={'projection':ccrs.SouthPolarStereo()})
for i, ax in enumerate(axs):
    data = uvic_siconc[i+1]["O_icefra"] 
    masked_data = xr.where((mask["G_mskt"][i+1]==1), data, np.nan) 
    anom = masked_data - control
    ax.set_extent([-180, 180,-90, -50], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    sm = plt.cm.ScalarMappable(cmap="PuOr", norm=norm)
    pcm = ax.pcolormesh(lons, lats, anom, cmap = "PuOr",transform= ccrs.PlateCarree(),norm=norm)
    ax.contour(lons, lats, masked_data, c="black", levels = [0.15], transform= ccrs.PlateCarree())
    ax.contour(lons, lats, control, c="black", levels = [0.15],linestyles="dashed", transform= ccrs.PlateCarree())

    sns.scatterplot(data = siconc, 
                x = siconc["Longitude"],
                y= siconc["Latitude"], 
                hue = siconc["siconc anom"], 
                marker = "^",
                palette="PuOr",
                transform= ccrs.PlateCarree(),
                hue_norm=norm,
                legend=False,
                edgecolors="black",
                ax=ax,
                s=80)


    ax.set_title(ref_is[i+1])

cb = fig.colorbar(axs[1].collections[1], ax=axs, orientation='horizontal', fraction=0.05, pad=0.04)

# %%
#write to file for rmse/tables

outpath = "/Users/mackenfr/proxy_compare_2025/processed_proxy/"
siconc.to_csv(outpath + "siconc.csv")

# %%
