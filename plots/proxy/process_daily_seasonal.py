#%%
import xarray as xr

path_to_uvic = "/Users/mackenfr/uvic_outputs/daily/timed/"

icesheets = ["prdpi","prdlig", "amnlig", "eqblig"]
ref_is = []
ref_orbs = []
seas = []
#seas_alt = []
#%%
ann = []
for ice in icesheets:
    filename = ice+"_timed.nc"
    data = xr.open_dataset(path_to_uvic +filename, decode_times=True, chunks = 120)
    if "pi" in ice:
        monlen = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]*100
    else: 
        monlen = [34, 30, 32, 29, 29, 27, 28, 29, 29, 32, 32, 34]*100

    new_monlen = data.time.dt.days_in_month.copy()
    new_monlen.values = monlen
    
    weights_seas = (new_monlen.groupby("time.season") / new_monlen.groupby("time.season").sum())
    data_seas = (data*weights_seas).groupby("time.season").sum(dim="time")

    # # if using JFM for summer:
    # weights_seas_alt = (new_monlen.groupby_bins("time.month", bins, labels = bin_labels) / new_monlen.groupby_bins("time.month", bins, labels = bin_labels).sum())
    # data_seas_alt = (data*weights_seas).groupby_bins("time.month", bins, labels = bin_labels).sum(dim="time")

    weights_ann = (new_monlen / new_monlen.sum())
    data_ann = (data*weights_ann).sum(dim="time")

    ann.append(data_ann)
    seas.append(data_seas)
    #seas_alt.append(data_seas_alt)
    ref_is.append(ice.upper())
#%%
outpath = "/Users/mackenfr/proxy_compare_2025/"

for data, ice in zip(ann, icesheets):
    data.to_netcdf(outpath + "processed_uvic/" + ice + "annual.nc")
for data, ice in zip(seas, icesheets):
    data.to_netcdf(outpath + "processed_uvic/" + ice + "seasonal.nc")

#%% september sea ice
simeans = []
for ice in icesheets:
    filename = ice+"_timed.nc"
    data = xr.open_dataset(path_to_uvic +filename, decode_times=True, chunks = 100)
    septembers = data["O_icefra"].loc[{"time": data.time.dt.month == 9}]
    mean = septembers.mean(dim="time")
    simeans.append(mean)

# %%
for data, ice in zip(simeans, icesheets):
    data.to_netcdf(outpath + "processed_uvic/" + ice + "sep_siconc.nc")
# %%
