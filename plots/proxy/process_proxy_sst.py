
#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import seaborn as sns
import scipy
#%% import data
path_to_data = "/Users/mackenfr/proxy_compare_2025/raw/"
hadisst_raw = xr.load_dataset(path_to_data + "HadISST_sst.nc")

chadwick_raw = pd.read_csv(path_to_data + "Chadwick_2021.csv")
capron_raw = pd.read_csv(path_to_data + "capron_2017.csv", header=12)
hoffman_raw = pd.read_csv(path_to_data + "hoffman_2017.csv", header=14)
chandler_raw = pd.read_csv(path_to_data + "chandler_langebroek.csv", header=76)
# %% #get pi temps from hadisst for each site
def get_T(data, longitude, latitude):
    #extract array of lats and lons
    lats=data["latitude"][:]
    lons=data["longitude"][:]   
    #find values nearest to given coords
    lat_index=(np.abs(lats - latitude)).argmin()
    lon_index=(np.abs(lons - longitude)).argmin()
    SST = data[lat_index, lon_index]
    return(float(SST))

hadisst_pi = hadisst_raw["sst"][:360]
monlen = hadisst_pi.time.dt.days_in_month

# need to bin months if using JFM means
bins = [0,3,6,9,12]
bin_labels = ["JFM", "MAJ", "JAS", "OND"]

weights_alt = (monlen.groupby_bins("time.month", bins, labels = bin_labels) / monlen.groupby_bins("time.month", bins, labels = bin_labels).sum())
#not weighting by month len gives difference < 0.1 deg
hadisst_seas_jfm = (hadisst_pi*weights_alt).groupby_bins("time.month", bins, labels = bin_labels).sum(dim="time")

weights = (monlen.groupby("time.season") / monlen.groupby("time.season").sum())
hadisst_seas = (hadisst_pi*weights).groupby("time.season").sum(dim="time")
hadisst_ann = hadisst_pi.mean(dim="time")
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
                            "SST sum [¬∞C] (Modern analog technique (MAT))"]]
                            .rename(columns={"Age [ka BP] (Age model, LR04 Lisiecki & Ra...)": "Age [ka BP]",
                                            "SST sum [¬∞C] (Modern analog technique (MAT))": "SST",
                                            "Event": "Site"})
)

subset_edc3 = (chadwick_raw[["Event",
                             "Latitude",
                             "Longitude",
                             "Age [ka BP] (Age model, EDC3 (EPICA Ice Do...)", 
                             "SST sum [¬∞C] (Modern analog technique (MAT))"
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

chadwick_127k["Season"] = "S"
chadwick_127k_sst = chadwick_127k.drop(columns=["Age [ka BP]"])
#%% capron
#need to select for 127ka PiAn (anomaly from HadISST - JFM)
# hadISST value or absolute temperature reconstruction T not given, but uncertainty on T is?

capron_127k = (capron_raw[["Station",
                    "Latitude",
                    "Longitude",
                    #"Uncertainty on temperature reconstruction [°C]",
                    "Type",
                    "127 ka Median PIAn [°C]",
                    "127 ka 2s PIAn [°C]"]]
                    .loc[(capron_raw["Type"] == "Summer SST") | (capron_raw["Type"] == "Annual SST")]
                    .rename(columns = {"Station": "Site",
                                       "127 ka Median PIAn [°C]": "SST anom",
                                       "127 ka 2s PIAn [°C]": "SST anom error"})
)
capron_127k["Season"] = capron_127k.apply(lambda row: "S" if row["Type"] == "Summer SST" else "A", axis=1)
capron_127k["Source"] = "Capron et al., 2017"
# SST anom is for JFM from HadISST, want DJF. Add JFM val from HadISST to get abs SST (dodgy? idk what this does to the error)
capron_127k["HadISST JFM"] = capron_127k.apply(lambda row: get_T(hadisst_seas_jfm.loc[{"month_bins": "JFM"}], row["Longitude"],row["Latitude"]), axis=1)
capron_127k["SST"] = capron_127k["HadISST JFM"] + capron_127k["SST anom"]
capron_127k = capron_127k.drop(columns=["Type", "SST anom", "HadISST JFM"])
# recalculate SST anom later from DJF values
#%% hoffman
hoffman_127k = (hoffman_raw[["Station",
                    "Latitude",
                    "Longitude",
                    "Interpolated 127 ka value  (°C)",
                    "Type",
                    "127 ka 2σ (°C)"]]
                    .rename(columns = {"Station": "Site",
                                       "Interpolated 127 ka value  (°C)": "SST",
                                       "127 ka 2σ (°C)": "SST error"})
)
hoffman_127k["Season"] = hoffman_127k.apply(lambda row: "S" if row["Type"] == "Summer SST" else "A", axis=1)
hoffman_127k["Source"] = "Hoffman et al. 2017 (Capron 2017)"
hoffman_127k = hoffman_127k.drop(columns=["Type"])

#%% chandler langebroek
# need to selct 128k bc time isn't evenly spaced so can't interpolate (following Gao 2024)
chandler_128k = (chandler_raw[["Site",
                              "Latitude",
                              "Longitude",
                              "Season (A annual, S summer (months JFM))",
                              "SST [¬∞C] (Reconstructed)",
                              "Age [ka BP] (Chronology follows Lisiecki a...)"
                              ]]
                              .loc[chandler_raw["Age [ka BP] (Chronology follows Lisiecki a...)"]==128]
                              .rename(columns={"Season (A annual, S summer (months JFM))": "Season",
                                               "SST [¬∞C] (Reconstructed)": "SST"})
                              .drop(["Age [ka BP] (Chronology follows Lisiecki a...)"], axis=1)
                              .drop_duplicates()
                              )

chandler_128k["Source"] = "Chandler et al., 2021"




#%% put em together

all_data = pd.concat([chadwick_127k_sst, capron_127k, hoffman_127k, chandler_128k])
annual = all_data.loc[all_data["Season"] == "A"]
summer = all_data.loc[all_data["Season"] == "S"]

#%%
#i suspect this isn't the best way to do this and it gives warnings about indexing but idk how else to do it
annual["HadISST"] = annual.apply(lambda row: get_T(hadisst_ann, row["Longitude"],row["Latitude"]), axis=1)
annual["SST anom"] = annual["SST"] - annual["HadISST"]

summer["HadISST"] = summer.apply(lambda row: get_T(hadisst_seas.loc[{"season": "DJF"}], row["Longitude"],row["Latitude"]), axis=1)
summer["SST anom"] = summer["SST"] - summer["HadISST"]


#%%
path_to_ref_data = "/Users/mackenfr/uvic_outputs/"

path_to_uvic = "/Users/mackenfr/proxy_compare_2025/processed_uvic/"
ref_data = xr.open_dataset(path_to_ref_data + "PRD-PI/280/tavg.nc" , decode_times=False)

icesheets = ["prdpi","prdlig", "amnlig", "eqblig"]
ref_is = []
ref_orbs = []
seas = []
#seas_alt = []
ann = []
for ice in icesheets:
    filename_ann = ice+ "annual.nc"
    filename_seas = ice+ "seasonal.nc" 
    data_ann = xr.open_dataset(path_to_uvic +filename_ann, decode_times=True)
    data_seas = xr.open_dataset(path_to_uvic +filename_seas, decode_times=True)
    ann.append(data_ann)
    seas.append(data_seas)
    ref_is.append(ice.upper())


def get_T_uvic(data, longitude, latitude):
    #extract array of lats and lons
    lats=ref_data["latitude"][:]
    lons=ref_data["longitude"][:]   
    lons_adj = np.where(lons>180, lons-360, lons)
    #find values nearest to given coords
    lat_index=(np.abs(lats - latitude)).argmin()
    lon_index=(np.abs(lons_adj - longitude)).argmin()
    SST = data[lat_index, lon_index]
    return(float(SST))

for i, ice in enumerate(ref_is):
    annual[ice] = annual.apply(lambda row: get_T_uvic(ann[i]["O_temp"][0,:,:], row["Longitude"],row["Latitude"]), axis=1)
    summer[ice] = summer.apply(lambda row: get_T_uvic(seas[i]["O_temp"][0,0,:,:], row["Longitude"],row["Latitude"]), axis=1)
    if i>0:
        annual[ice +" anom"] = annual[ice]-annual["PRDPI"]
        summer[ice +" anom"] = summer[ice]-summer["PRDPI"]

#%% now pmip
models = ['CESM2',
 'NorESM1',
 'FGOALS-f3',
 'FGOALS-g3',
 'NESM3',
 'NorESM2',
 'HadGEM3',
 'ACCESS',
 'GISS']
path_to_pmip = "/Users/mackenfr/proxy_compare_2025/pmip_data/anoms_processed/"
pmip_ann = []
pmip_summ = []
for model in models:
    filename_ann = model + "ann_anom.nc"
    filename_seas = model + "summer_anom.nc" 
    data_ann = xr.open_dataset(path_to_pmip +filename_ann, decode_times=True)
    data_seas = xr.open_dataset(path_to_pmip +filename_seas, decode_times=True)
    pmip_ann.append(data_ann)
    pmip_summ.append(data_seas)

pmip_mmm_ann = xr.open_dataset(path_to_pmip +"ensmean_ann.nc")
pmip_mmm_summ = xr.open_dataset(path_to_pmip +"ensmean_summer.nc")

pmip_ann.append(pmip_mmm_ann)
pmip_summ.append(pmip_mmm_summ)
models.append("MMM")

for i, model in enumerate(models):
     annual[model + " anom"] = annual.apply(lambda row: get_T_uvic(pmip_ann[i]["tos"], row["Longitude"],row["Latitude"]), axis=1)
     summer[model + " anom"] = summer.apply(lambda row: get_T_uvic(pmip_summ[i]["tos"], row["Longitude"],row["Latitude"]), axis=1)

#%%

outpath = "/Users/mackenfr/proxy_compare_2025/processed_proxy/"
annual.to_csv(outpath + "annual_allproxy.csv")
summer.to_csv(outpath + "summer_allproxy.csv")



# %%

