#%%
import xarray as xr
import cftime
import datetime

#%%
def convert_time(data):
    newdata = data.copy()
    newdata["time"] = data["time"] - data["time"][0] 
    if "pi" in data.attrs["history"]:
        units = 'common_years since 2000-1-1 00:00:00'
    else:
        units = 'common_years since 1999-12-27 00:00:00'#start date needs to reflect pcal 
    time_365 = cftime.num2date(newdata.time, units = units, calendar = '365_day')
    newtime = xr.DataArray(time_365, coords = [time_365], dims = 'time', name = 'data')
    newdata["time"] = newtime
    newdata["time"] = newdata.indexes["time"].to_datetimeindex(time_unit = 's')
    #newdata["time"] = newdata.indexes["time"]
    return(newdata)



path_to_data = "/Users/mackenfr/uvic_outputs/daily/adjusted/"#"/Users/mackenfr/uvic_outputs/daily/raw/" 
prefixes =  ["prdpi", "prdlig", "amnlig", "eqblig"]
datasets = []
for exp in prefixes:
    data = xr.open_dataset(path_to_data + "combined_" +  exp + "_adj.nc", decode_times=False, chunks={"time":100})
    timed = convert_time(data)
    datasets.append(timed)
#%%
#%%

outpath = "/Users/mackenfr/uvic_outputs/daily/timed/"#"/Users/mackenfr/uvic_outputs/daily/timed/"
for i, data in enumerate(datasets):
    data.to_netcdf(outpath + prefixes[i] +"_timed.nc")
# %%
"""
then  cdo setcalendar,365_day
then need to use /Users/mackenfr/PaleoCalAdjust/extract_vars.sh to split each file into individual variables
then run scripts/paleocal_make_info.py to make the info file for paleocal
then when paleocal adjust has run, cdo -delname,rmonlen,rmonmid,rmonbeg,rmonend -merge *_[name]_adj.nc combined_[name].nc
"""