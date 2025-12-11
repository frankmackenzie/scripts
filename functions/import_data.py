import xarray as xr
import pandas as pd

def import_data(path_to_data, co2s = ["280", "300", "350", "400", "450"], config_names=["PRD-PI", "AMN-PI"]):
    config_names = config_names # "PRD-LIG","AMN-LIG"
    #
    tavgs = []
    names = []
    co2_all=[]
    ices = []
    orbs = []
    for config in config_names:
        for co2 in co2s:
            filename_tavg = path_to_data + "/" + config + "/" + co2 + "/tavg.nc"
            data = xr.open_dataset(filename_tavg, decode_times=False, engine = "netcdf4")
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