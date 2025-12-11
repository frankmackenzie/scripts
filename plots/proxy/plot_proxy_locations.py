#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import seaborn as sns
import scipy
#%%
path = "/Users/mackenfr/proxy_compare_2025/processed_proxy/"
annual = pd.read_csv(path + "annual_allproxy.csv", index_col=0)
annual["Type"] = "Annual SST"
summer = pd.read_csv(path + "summer_allproxy.csv", index_col=0)
summer["Type"] = "Summer SST"
siconc = pd.read_csv(path + "siconc.csv", index_col=0)
siconc["Type"] = "Sea ice conc."
sat = pd.read_csv(path + "sat.csv", index_col=0)
sat["Type"] = "SAT"
sat["Source"] = "Capron et al., 2017"
all = pd.concat([annual, summer, siconc, sat])
all_shorter = (all[["Latitude", "Longitude", "Source", "Type"]]
               .loc[all["Source"]!="Chadwick et al., 2021 (EDC3)"])
all_shorter['Source'] = all_shorter['Source'].str.replace('Chadwick et al., 2021 (LR04)', 'Chadwick et al., 2021')

# %% plot locations
fig, ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree()}, figsize= (9, 10))
ax.set_extent([-180, 180,-90, 90], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
sns.scatterplot(data = all_shorter, 
                x = all_shorter["Longitude"],
                y= all_shorter["Latitude"], 
                hue = all_shorter["Type"], 
                style=all_shorter["Source"], transform = ccrs.PlateCarree())

sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.35), ncol=2)

# %%
