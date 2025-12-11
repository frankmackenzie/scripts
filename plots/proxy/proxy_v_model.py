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
summer = pd.read_csv(path + "summer_allproxy.csv", index_col=0)

#%% plots
fig, ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([-180, 180,-90, 90], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
c = ax.scatter(annual["Longitude"], annual["Latitude"], c = annual["SST"], cmap = "Spectral_r", transform=ccrs.PlateCarree())
cb = fig.colorbar(c, ax=ax, orientation='horizontal', fraction=0.05, pad=0.04)

# %%
fig, ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([-180, 180,-90, 90], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
sns.scatterplot(data = annual, x = annual["Longitude"],y= annual["Latitude"], hue = annual["Source"])
sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.3), ncol=3)
# %%

labels_ann = ['Capron et al., 2017', 'Hoffman et al., 2017', 'Chandler et al., 2021']
refs_ann = ['Capron et al., 2017','Hoffman et al. 2017 (Capron 2017)', 'Chandler et al., 2021']
markers = ["D", "v", "^","X"]
norm = plt.Normalize(-5, 5)

fig, ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([-180, 180,-90, 90], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)

for i, ref in enumerate(refs_ann):
    subset_ann = summer.loc[summer["Source"] == ref]
    sns.scatterplot(data = subset_ann, 
                x = subset_ann["Longitude"],
                y= subset_ann["Latitude"], 
                hue = subset_ann["SST anom"], 
                marker = markers[i],
                palette="RdBu_r",
                label = labels_ann[i],
                transform= ccrs.PlateCarree(),
                hue_norm=norm,
                legend=False,
                s=80)

hand, labl = ax.get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    lablout.append(l)
    handout.append(h)

leg = ax.legend(handout, lablout, loc = "lower center", bbox_to_anchor=(0.5, -0.6))
for i in list(range(len(hand))):
    leg.legend_handles[i].set_facecolor('black')

cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.05, pad=0.04)
ax.set_title("annual SST anom")


#%%
norm = plt.Normalize(-5, 5)
labels_summ = ['Capron et al., 2017', 'Hoffman et al. 2017', 'Chandler et al., 2021', 'Chadwick et al., 2021']
refs_summ = ['Capron et al., 2017', 'Hoffman et al. 2017 (Capron 2017)', 'Chandler et al., 2021', 'Chadwick et al., 2021 (EDC3)']
fig, ax = plt.subplots(subplot_kw={'projection':ccrs.SouthPolarStereo()})
ax.set_extent([-180, 180,-90, -30], ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(linestyle='--', draw_labels=False, color="gray")
sm = plt.cm.ScalarMappable(cmap="RdBu_r",norm=norm)
for i, ref in enumerate(refs_summ):
    subset_summer = summer.loc[summer["Source"] == ref]
    sns.scatterplot(data = subset_summer, 
                x = subset_summer["Longitude"],
                y= subset_summer["Latitude"], 
                hue = subset_summer["SST anom"], 
                marker = markers[i],
                palette="RdBu_r",
                label = labels_summ[i],
                transform= ccrs.PlateCarree(),
                hue_norm=norm,
                legend=False,
                s=100)

hand, labl = ax.get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    lablout.append(l)
    handout.append(h)
leg = ax.legend(handout, lablout, loc = "lower center", bbox_to_anchor=(0.5, -0.6))
for i in list(range(len(hand))):
    leg.legend_handles[i].set_facecolor('black')

cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.05, pad=0.04)
ax.set_title("summer SST anom")
# %% load uvic data
path_to_uvic = "/Users/mackenfr/proxy_compare_2025/processed_uvic/"
path_to_ref_data = "/Users/mackenfr/uvic_outputs/"

ref_data = xr.open_dataset(path_to_ref_data + "PRD-PI/280/tavg.nc" , decode_times=False)
ref_data = xr.open_dataset(path_to_ref_data + "PRD-PI/280/tavg.nc" , decode_times=False)
mask = xr.open_dataset(path_to_ref_data+"masks.nc", decode_times=False, engine = "netcdf4")

masks = [mask["G_mskt"][0]]
for i in range(3):
    masks.append(mask["G_mskt"][i])
#%%
control_ann = xr.open_dataset(path_to_uvic +"prdpiannual.nc", decode_times=True)
control_seas = xr.open_dataset(path_to_uvic +"prdpiseasonal.nc", decode_times=True)

icesheets = ["prdlig", "amnlig", "eqblig"]
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
    ann.append(data_ann["O_temp"][0,:,:]-control_ann["O_temp"][0,:,:])
    seas.append(data_seas["O_temp"][0,0,:,:]-control_seas["O_temp"][0,0,:,:])
    ref_is.append(ice.upper())

path_to_pmip= "/Users/mackenfr/proxy_compare_2025/pmip_data/anoms_processed/"

pmip_mmm_ann = xr.open_dataset(path_to_pmip +"ensmean_ann.nc")
pmip_mmm_summ = xr.open_dataset(path_to_pmip +"ensmean_summer.nc")

ann_subset = ann.copy()
ann_subset.append(pmip_mmm_ann["tos"])

seas_subset = seas.copy()
seas_subset.append(pmip_mmm_summ["tos"])

is_labels = ref_is.copy()
is_labels.append("PMIP4 MMM")
#%%
models = ['CESM2',
          'NorESM1',
          'FGOALS-f3',
          'FGOALS-g3',
          'NESM3',
          'NorESM2',
          'HadGEM3',
          'ACCESS',
          'GISS']

for model in models:
    filename_ann = model + "ann_anom.nc"
    filename_seas = model + "summer_anom.nc" 
    data_ann = xr.open_dataset(path_to_pmip +filename_ann)
    data_seas = xr.open_dataset(path_to_pmip +filename_seas)
    ann.append(data_ann["tos"])
    seas.append(data_seas["tos"])
    ref_is.append(model)




# %% annual
lons = ref_data["longitude"]
lats = ref_data["latitude"]
norm = plt.Normalize(-3, 3)

fig, axes = plt.subplots(4, 3, figsize=(10,12), subplot_kw={'projection':ccrs.SouthPolarStereo()})
for i, ax in enumerate(axes.flat):
    ax.set_extent([-180, 180,-90, -35], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    values = ann[i]
    if i <= 3:
        values = xr.where((masks[0]==1)&(masks[i]==1), values, np.nan) 
    ax.pcolormesh(lons, lats, values, cmap = "coolwarm", norm=norm,transform= ccrs.PlateCarree())

    for j, ref in enumerate(refs_ann):
        subset_ann = annual.loc[annual["Source"] == ref]
        sns.scatterplot(data = subset_ann, 
                    x = subset_ann["Longitude"],
                    y= subset_ann["Latitude"], 
                    hue = subset_ann["SST anom"], 
                    marker = markers[j],
                    palette="coolwarm",
                    label = labels_ann[j],
                    transform= ccrs.PlateCarree(),
                    hue_norm=norm,
                    legend=False,
                    ax=ax,
                    s=80)
    ax.set_title(ref_is[i], size=8)

hand, labl = axes[0,0].get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    lablout.append(l)
    handout.append(h)

leg = fig.legend(handout, lablout, loc = "lower center")
for k in list(range(len(hand))):
    leg.legend_handles[k].set_facecolor('black')

cb = fig.colorbar(axes[0,0].collections[1], ax=axes, orientation='horizontal', fraction=0.02, pad=0.01)
fig.suptitle("Annual SST", y=0.92)
# %% summer anom

fig, axes = plt.subplots(4, 3, figsize=(10,12), subplot_kw={'projection':ccrs.SouthPolarStereo()})
for i, ax in enumerate(axes.flat):
    ax.set_extent([-180, 180,-90, -35], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    values = seas[i]
    if i <= 3:
        values = xr.where((masks[0]==1)&(masks[i]==1), values, np.nan) 
    ax.pcolormesh(lons, lats, values, cmap = "coolwarm", norm=norm,transform= ccrs.PlateCarree())

    for j, ref in enumerate(refs_summ):
        subset_summer = summer.loc[summer["Source"] == ref]
        sns.scatterplot(data = subset_summer, 
                    x = subset_summer["Longitude"],
                    y= subset_summer["Latitude"], 
                    hue = subset_summer["SST anom"], 
                    marker = markers[j],
                    palette="coolwarm",
                    label = labels_summ[j],
                    transform= ccrs.PlateCarree(),
                    hue_norm=norm,
                    legend=False,
                    ax=ax,
                    s=80)
    ax.set_title(ref_is[i], size=8)

hand, labl = axes[0,0].get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    lablout.append(l)
    handout.append(h)

leg = fig.legend(handout, lablout, loc = "lower center")
for k in list(range(len(hand))):
    leg.legend_handles[k].set_facecolor('black')

cb = fig.colorbar(axes[0,0].collections[1], ax=axes, orientation='horizontal', fraction=0.02, pad=0.01)
fig.suptitle("Summer (DJF) SST", y=0.92)

# %% annual sst corr
norm = plt.Normalize(-90, 90)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)

for i in ref_is:
    mask = ~np.isnan(annual[i]) & ~np.isnan(annual["SST"])
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(annual[i][mask],annual["SST"][mask])
    fig, ax = plt.subplots()
    sns.regplot(x=annual[i],
                y=annual["SST"],
                ax=ax,
                ci=None,
                scatter=False,
                label=None)
    for j, name in enumerate(refs_ann):
        subset=annual.loc[annual["Source"]==name]
        sns.scatterplot(data = subset,
                        x=subset[i], 
                        y = subset["SST"], 
                        hue=subset["Latitude"],
                        marker=markers[j],
                        palette="viridis",
                        label = labels_ann[j],
                        hue_norm=norm,
                        legend=False)
    
    cb = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.04)

    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    ax.annotate("r = " + str(round(r_value, 3)), (5,4))
    ax.annotate("p= " + f"{p_value:.3e}", (5,3))
    ax.patch.set_alpha(0)   
    
    hand, labl = ax.get_legend_handles_labels()
    handout=[]
    lablout=[]
    for h,l in zip(hand,labels_ann[:]):
        lablout.append(l)
        handout.append(h)

    leg = ax.legend(handout, lablout, loc = "lower right")
    for k in list(range(len(hand))):
        leg.legend_handles[k].set_facecolor('black')
    ax.set_title("Annual SST")

# %% summer sst

for i in ref_is:
    mask = ~np.isnan(summer[i]) & ~np.isnan(summer["SST"])
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(summer[i][mask],summer["SST"][mask])
    fig, ax = plt.subplots()
    sns.regplot(x=summer[i],
                y=summer["SST"],
                ax=ax,
                ci=None,
                scatter=False,
                label=None)
    for j, name in enumerate(refs_summ):
        subset=summer.loc[summer["Source"]==name]
        sns.scatterplot(data = subset,
                        x=subset[i], 
                        y = subset["SST"], 
                        hue=subset["Latitude"],
                        marker=markers[j],
                        palette="viridis",
                        label = labels_summ[j],
                        hue_norm=norm,
                        legend=False)
    
    cb = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.04)

    ax.set_xlim(-1, 20)
    ax.set_ylim(-1, 20)
    ax.annotate("r = " + str(round(r_value, 3)), (5,4))
    ax.annotate("p= " + f"{p_value:.3e}", (5,3))
    ax.patch.set_alpha(0)   
    
    hand, labl = ax.get_legend_handles_labels()
    handout=[]
    lablout=[]
    for h,l in zip(hand,labels_summ[1:]):
        lablout.append(l)
        handout.append(h)

    leg = ax.legend(handout, lablout, loc = "lower right")
    for k in list(range(len(hand))):
        leg.legend_handles[k].set_facecolor('black')
    ax.set_title("summer SST")

# %% annual sst anom
fig, axes = plt.subplots(4, 3, figsize=(12,14), layout='constrained')
for i, ax in enumerate(axes.flat):
    mask = ~np.isnan(annual[ref_is[i]+" anom"]) & ~np.isnan(annual["SST anom"])
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(annual[ref_is[i]+" anom"][mask],annual["SST anom"][mask])
    sns.regplot(x=annual[ref_is[i]+" anom"],
                y=annual["SST anom"],
                ax=ax,
                ci=None,
                scatter=False,
                label=None)
    for j, name in enumerate(refs_ann):
        subset=annual.loc[annual["Source"]==name]
        sns.scatterplot(data = subset,
                        x=subset[ref_is[i] +" anom"], 
                        y = subset["SST anom"], 
                        hue=subset["Latitude"],
                        marker=markers[j],
                        palette="viridis",
                        label = labels_ann[j],
                        hue_norm=norm,
                        ax=ax,
                        legend=False)
    
    

    ax.set_xlim(-3, 5)
    ax.set_ylim(-3, 5)
    #ax.annotate("r = " + str(round(r_value, 3)), (3,4))
    #ax.annotate("p= " + f"{p_value:.3e}", (3,3))
    ax.patch.set_alpha(0)   
    ax.set(ylabel=None)
    ax.set(ylabel=None)
fig.supylabel('SST anom (proxy)')
fig.suptitle("Annual SST anomaly (LIG-PI)")
hand, labl = ax.get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labels_ann[:]):
    lablout.append(l)
    handout.append(h)

leg = fig.legend(handout, lablout, loc = "upper right", bbox_to_anchor=(0.91, 0.99))
for k in list(range(len(hand))):
    leg.legend_handles[k].set_facecolor('black')

cb = fig.colorbar(sm, ax=axes[:, 2], orientation='vertical', shrink=0.2,label="latitude" )
# %% summer anom
fig, axes = plt.subplots(4, 3, figsize=(12,14), layout='constrained',sharey=True,sharex=True)
for i, ax in enumerate(axes.flat):
    SH_summer = summer.loc[summer["Latitude"]<0]
    mask = ~np.isnan(summer[ref_is[i]+" anom"]) & ~np.isnan(summer["SST anom"])
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(summer[ref_is[i]+" anom"][mask],summer["SST anom"][mask])
    sns.regplot(x=SH_summer[ref_is[i]+" anom"],
                y=SH_summer["SST anom"],
                ax=ax,
                ci=None,
                scatter=False,
                label=None)
    for j, name in enumerate(refs_summ):
        subset=SH_summer.loc[SH_summer["Source"]==name]
        sns.scatterplot(data = subset,
                        x=subset[ref_is[i] +" anom"], 
                        y = subset["SST anom"], 
                        hue=subset["Latitude"],
                        marker=markers[j],
                        palette="viridis",
                        label = labels_summ[j],
                        hue_norm=norm,
                        ax=ax,
                        legend=False)
    
    

    ax.set_xlim(-5, 7)
    ax.set_ylim(-5, 7)
    #ax.summotate("r = " + str(round(r_value, 3)), (3,4))
    #ax.summotate("p= " + f"{p_value:.3e}", (3,3))
    ax.patch.set_alpha(0)   
    ax.set(ylabel=None)
    ax.set(ylabel=None)
fig.supylabel('SST anom (proxy)')
fig.suptitle("DJF SST anomaly (LIG-PI)")

hand, labl = ax.get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labels_summ[:]):
    lablout.append(l)
    handout.append(h)

leg = fig.legend(handout, lablout, loc = "upper right", bbox_to_anchor=(0.91, 0.99))
for k in list(range(len(hand))):
    leg.legend_handles[k].set_facecolor('black')

cb = fig.colorbar(sm, ax=axes[:, 2], orientation='vertical', shrink=0.2,label="latitude" )
# %% 
all_pmip = xr.open_dataset(path_to_pmip +"catted_summ.nc")

# %%
norm = plt.Normalize(-3, 3)

fig, axes = plt.subplots(2,2, figsize=(10,10), subplot_kw={'projection':ccrs.SouthPolarStereo()})
for i, ax in enumerate(axes.flat):
    ax.set_extent([-180, 180,-90, -35], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    values = ann_subset[i]
    if i <= 3:
         values = xr.where((masks[0]==1)&(masks[i]==1), values, np.nan) 
    ax.pcolormesh(lons, lats, values, cmap = "coolwarm", norm=norm,transform= ccrs.PlateCarree())

    for j, ref in enumerate(refs_ann):
        subset_ann = annual.loc[annual["Source"] == ref]
        sns.scatterplot(data = subset_ann, 
                    x = subset_ann["Longitude"],
                    y= subset_ann["Latitude"], 
                    hue = subset_ann["SST anom"], 
                    marker = markers[j],
                    palette="coolwarm",
                    label = labels_ann[j],
                    transform= ccrs.PlateCarree(),
                    hue_norm=norm,
                    legend=False,
                    ax=ax,
                    s=80)
    ax.set_title(is_labels[i], size=12)

hand, labl = axes[0,0].get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    lablout.append(l)
    handout.append(h)

leg = fig.legend(handout, lablout, loc = "lower center")
for k in list(range(len(hand))):
    leg.legend_handles[k].set_facecolor('black')

cb = fig.colorbar(axes[0,0].collections[1], ax=axes, orientation='horizontal', fraction=0.02, pad=0.01)
fig.suptitle("Annual SST", y=0.92)

# %%

fig, axes = plt.subplots(2,2, figsize=(10,10), subplot_kw={'projection':ccrs.SouthPolarStereo()})
for i, ax in enumerate(axes.flat):
    ax.set_extent([-180, 180,-90, -35], ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(linestyle='--', draw_labels=False, color="gray")
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    values = seas_subset[i]
    if i <= 3:
         values = xr.where((masks[0]==1)&(masks[i]==1), values, np.nan) 
    ax.pcolormesh(lons, lats, values, cmap = "coolwarm", norm=norm,transform= ccrs.PlateCarree())

    for j, ref in enumerate(refs_summ):
        subset_summ = summer.loc[summer["Source"] == ref]
        sns.scatterplot(data = subset_ann, 
                    x = subset_summ["Longitude"],
                    y= subset_summ["Latitude"], 
                    hue = subset_summ["SST anom"], 
                    marker = markers[j],
                    palette="coolwarm",
                    label = labels_summ[j],
                    transform= ccrs.PlateCarree(),
                    hue_norm=norm,
                    legend=False,
                    ax=ax,
                    s=80)
    ax.set_title(is_labels[i], size=12)

hand, labl = axes[0,0].get_legend_handles_labels()
handout=[]
lablout=[]
for h,l in zip(hand,labl):
    lablout.append(l)
    handout.append(h)

leg = fig.legend(handout, lablout, loc = "lower center", bbox_to_anchor=(0.5,-0.01))
for k in list(range(len(hand))):
    leg.legend_handles[k].set_facecolor('black')

cb = fig.colorbar(axes[0,0].collections[1], ax=axes, orientation='horizontal', fraction=0.02, pad=0.01)
fig.suptitle("Summer SST", y=0.92)
# %%
