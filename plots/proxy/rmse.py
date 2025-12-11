#%%
import pandas as pd
import numpy as np
path = "/Users/mackenfr/proxy_compare_2025/processed_proxy/"
annual = pd.read_csv(path + "annual_allproxy.csv", index_col=0)
summer = pd.read_csv(path + "summer_allproxy.csv", index_col=0)
siconc = pd.read_csv(path + "siconc.csv", index_col=0)
sat = pd.read_csv(path + "sat.csv", index_col=0)

#several pmip models and hadisst have missing data as nearest cell to proxy site is masked - could fix w interpolation
annual = annual[['Site', 'Latitude', 'Longitude', 'Source', 'Season',
       'HadISST', 'SST anom', 'PRDPI', 'PRDLIG',
       'PRDLIG anom', 'AMNLIG', 'AMNLIG anom', 'EQBLIG', 'EQBLIG anom',
       'CESM2 anom', 'NorESM1 anom', 'FGOALS-f3 anom', 'FGOALS-g3 anom',
       'NESM3 anom', 'NorESM2 anom', 'HadGEM3 anom', 'ACCESS anom',
       'GISS anom', 'MMM anom']].dropna(how='any',axis=0)
summer = summer[['Site', 'Latitude', 'Longitude', 'Source', 'Season',
       'HadISST', 'SST anom', 'PRDPI', 'PRDLIG',
       'PRDLIG anom', 'AMNLIG', 'AMNLIG anom', 'EQBLIG', 'EQBLIG anom',
       'CESM2 anom', 'NorESM1 anom', 'FGOALS-f3 anom', 'FGOALS-g3 anom',
       'NESM3 anom', 'NorESM2 anom', 'HadGEM3 anom', 'ACCESS anom',
       'GISS anom','MMM anom']].dropna(how='any',axis=0)

# %%

#theres some nans in sstanom bc of hadisst grid - ignore these for now, come back and fix (interp?)
filtered_ann = annual[["AMNLIG","SST anom"]].loc[(annual["Source"]=="Capron et al., 2017")&(annual["Latitude"]<0)]
rmse = np.sqrt(sum((0-filtered_ann["SST anom"])**2)/filtered_ann["SST anom"].count())
# %%
refs = ['Capron et al., 2017', 'Hoffman et al. 2017 (Capron 2017)', 'Chandler et al., 2021', 'Chadwick et al., 2021 (EDC3)']
labels = ['Capron2017', 'Hoffman2017', 'Chandler2021', 'Chadwick2021']

models = ["PRDLIG", 
          "AMNLIG",
          "EQBLIG",
          "CESM2",
          "NorESM1",
          "FGOALS-f3",
          "FGOALS-g3",
          "NESM3",
          "NorESM2",
          "HadGEM3",
          "ACCESS",
          "GISS",
          "MMM",
          "null"]

annual_rmse = pd.DataFrame(index=models, columns=labels)
summer_rmse = pd.DataFrame(index=models, columns=labels)


for i, ref in enumerate(refs):
    rmses_ann = []
    rmses_summer = []
    for model in models:
        subset_ann = annual.loc[(annual["Source"]==ref)]
        subset_summ = summer.loc[(summer["Source"]==ref)]
        if model == "null":
            rmse_ann = np.sqrt(sum((0-subset_ann["SST anom"])**2)/subset_ann["SST anom"].count())
            rmse_summ = np.sqrt(sum((0-subset_summ["SST anom"])**2)/subset_summ["SST anom"].count())

        else:
            rmse_ann = np.sqrt(sum((subset_ann[model+" anom"]-subset_ann["SST anom"])**2)/subset_ann["SST anom"].count())
            rmse_summ = np.sqrt(sum((subset_summ[model+" anom"]-subset_summ["SST anom"])**2)/subset_summ["SST anom"].count())
        rmses_summer.append(rmse_summ)
        rmses_ann.append(rmse_ann)

    annual_rmse[labels[i]]=rmses_ann
    summer_rmse[labels[i]]=rmses_summer
# %%
uvic = ["PRDLIG", 
          "AMNLIG",
          "EQBLIG",
          "null"]

siconc_rmse = pd.DataFrame(index=uvic, columns=["Chadwick 2021 SI conc"])
sat_rmse = pd.DataFrame(index=uvic, columns=["Capron 2017 SAT"])

siconc_rmses = []
sat_rmses = []
for model in uvic:
    if model == "null":
        rmse_siconc = np.sqrt(sum((0-siconc["siconc anom"])**2)/siconc["siconc anom"].count())
        rmse_sat = np.sqrt(sum((0-sat["SAT anom"])**2)/sat["SAT anom"].count())
    else:
        rmse_siconc = np.sqrt(sum((siconc[model+" anom"]-siconc["siconc anom"])**2)/siconc["siconc anom"].count())
        rmse_sat = np.sqrt(sum((sat[model+" anom"]-sat["SAT anom"])**2)/sat["SAT anom"].count())   
    siconc_rmses.append(rmse_siconc)
    sat_rmses.append(rmse_sat)
siconc_rmse["Chadwick 2021 SI conc"]=siconc_rmses
sat_rmse["Capron 2017 SAT"]=sat_rmses
# %%
