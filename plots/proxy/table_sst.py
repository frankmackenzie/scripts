#%%
import pandas as pd
import numpy as np
path = "/Users/mackenfr/proxy_compare_2025/processed_proxy/"
annual = pd.read_csv(path + "annual_allproxy.csv", index_col=0)
summer = pd.read_csv(path + "summer_allproxy.csv", index_col=0)


# %%
models = ["PRDLIG anom", 
          "AMNLIG anom",
          "EQBLIG anom",
          "MMM"]




annual["combined error"] = np.where(pd.isnull(annual["SST anom error"]), annual["SST error"], annual["SST anom error"])
annual["MMM"] = annual.iloc[:, 17:26].mean(axis=1)
annual_SO = annual.loc[annual["Latitude"]<0, ["Source", "SST anom", "combined error"] + models]

#%%
summer["combined error"] = np.where(pd.isnull(summer["SST anom error"]), summer["SST error"], summer["SST anom error"])
summer["MMM"] = summer.iloc[:, 17:26].mean(axis=1)
summer_SO = summer.loc[summer["Latitude"]<0, ["Source", "SST anom", "combined error"] + models]

#%%
samesign_ann = []
samesign_summ = []
within_err_ann = []
within_err_summ = []


n_ann = len(annual_SO) #includes site duplicates
n_summ = len(summer_SO)
n_ann_err = len(annual_SO.loc[~pd.isnull(annual_SO["combined error"])])
n_summ_err = len(summer_SO.loc[~pd.isnull(summer_SO["combined error"])])
for model in models:
    n_samesign_annual =  len(annual_SO.loc[((annual_SO["SST anom"]>0)&(annual_SO[model]>0))|((annual_SO["SST anom"]<0)&(annual_SO[model]<0))])
    n_samesign_summer =  len(summer_SO.loc[((summer_SO["SST anom"]>0)&(summer_SO[model]>0))|((summer_SO["SST anom"]<0)&(summer_SO[model]<0))])
    n_withinerr_ann = len(annual_SO.loc[((annual_SO["SST anom"]+annual_SO["combined error"])>annual_SO[model])
                                    &((annual_SO["SST anom"]-annual_SO["combined error"])<annual_SO[model])])
    n_withinerr_summ = len(summer_SO.loc[((summer_SO["SST anom"]+summer_SO["combined error"])>summer_SO[model])
                                    &((summer_SO["SST anom"]-summer_SO["combined error"])<summer_SO[model])])
    samesign_ann.append(n_samesign_annual)
    samesign_summ.append(n_samesign_summer)
    within_err_ann.append(n_withinerr_ann)
    within_err_summ.append(n_withinerr_summ)


counts = pd.DataFrame({"Model": models,
                         "same sign Ann": samesign_ann,
                         "within error Ann": within_err_ann,
                         "same sign Summ": samesign_summ,
                         "within error Summ": within_err_summ},
                         )

percentages = pd.DataFrame({"Model": models,
                         "same sign Ann n = " +str(n_ann): (counts["same sign Ann"]/n_ann)*100,
                         "within error Ann n = " +str(n_ann_err):  (counts["within error Ann"]/n_ann_err)*100,
                         "same sign Summ n = " +str(n_summ):  (counts["same sign Summ"]/n_summ)*100,
                         "within error Summ n = "+str(n_summ_err): (counts["within error Summ"]/n_summ_err)*100},
                         )
# %%
