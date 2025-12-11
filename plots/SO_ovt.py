import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
from import_data import import_data
#%%
path_to_data = "/Users/mackenfr/uvic_outputs/" 
tavgs, ref_df = import_data(path_to_data, co2s=["280"], config_names=["PRD-PI", "AMN-PI", "TRS-PI", "EQB-PI"])

#%%
def calculate_density(data):
    lon = data["longitude"]
    lat = data["latitude"]
    SP = data["O_sal"] #practical salinity
    pt = data["O_temp"] #potential temp
    p = data["depth"] #sea pressure (approximate depth in m ~= pressure in dbar)
    SA = gsw.conversions.SA_from_SP(SP, p, lon, lat) #convert practical salinity into absolute salinity
    CT = gsw.conversions.CT_from_pt(SA, pt) #potential temperature to conservative temperature
    sigma2 = gsw.density.sigma2(SA, CT) #calculate potential density anomaly (anomaly from surface)
    return(sigma2)
#%%
data = tavgs[0]
lons= data["longitude_V"]
lats = data["latitude_V"] 

velY= data['O_velY'][0,:,:,:] #meridional velocity in m s^-1
dxu = data["G_dxu"] #u grid width in m
dyu = data["G_dyu"]#u grid length  in m
z = data["G_dzt"]

nbins, sigmin, sigstp = 100,30,0.1

# make the density vector
sigma = sigmin + (np.linspace(1, nbins, nbins) - 0.5) * sigstp
#%%
# calculate sigma2
sigma2 = calculate_density(data)[0,:,:,:]
#interpolate density onto velocity grid 
sigma2 = sigma2.interp(latitude = velY["latitude_V"])
sigma2 = sigma2.interp(longitude = velY["longitude_V"])
sigma2 = sigma2.drop_vars(["latitude", "longitude"])

lat=lats#.expand_dims({"depth": z,"longitude_V": lons}, axis=[0,2])
uX=dxu#.expand_dims({"depth": z}, axis=[0])
Z=z#.expand_dims({"latitude_V": lats, "longitude_V": lons}, axis=[1,2])
#%%
#use sigma as mask so vector have the same length

stacked = xr.Dataset({
    "sigma2": sigma2,
    "lat": lat,
    "velY": velY,
    "uX": uX,
    "uZ": Z,
}).stack(points=sigma2.dims)

sigma1d = stacked['sigma2'].values
lat1d = stacked['lat'].values
weights1d = stacked['velY'].values * stacked['uX'].values * stacked['uZ'].values 
#%%
#latitude vector
lat_ax=data["latitude_V_edges"]
# do 2d historgam with volume transports acuumulated
hist2d_tmp=np.histogram2d(sigma1d,lat1d,bins=(sigma,lat_ax),weights=weights1d)
#set zeros  nan
hist2d=np.where(hist2d_tmp[0]==0,np.nan,hist2d_tmp[0])
#hist2d_ma = np.ma.masked_invalid(hist2d)

#generate streamfunction
streamfunction = np.nancumsum(np.flipud(hist2d), axis=0) / 1e6 #integrate from dense to light at each latitude -> Sv
#AABW_export= streamfunction[136:,29:31].max() #get maximum at 36.9:33.3S and denser than 36.8 which is the AABW export Sv
#AABW_overturning= streamfunction[136:,:30].max() #get maximum at 90S:35S and denser than 36.8 which is  AABW in ORAS5 Sv
# return AABW_export,AABW_overturning
#get_sigma_overturning(tavgs[0])
x, y = np.meshgrid(lat_ax, sigma[::-1])
fig, ax= plt.subplots()
ax.set_ylim((30,39))
ax.invert_yaxis()
ax.set_xlim((-90,90))

value = streamfunction
c = ax.pcolormesh(x, y, streamfunction,
                vmin=-30, vmax=30, 
                cmap="Spectral_r", shading="auto")
fig.colorbar(c)



#%% alt attempt
#convert velocity to volume
data["sigma2"]=sigma2
transport = data["O_velY"] * data["G_dxu"] * data["G_dzt"]
#zonalsum = binned_transport.sum(dim="stacked_depth_latitude_V_longitude_V")
stacked = xr.Dataset({
    "sigma2": sigma2,
    "lat": lat,
    "velY": velY,
    "uX": uX,
    "uZ": Z,
}).stack(points=sigma2.dims)

binned = stacked.groupby_bins(
    stacked["sigma2"],
    bins=sigma,
    restore_coord_dims=True
)
#%%
# integrate from bottom
integ_layers_from_bottom = binned.cumsum() 
# the result of the integration over layers is evaluated at the interfaces
# with psi = 0 as the bottom boundary condition for the integration
bottom_condition = xr.zeros_like(integ_layers_from_bottom.isel({"depth": -1}))
psi_raw = xr.concat([integ_layers_from_bottom, bottom_condition], dim="depth")
psi_raw = psi_raw.chunk({"depth": len(psi_raw["depth"])})  # need to rechunk to new size
# convert kg.s-1 to Sv (1e6 m3.s-1)
psi_Sv = psi_raw / 1.0e6

# %%
x, y = np.meshgrid(lat_ax, sigma[::-1])
fig, ax= plt.subplots()
ax.invert_yaxis()
ax.set_xlim((-90,90))
value = streamfunction
c = ax.pcolormesh(x, y, integ_layers_from_bottom["velY"],
                vmin=-30, vmax=30, 
                cmap="Spectral_r", shading="auto")
fig.colorbar(c)

# %%
# %%# %%
