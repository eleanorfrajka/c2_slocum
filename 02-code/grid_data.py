"""
Script to grid glider (slocum) data and make various calculations on the profiles.

To run it:
- In the terminal, run $ conda activate env_c2slocum (if not already activated)
- In the terminal, run $ python grid_glider_data.py

1. Works locally on the processed netcdf glider time series *o2.nc
2. Grid glider data
3. Calculate MLD
4. Load Argo data

Runs offline, using the netcdf files created in process_glider_tseries.py
"""

import numpy as np
import xarray as xr
import os
import glob
import datetime as dt
from pathlib import Path
# Own packages of code
from setdir import *
from parseglider import *
from calc_glider import *
from scipy.interpolate import interp1d

# Choice of grid interval (pressure in dbar)
dp=10

## CHANGE TO A CONFIG FILE WITH USER DEFINABLE PARAMETERS
# Slocum gliders: A dictionary with the key as the serial number ('unit_398') 
# and then the plain text name, "Churchill"
glider_names = {
    'unit_398': 'Churchill',
    'unit_409': 'Grease',
}

# List of glider serial numbers for API
unit_list = [(k) for k in glider_names.keys()]



#----------------------------------------------------------------------
# GRIDDING DATA AND CALCULATIONS ON THE GRIDDED DATA 
# Grid the data and make some calculations on the gridded data (MLD)
#----------------------------------------------------------------------
for uname in unit_list:
    fname = uname+'*_data_o2.nc'
    
    # Extract a list with the names of existing interim data files
    existing_files = glob.glob(cat_interim_path(fname))
    
    # Check whether there are any files
    if len(existing_files) > 0:
        # Extract the most recent filename
        existing_files = sorted(existing_files)
        latest_file = existing_files[-1]
        
        # Open the dataset
        data_ds = xr.open_dataset(latest_file)
         
        if 0:
            # Check whether a gridded file has already been created
            # Not yet implemented
            proc_files = glob.glob(cat_interim_path(fname))
            if not len(proc_files) > 0:
                print('No processed files for that glider')
   
        #--------------------------------------------------------------
        # Grid data onto a regular pressure grid (intervals given by dp)
        # - Grid data into a 2d matrix against profile index & pressure grid 
        #    NOTE: Gridding is rough and *not* science quality
        #--------------------------------------------------------------
        grid_ds = bin_dp(data_ds, data_ds.attrs['Serial number'], dp)
       
        # EFW: I think closing these helps with file management & permission 
        # denied problems? 
        data_ds.close()


        #------------------------------------------
        # ADD EXTRA COORDINATES (length divenum)
        #------------------------------------------
        # Simplifies plotting later to plot against time or distance
        mtime = grid_ds.time.mean(dim='pressure').values
        mlon = grid_ds.m_lon.mean(dim='pressure').values
        mlat = grid_ds.m_lat.mean(dim='pressure').values

        # Interpolate over lat and long values
        divenum = grid_ds.divenum.values

        # Lon
        idxnan = (~np.isnan(mlon))
        divenum_nonnan = divenum[idxnan]
        mlon_nonnan = mlon[idxnan]
        flon = interp1d(divenum_nonnan, mlon_nonnan,
                        kind='linear', fill_value="extrapolate")
        mlon_full = flon(divenum)

        # Lat
        idxnan = (~np.isnan(mlat))
        divenum_nonnan = divenum[idxnan]
        mlat_nonnan = mlat[idxnan]
        flat = interp1d(divenum_nonnan,mlat_nonnan,
                        kind='linear', fill_value="extrapolate")
        mlat_full = flat(divenum)

        # Calculate distances from the interpolated lat/lon positions
        dist_km = gsw.distance(mlat_full, mlon_full, 0, axis=0)/1000
        dist_km_pad = np.append(0, dist_km)
        # Cumsum is a problem, need to do something about NaN?
        dist_along_track = np.cumsum(dist_km_pad)

        # Create data array versions
        DAT_2 = xr.DataArray(dist_along_track, 
                             coords={"divenum": grid_ds.divenum},
                             attrs=dict(long_name="Distance", units="km"))
        TIME_2 = xr.DataArray(mtime, 
                              coords={"divenum": grid_ds.divenum},
                             attrs=dict(long_name="Date"))
        LAT_2 = xr.DataArray(mlat_full, 
                             coords={"divenum": grid_ds.divenum},
                            attrs=dict(long_name="Latitude"))
        LON_2 = xr.DataArray(mlon_full, 
                             coords={"divenum": grid_ds.divenum},
                            attrs=dict(long_name="Longitude"))



        grid_ds["dist_along_track"] = DAT_2
        grid_ds["timevec"] = TIME_2
        grid_ds["lonvec"] = LON_2
        grid_ds["latvec"] = LAT_2

        # Change the variables to coordinates
        grid_ds = grid_ds.set_coords(['dist_along_track','timevec',
                                      'lonvec','latvec'])

        #-------------------------------------------------
        # Calculate mixed layer depth
        #-------------------------------------------------
        grid_ds = calc_MLD(grid_ds)

        #-------------------------------------------------
        # Save gridded to 01-data/03-processed/*_bin10m.nc
        #-------------------------------------------------
        # Filename as 'unit_409_YYYYMMDD_bin10m.nc'
        uname = data_ds.attrs['Serial number']
        maxtimestr = data_ds.attrs['End Time']
        outfile = uname+'_'+maxtimestr+'_bin10m.nc'
        print('Saving processed to '+cat_proc_path(outfile))
        grid_ds.to_netcdf(cat_proc_path(outfile), mode='w')
        
        # EFW: I think closing these helps with file management & permission 
        # denied problems? 
        grid_ds.close()    
        
        
#----------------------------------------------------------------------
# Load and save Argo data using argopy
# Calculate Argo MLD
#----------------------------------------------------------------------

# Instantiate default argopy data fetcher
from argopy import DataFetcher as ArgoDataFetcher
argo_loader = ArgoDataFetcher()

date_now = datetime.datetime.now()
date_1mth = date_now+datetime.timedelta(days=-2*30)
# Request data for a specific space/time domain
ag_points = argo_loader.region([-66,-44,45,68,0,1000,
                                date_1mth.strftime("%Y-%m-%d"),
                                date_now.strftime("%Y-%m-%d")]).to_xarray() # 45-68N - 66-44W
# Save by profiles
ag = ag_points.argo.point2profile()
# TEOS-10 variables
ag.argo.teos10(['SA', 'CT', 'PV'])

# DATA_MODE ( R for real time data, D for delayed mode data, A for real time adjusted data )
DATA_MODE = ['Real time','Delayed mode','Adjusted']
print( str(ag.N_PROF.values.shape[0]) + ' profiles and ' + str(np.unique(ag.PLATFORM_NUMBER.values).shape[0]) + ' Argo floats found in the last 2 months')
for i in range(len(DATA_MODE)):
    print('Profiles ' + DATA_MODE[i] + ':' + str( (np.where( ag.DATA_MODE.values == DATA_MODE[i][0] )[0]).shape[0]) )

#-------------------------------------------------
# Calculate MLD
#-------------------------------------------------
ag['MLD']=xr.DataArray( np.full(ag.N_PROF.shape[0],np.nan), coords={"N_PROF": ag.N_PROF})
for i in range(ag.N_PROF.shape[0]):
    if np.nanmax(ag.PRES.values[i,:])>10:
        ag['MLD'][i]=MLD_i( gsw.sigma0(ag.SA.values[i,:],ag.CT.values[i,:]), gsw.z_from_p(ag.PRES.values[i,:],ag.LATITUDE.values[i]))


#-------------------------------------------------
# Save Argo in nc file
#-------------------------------------------------
outfile = 'Argo_'+ date_now.strftime("%Y-%m-%d") +'.nc'
if Path(cat_proc_path(outfile)).is_file()==0:
    print('Saving Argo file to '+cat_proc_path(outfile))
    ag.to_netcdf(cat_proc_path(outfile), mode='w')
    ag.close()    
    
    