"""
Script to grid glider (slocum) data and make various calculations on the profiles.

To run it:
- In the terminal, run $ conda activate env_c2slocum (if not already activated)
- In the terminal, run $ python grid_glider_data.py

1. Works locally on the processed netcdf glider time series *o2.nc
2. Grid glider data
3. Calculate MLD

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

