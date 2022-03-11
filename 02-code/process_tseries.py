"""
Script to process slocum glider downloaded from C2

To run it:
- In the terminal, run $ conda activate env_c2slocum (if not already activated)
- In the terminal, run $ python process_glider_tseries.py

1. Works locally on the downloaded raw netcdf files.
2. Assign profile index, separate dives and climbs
3. Calculate oxygen

FEATURE TO ADD: 
- Will probably need to calculate TEOS-10 variables
- Any handling of wetlabs data?

Runs offline, using the netcdf files created in download_raw_data.py

Next step: process_glider_tseries.py which works on the raw data and calculates 
some things (profile index, oxygen concentration)
"""

import numpy as np
import xarray as xr
import os
import glob
import datetime as dt
# Own packages of code
from setdir import *
from parseglider import *
from calc_oxy import *

# Choice of grid interval (pressure in dbar)
dp=10

## CHANGE TO A CONFIG FILE WITH USER DEFINABLE PARAMETERS
# Slocum gliders: A dictionary with the key as the serial number ('unit_398') 
# and then the plain text name, "Churchill"
glider_names = {
    'unit_398': 'Churchill',
    'unit_409': 'Grease',
}

sensor_sn = {
    'unit_398': {"optode SN": "232"},
    'unit_409': {"optode SN": "268"},
}
# Dictionary keys MUST match the serial number format used in the API.  


# Choose name for new DataArrays to index the profiles:
idxname = 'profile_index'
# Choose name for new DataArray for pressure in dbar:
presname = 'pressure_dbar'

# List of glider serial numbers for API
unit_list = [(k) for k in glider_names.keys()]


#--------------------------------------------------------------
# DATA PROCESSING:
# - Assign a profile index to separate dives and climbs
#
# Save new files to 01-data/01-raw/ as
#   UNIT_YYYYMMDD_data.nc for the full data as a vector
#--------------------------------------------------------------
for uname in unit_list:
    fname = uname+'*_data.nc'
    
    # Extract a list with the names of existing raw data files
    existing_files = glob.glob(cat_raw_path(fname))

    # Check whether there are any files
    if len(existing_files) > 0:
        # Extract the end date from the filename
        existing_files = sorted(existing_files)
        latest_file = existing_files[-1]
        
        # Open the dataset
        data_ds = xr.open_dataset(latest_file)
        
        #--------------------------------------------------------------
        # Assign profile index (separates dives and climbs)
        # This should actually be moved to process_data.py
        #--------------------------------------------------------------
        # where 20.0 means the twentieth dive (downward profile)
        # and 20.5 is the twentieth climb (upward profile)
        data_ds, _, _ = dive_index(data_ds, presname, idxname)
                 
        #--------------------------------------------------------------
        # Calculate oxygen 
        #--------------------------------------------------------------
        sensorsn1 = sensor_sn[uname]
        data_ds = data_ds.assign_attrs(sensorsn1)
        data_ds = calc_o2conc_cal(data_ds)
        fname2 = latest_file[0:-3]+'_o2.nc'
        fname2 = os.path.basename(fname2)
        
        print('Saving to '+cat_interim_path(fname2))
        data_ds.to_netcdf(cat_interim_path(fname2), mode='w')
        data_ds.close()
        
