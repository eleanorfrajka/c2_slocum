"""
Script to visualise glider data (run to update basic plots)

* conda activate env_c2slocum
* python visualise_glider.py

Script to grid glider (slocum) data and make various calculations on the profiles.

1. Works locally on the processed netcdf glider time series *o2.nc
2. Grid glider data
3. Calculate MLD

Runs offline, using the netcdf files created in process_glider_tseries.py
"""
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import datetime
from setdir import *
from plotglider_grid import *
from plotglider_tseries import *
from niceplotting import *
from scipy.io import loadmat # to load bathymetry
import xarray as xr
#import seaborn as sns
import glob
from scipy import stats
import cmocean

#----------------------------------------------------------------------
# Load latest full dataset and gridded profiles
#----------------------------------------------------------------------
glider_names = {
    'unit_409': 'Grease',
    'unit_398': 'Churchill',
}
unit_list = [(k) for k in glider_names.keys()]

for uname in unit_list:
    # Glider data
    fname = uname+'*_data_o2.nc'

    # Extract a list with the names of existing raw data files
    existing_files = glob.glob(cat_interim_path(fname))

    # Check whether there are any
    if len(existing_files) > 0:
        # Extract the end date from the filename
        existing_files = sorted(existing_files)
        latest_file = existing_files[-1]
        # Open the dataset
        data_ds = xr.open_dataset(latest_file)
        
        if uname=='unit_409':
            exec('unit409=data_ds.copy()')
        elif uname=='unit_398':
            exec('unit398=data_ds.copy()')
        
    #Glider positions
    fname = uname+'*_position.nc'

    # Extract a list with the names of existing raw data files
    existing_files = glob.glob(cat_raw_path(fname))

    if len(existing_files) >0:
        existing_files = sorted(existing_files)
        latest_file = existing_files[-1]
        
        pos_ds = xr.open_dataset(latest_file)
        
        if uname=='unit_409':
            exec('unit409pos=pos_ds.copy()')
        elif uname=='unit_398':
            exec('unit398pos=pos_ds.copy()')
            
    # Gridded glider data
    fname = uname+'*_bin10m.nc'
    
    # Extract a list with the names of existing interim data files
    existing_files = glob.glob(cat_proc_path(fname))
    
    # Check whether there are any files
    if len(existing_files) > 0:
        # Extract the most recent filename
        existing_files = sorted(existing_files)
        latest_file = existing_files[-1]
        
    grid_ds = xr.open_dataset(latest_file)
    
    if uname=='unit_409':
        exec('grid409=grid_ds.copy()')
    elif uname=='unit_398':
        exec('grid398=grid_ds.copy()')
        
            
# Location for bathymetry file
matlab_file = 'labsea_66.44W_45.68N_5min.mat'
input_bathy_file = cat_proc_path(matlab_file)
mat_data = loadmat(input_bathy_file)
bathy_data = mat_data['bathy']
bathylat = bathy_data['lat'][0][0].flatten()
bathylon = bathy_data['lon'][0][0].flatten()
bathy = bathy_data['depth'][0][0]
#========================================================================



#----------------------------------------------------------------------
# Initialise the figure directory
#----------------------------------------------------------------------
figdir = create_figdir()

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

# Choice of grid interval (pressure in dbar)
dp=10
#========================================================================



#----------------------------------------------------------------------
# Map the bathymetry and glider tracks using position files
#----------------------------------------------------------------------
# Ok, the glider tracks are spotty: Might mean that the m_gps_lon and m_gps_lat 
# aren't the right variables to download!
map_tracks_pos(bathylon,bathylat,bathy,unit409pos,unit398pos)
#========================================================================


#----------------------------------------------------------------------
# Simple colour sections
#----------------------------------------------------------------------
# glider 409:  Pick some standard variables
varlist2 = ['derived_salinity', 'sci_water_temp', 'o2conc_cal',
            'sci_flbbcd_chlor_units', 'sci_flbbcd_cdom_units', 
           'sci_flbbcd_bb_units'
            ]
plot_sxn(grid409, varlist2,'dist_along_track')
plot_sxn(grid409, varlist2,'timevec')


# glider 398: Pick some standard variables
varlist3 = ['derived_salinity', 'sci_water_temp', 'o2conc_cal',
              'sci_bb2flsv9_chl_scaled', 'sci_bb2flsv9_b532_scaled', 
           'sci_bb2flsv9_b700_scaled'
          ]
plot_sxn(grid398, varlist3, 'dist_along_track')
plot_sxn(grid398, varlist3, 'timevec')
#========================================================================


#----------------------------------------------------------------------
# Recent profile data
#----------------------------------------------------------------------
ndays = 15 # number of recent days to plot
varlist = ['derived_salinity', 'sci_water_temp', 'sci_flbbcd_chlor_units',
          ]

varlist2 = ['derived_salinity', 'sci_water_temp', 'sci_bb2flsv9_chl_scaled']

# Makes some simple line plots with color, but needs significant updating
# to have an accurate and appropriate colorbar
plot_gridprof(grid409,ndays,varlist,'grid409', bathylon,bathylat, bathy)
plot_gridprof(grid398,ndays,varlist2,'grid398', bathylon,bathylat, bathy)
#========================================================================


