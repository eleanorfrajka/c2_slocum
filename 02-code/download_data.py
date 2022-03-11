"""
Script to download glider (slocum) data from the C2/API web services at NOC MARS.  

0. Activate the environment in the terminal window: conda activate env_c2slocum
1. Update your token in myToken.txt: https://api.c2.noc.ac.uk/charon/tokens/issue
2. Run update_data.py at the command line, python download_raw_data.py

Creates netcdf files for later processing of glider data

Requires an up-to-date token in the text file: myToken.txt

FEATURE TO ADD:
Change to user-defined inputs in a config file including
- the glider to download
- variables to include (and serial numbers for sensors?)
- mission start date
- more?
"""
# March 2022 - reduced packages (no plotting or parsing of the glider)
import numpy as np
import pandas as pd
import xarray as xr
import os
import glob
import datetime as dt
import requests
import json
from io import StringIO
import ast # To handle the string conversion when loading json filee
from pathlib import Path
# Own packages of code
from setdir import *


#----------------------------------------------------------------------------
## ISSUE: CHANGE TO A CONFIG FILE WITH USER DEFINABLE PARAMETERS
#----------------------------------------------------------------------------
# Slocum gliders: A dictionary with the key as the serial number ('unit_398') 
# and then the plain text name, "Churchill"
glider_names = {
    'unit_398': 'Churchill',
    'unit_409': 'Grease',
}

# Note: Dictionary keys in glider_names MUST match the serial number format 
# used in the c2 API.  

#----------------------------------------------------------------------------
# Specify the startdate for the download
# ISSUE: In a future edit, only the recent data could be downloaded (rather than 
# redownloading the whole dataset each time)
#----------------------------------------------------------------------------
# Choose a start date in YYYY-MM-DD.  
# Earliest valid data for TERIFIC was 2021-12-12, but there were some in-air 
# measurements before
mission_startdate = '2021-12-12'

#----------------------------------------------------------------------------
# CHOOSE VARIABLE NAMES
#----------------------------------------------------------------------------
# Check the Slocum master data list 8.2 for a range of options
# ISSUE: This was a bit of hit-or-miss to find out which variables existed
var_physics = ['sci_water_pressure', 'sci_water_temp', 'sci_water_cond',
            'derived_salinity', 'derived_potential_density', 'derived_potential_temperature',
           ]

var_other = ['m_final_water_vx', 'm_final_water_vy',
             'm_final_water_vx_at_surface',
             'm_final_water_vy_at_surface',
             'm_water_vx', 'm_water_vy',
             'm_gps_lon', 'm_gps_lat',
             'm_lat', 'm_lon',
            ]

var_oxy = ['sci_oxy4_oxygen',
           'sci_oxy4_calphase',
           'sci_oxy4_temp',
          ]

#----------------------------------------------------------------------------
# CHOOSE VARIABLE NAMES - here's the hit-or-miss part.  
# Try everything and see what sticks!
#----------------------------------------------------------------------------
# ISSUE: Need to know which variables actually exist in the dataset, and only
# pick the right ones.  Otherwise, this might slowdown the download.
# Wetlabs on Unit_398: these parameters seem to work (bb2flsv9)
# Wetlabs on unit_409: these parameters seem to work (flbbcd)
# POSSIBLE TO JUST ADD MORE VARIABLES AND THEY WILL BE REMOVED LATER
# BUT PROBABLY SLOWS THINGS DOWN...
var_bio = ['sci_bb2flsv9_b532_scaled', # units ug/l  - blue?? or green
            'sci_bb2flsv9_b700_scaled', # units ug/l - red
            'sci_bb2flsv9_chl_scaled', # units ug/l
            'sci_bb2flsv9_b532_sig', # units ug/l  - blue?? or green
            'sci_bb2flsv9_b700_sig', # units ug/l - red
            'sci_bb2flsv9_chl_sig', # units ug/l
            'sci_bb2flsv9_b532_ref', # units ug/l  - blue?? or green
            'sci_bb2flsv9_b700_ref', # units ug/l - red
            'sci_bb2flsv9_chl_ref', # units ug/l
            'sci_flbbcd_cdom_units', # ppb - 409
            'sci_flbbcd_chlor_units', # ug/l
            'sci_flbbcd_bb_units', # ??? is this blue backscatter?
            'sci_flbbcd_cdom_scaled', # ppb - 409
            'sci_flbbcd_chlor_scaled', # ug/l
            'sci_flbbcd_bb_scaled', # ??? is this blue backscatter?
            'sci_flbbcd_cdom_sig', # ppb - 409
            'sci_flbbcd_chlor_sig', # ug/l
            'sci_flbbcd_bb_sig', # ??? is this blue backscatter?
          ]

# Some details for the attributes in the netcdf file.
platform_type = 'slocum' # Must be in this format to work with the API
project_name = 'TERIFIC'
institution_name = 'National Oceanography Centre, UK'


# ==============================================================================
# Should not need to edit below here
# ==============================================================================

#----------------------------------------------------------------------------
# Check that output directories exist - if not, then exit
# ISSUE: Is there a better way to do this with error handling? [EFW]
#----------------------------------------------------------------------------
outpath = cat_interim_path('')
if not os.path.isdir(outpath):
    print('Directory not found: '+outpath+'    -- EXITING')
    sys.exit(1)
outpath = cat_proc_path('')
if not os.path.isdir(outpath):
    print('Directory not found: '+outpath+'    -- EXITING')
    sys.exit(1)

#----------------------------------------------------------------------------
# Token and API choices/formatting.
#----------------------------------------------------------------------------
# https://api.c2.noc.ac.uk/charon/tokens/issue
#
# Need to copy and paste the token you generated by logging in at the website
# above into the file 02-code/myToken.txt
with open("myToken.txt", "r") as myfile:
    myToken = myfile.read().replace('\n', '')

from requests.structures import CaseInsensitiveDict
headers = CaseInsensitiveDict()
headers["Accept"] = "application/json"
headers["Authorization"] = f'Bearer {myToken}'

# List of glider serial numbers for API
unit_list = [(k) for k in glider_names.keys()]

# URL for the data
api_root = 'https://api.c2.noc.ac.uk/'

# Platform type for API
platform = platform_type

# Format the time string
time_strf = '%Y%m%d'

# Used to chop data before this date
tstart = pd.Timestamp(mission_startdate+'T00')

# Change this to a later value to download only a subset of the data
download_startdate = mission_startdate+'T00%3A00%3A00'

# Date created (for attributes in netcdf file)
date_created = dt.datetime.now().strftime(time_strf)

#----------------------------------------------------------------------------
# API choices specific for data (not positioning)
#----------------------------------------------------------------------------
# Choice of API website for glider sensor data: timeseries/
api_choice = 'timeseries/observations/'

# Specify format for downloaded file:
# Using csv_combined_transposed rather than csv_combined since it seems to 
# help with getting all the data (not just when the wetlabs was on)
format_choice = 'csv_combined_transposed?'

# Format the variable list for the API
var_list = var_physics+var_bio+var_oxy+var_other
var_str = ''
for i in var_list:
    var_str = var_str+'variable='+i+'&'


#--------------------------------------------------------------
# SCIENCE DATA DOWNLOAD: 
# Loop through the individual gliders & download glider data
# - Format the time variable
# - Change units on pressure from bar to dbar
# - Remove negative derived_salinity values
#--------------------------------------------------------------
for uname in unit_list:
    #--------------------------------------------------------------
    # Format the request
    #--------------------------------------------------------------
    if 0:
        # Time limited - can use this to make the dataset smaller when
        # testing changes.  Just use a later value for 'dstart'
        opt_str = f'from={download_startdate}&'\
        f'{var_str}platform_type={platform}'\
        f'&platform_serial{uname}&reverse_order=false&skip_nulls=false'\
        '&cached=false'
        start_yyyymmdd = '_'+str.replace(dstart[0:10],'-','')
        
    # No time limiting - Download everything
    opt_str = f'{var_str}platform_type={platform}'\
    f'&platform_serial={uname}&reverse_order=false&skip_nulls=false'\
    '&cached=false'
    start_yyyymmdd = ''

    # Concatenate request string
    url = api_root+api_choice+format_choice+opt_str

    #--------------------------------------------------------------
    # Request the data - save as text in variable 'resp'
    #--------------------------------------------------------------
    resp = requests.get(url, headers=headers)

    # Check the response code 
    # (200 is good.  If you get something else, token may need refreshing)
    if not resp.status_code==200:
        print(uname+' - [ resp '+str(resp.status_code)+' ] '\
              'Cannot access data - May need to refresh token? or check URL variable')
    else:
        print(uname+' - [ resp '+str(resp.status_code)+' ] '\
              'Good response code - parsing data')

        # Parse the 'resp' string into a dataFrame
        aa = resp.content.decode("utf-8") 
        data_df = pd.read_csv(StringIO(aa)) # Get rid of the leading b
        data_df = data_df.sort_values(['timestamp']) # Sort by time
 
        # Print a little table to the screen
        #    print(data_df.head(3))

        #--------------------------------------------------------------
        # Clean up time format and convert pressure units to dbar
        #--------------------------------------------------------------
        data_df['time'] = data_df.timestamp.apply(lambda x:
                                    dt.datetime.fromtimestamp(x*0.001))
        data_df = data_df.drop(columns='timestamp')
        # Cut data to post deployment
        data_df_2021 = data_df[data_df.time>=tstart].copy()

        # Change pressure from bars to dbars
        data_df_2021['pressure_dbar'] =  data_df_2021.sci_water_pressure * 10

        # Remove negative salinities
        df1 = data_df_2021['derived_salinity']
        df2 = df1.where(df1>0)
        data_df_2021['derived_salinity'] = df2
        
        
        #--------------------------------------------------------------
        # Format the output file name (and path in ../01-data/01-raw/
        #--------------------------------------------------------------
        # Prepare to convert to xarray
        data_df2 = data_df_2021
        data_df2 = data_df2.set_index("time")
        data_df2 = data_df2.drop(columns="sci_water_pressure")

        # Convert to xarray
        ds_2021 = data_df2.to_xarray()

        # Set some attributes
        maxtimestr = pd.to_datetime(
            ds_2021.time.values.max()).strftime(time_strf)

        # Create a dictionary of attributes
        attr_dict = {"Platform": platform,
                     "End Time": maxtimestr,
                     "Project": project_name,
                     "Institution": institution_name,
                     "Date created": date_created, 
                     "Serial number": uname,
                     "Platform name": glider_names[uname],
                }

        ds_2021 = ds_2021.assign_attrs(attr_dict)

        #--------------------------------------------------------------
        # Clean up Xarray dataset:
        # Remove data ararys that have no non-nan values
        #--------------------------------------------------------------
        for varname in ds_2021.keys():
            data_val = ds_2021[varname].values
            num = np.count_nonzero(~np.isnan(data_val))
            if num==0:
                ds_2021 = ds_2021.drop(varname)
                
        # Save a netcdf file
        outfile = uname+'_'+maxtimestr+'_data.nc'
        outfile_with_path = cat_raw_path(outfile)

        ds_2021.to_netcdf(outfile_with_path)

#--------------------------------------------------------------
# GLIDER POSITIONS DOWNLOAD: 
# Loop through the individual gliders & download glider tracks
#--------------------------------------------------------------
# API website to use
api_choice = 'positions/'
# Data format to request
format_choice = 'positions?'

# Loop through the glider list
for uname in unit_list:
    #--------------------------------------------------------------
    # Format the request
    #--------------------------------------------------------------
    if 0:
        # Time limited
        opt_str = f'from={dstart}&platform_type={platform}'\
        f'&platform_serial={uname}&source_type=internal&time_order=descending'
        start_yyyymmdd = '_'+str.replace(dstart[0:10],'-','')
    
    # No time limiting
    opt_str = f'platform_type={platform}'\
    f'&platform_serial={uname}&source_type=internal&time_order=descending'
    start_yyyymmdd = ''

    # Concatenate request string
    url = api_root+api_choice+format_choice+opt_str

    #--------------------------------------------------------------
    # Request the data - save as text in variable 'resp'
    #--------------------------------------------------------------
    resp = requests.get(url, headers=headers)

    # Check the response code 
    # (200 is good.  If you get something else, token may need refreshing)
    if not resp.status_code==200:
        print(uname+' - [ resp '+str(resp.status_code)+' ] '\
              'Cannot access positions - May need to refresh token? Or check the URL variable')
    else:
        print(uname+' - [ resp '+str(resp.status_code)+' ] '\
              'Good response code - parsing positions')

        # Parse the string into a dataFrame
        json_string = resp.content.decode("utf-8") # Get rid of the leading b
        bb = ast.literal_eval(json_string)[0]
        data_df = pd.DataFrame(bb['positions']['internal'])
        data_df.head()

        #--------------------------------------------------------------
        # Format the output file name (and path in ../01-data/01-raw/
        #--------------------------------------------------------------
        # Prepare to convert to xarray
        data_df2 = data_df
        data_df2["time"] = data_df["time"].astype('datetime64').dt.round('1s')
        data_df2["time_received"] = data_df["time"].astype('datetime64').dt.round('1s')
        data_df2 = data_df2.set_index("time")
        data_df2 = data_df2.drop(columns="source")
        ds_pos = data_df2.to_xarray()

        # Latest date in position file.
        maxtimestr = pd.to_datetime(
            ds_pos.time.values.max()).strftime(time_strf)

       # Create a dictionary of attributes
        attr_dict = {"Platform": platform+' glider',
                     "End Time": maxtimestr,
                     "Project": project_name,
                     "Institution": institution_name,
                     "Date created": date_created, 
                     "Serial number": uname,
                      "Platform name": glider_names[uname],
                }


        ds_pos = ds_pos.assign_attrs(attr_dict)

    #--------------------------------------------------------------
    # Format the output file name (and path in ../01-data/01-raw/
    #--------------------------------------------------------------
    outfile = uname+'_'+maxtimestr+'_position.nc'
    outfile_with_path = cat_raw_path(outfile)

    ds_pos.to_netcdf(outfile_with_path, 'w')


    