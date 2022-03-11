# c2_slocum
Slocum glider tools working with C2

Objective of code:
- Download near real-time science data and positions from the C2/web services at NOC for a slocum glider
- Parse the data into a netcdf for each glider
- Run some basic mission diagnostics including, e.g.
    - Calculate horizonta resoltuion (spacing between dive/climb profiles)
    - Identify any anomalous signals of e.g. freshwater (useful for TERIFIC)
    - Plot recent positions so waypoints can be picked (manually)


## Getting started 

1. **update_raw_data.py** Downloads raw data files for 398 and 409 from c2

    This takes place in the file update_raw_data.py.  The *only* thing this file does is to download the science and position data.  It takes a while to run, so nothing else should be bundled in here.  (I.e., anything that can be calculated offline, should be.)
    
    To run it, 
    
    - In the terminal window, first activate your environment
    ```$ conda activate env_c2slocum```
    
    - In a web browser, update your token: https://api.c2.noc.ac.uk/charon/tokens/issue
    
    - Copy the new token and paste it into the file 02-code/myToken.txt
    
    - Then in the terminal window, run the update script:
    ```$ python update_raw_data.py```
    This may take some 5-10 minutes.  The netcdf files are then stored in ```01-data/01-raw/``` and have the format ```unit_409_20220311_data.nc``` for the science data, and ```unit_409_20220311_position.nc``` for the position data.

2. **02-process_data.ipynb** Processes the downloaded files by

   - Calculating oxygen concentration (uses calc_o2_conc_cal from calc_oxy.py).  Saves the result into ```01-data/02-intermediate/``` with the phrase ```_o2``` appended to the filename.
   - Grids the data to a 10m vertical grid (this takes a *long* time) using ```bin_dp``` from ```parseglider.py```.  Saves the result into ```01-data/03-processed/``` with the file name like ```unit_409_20220311_bin10m.nc```


3. **04-visualise.ipynb** Creates some figures (Eleanor's version - probably Louis' version should be separate to avoid conflicts!). 

    The scripts to actually make figures are stored in ```plotglider.py```.  This notebook loads the time series data and the gridded data, and makes various maps, timeseries plots, profiles and sections.  Very basic plots so far for quicklooks at realtime data.
    
    
    
Various packages

**parseglider.py** - Takes existing netcdf glider data and assigns dive index, bins data.  This basically reformats existing data without calculating anything new.

**calc_oxy.py** - Functions to calculate the oxygen following a script from Nicolai Bronikowski.  Requires both GSW and seawater toolboxes.

**plotglider.py** - Functions to plot glider data based on the time series or the gridded files.




    
