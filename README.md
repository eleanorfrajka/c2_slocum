# c2_slocum
Slocum glider tools working with C2

Objective of code:
- Download near real-time science data and positions from the C2/web services at NOC for a slocum glider
- Parse the data into a netcdf for each glider
- Run some basic mission diagnostics including, e.g.
    - Calculate horizonta resoltuion (spacing between dive/climb profiles)
    - Identify any anomalous signals of e.g. freshwater (useful for TERIFIC)
    - Plot recent positions so waypoints can be picked (manually)


## How to update the glider dataset for TERIFIC

### Update tokens

In a web browser
- Update your token https://api.c2.noc.ac.uk/charon/tokens/issue 

In 02-code/myToken.txt - input your up-to-date token

### Run the basic downloading, processing, gridding and visualisation scripts.

In the terminal
1. conda activate env_c2slocum
2. python download_data.py

    The download step downloads the glider science data and position data from the API, and saves these as netcdf files on your computer in 01-data/01-raw/

3. python process_tseries.py

    The process step processes the glider time series data, assigning a dive index and calculating the first-look calibrated oxygen.

4. python grid_data.py
    
    The gridding step creates a gridded file of all the data with 10m vertical resolution, and then calculates mixed layer depths.

5. python visualise_glider.py

    Creates some status plots to see how the gliders are doing lately.


## There are also various modules

### calc_glider.py

Only calculations MLD so far

### calc_oxy.py

Calculates the oxygen based on the serial number/calibration coefficients of sensors

### glider_varlist.py

Downloads a variable list from the API

### niceplotting.py

Little snippets of code to make nicer plots (not glider-specific)

### plotglider_grid.py

Plotting functions using the gridded data as the main input variables.


### plotglider_tseries.py

Plotting functions using the time series data as the main input variables

### setdir.py

Sets some paths

### waypoint_picker.ipynb

Probably out-of-date now, but was used to pick waypoints from the glider data.

