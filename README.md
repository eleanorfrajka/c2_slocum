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

1. Download raw data files for 398 and 409 from c2

    This takes place in the file update_raw_data.py.  The *only* thing this file does is to download the science and position data.  It takes a while to run, so nothing else should be bundled in here.  (I.e., anything that can be calculated offline, should be.)
    
    To run it, in the terminal window, first activate your environment
    ```$ conda activate env_c2slocum```
    
    In a web browser, update your token: https://api.c2.noc.ac.uk/charon/tokens/issue
    
    Copy the new token and paste it into the file 02-code/myToken.txt
    
    Then in the terminal window, run the update script:
    ```$ python update_raw_data.py```
    This may take some 5-10 minutes.
