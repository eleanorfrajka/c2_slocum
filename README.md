# c2_slocum
Slocum glider tools working with C2

Objective of code:
- Download near real-time science data and positions from the C2/web services at NOC for a slocum glider
- Parse the data into a netcdf for each glider
- Run some basic mission diagnostics including, e.g.
    - Calculate horizonta resoltuion (spacing between dive/climb profiles)
    - Identify any anomalous signals of e.g. freshwater (useful for TERIFIC)
    - Plot recent positions so waypoints can be picked (manually)
