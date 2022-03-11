"""
Calculate stuff from the glider data

- Mixed layer depth 
"""
import numpy as np
from setdir import *
import xarray as xr
import gsw

def calc_MLD( grid, ref_p=15, drho=0.01, bot_val=1):
    """ Calculate MLD of 2D dens/temp, with ref_depth=10m and drho=0.01 km/m^3 by default.
    Parameters
    ----------
    grid : gridded xarray.Dataset for each glider.
    ref_depth : float
                Reference depth [m].
    drho : float
           density threshold.
    bot_val: 0/1
             If no MLD is detected and MLD is below ~1000 m, MLD=deepest point (bot_val=1), otherwise (bot_val=0), MLD=nan.
    Returns
    -------
    grid : gridded xarray.Dataset with grid['MLD'].
    """
    
    pres = grid['pressure'].values.copy()
    dens = grid['derived_potential_density'].values.copy()
    grid['MLD'] = xr.DataArray(np.full(dens.shape[1],np.nan), dims=['divenum'])
    lat = np.nanmean(grid['m_lat'].values)
    depth = gsw.z_from_p(pres,lat)
    ref_z = gsw.z_from_p(ref_p,lat)
    
    for k in range(dens.shape[1]):
        # No MLD without shallowest point.
        i = np.nanargmin(np.abs(depth - ref_z))
        if np.isnan(dens[i,k]):
            grid['MLD'][k] = np.nan
        # No MLD without enough data.
        elif np.sum(~np.isnan(dens[:,k]))<2:
            grid['MLD'][k] = np.nan
        else:
            dens_diff = dens[:,k] - dens[i,k]
            dens_diff[depth > ref_z] = np.nan
            depth_idx = np.nanargmin(abs(dens_diff - drho))
            grid['MLD'][k] = depth[depth_idx]
        
            # MLD=nan if deeper than ~1000m.
            if bot_val==0:
                if depth.shape[0]-1==depth_idx:
                    grid['MLD'][k] = np.nan
            
    return grid