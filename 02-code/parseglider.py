from scipy.io import loadmat # to load bathymetry
import numpy as np
from setdir import *
import xarray as xr

def dive_index(unit):
    """ Define Profile index. 
    Use initial vertical sampling dp=5-10dbar loaded by Iridium. Check later with higher vertical sampling.
    Parameters
    ----------
    unit : xarray.Dataset for each glider.
    Returns
    -------
    unit : xarray.Dataset
           Profile index in unit['profile_index'] with .0 for dives and .5 for climbs.
           Localise dives indices (idx_d) and climbs (idx_c).
    """

    pres=unit['pressure_dbar'].values.copy()
    unit['profile_index'] = xr.DataArray(np.full(pres.shape[0],np.nan), dims=['time'])
    # remove surface values (above .5 m)
    pres[ np.where(pres<0.5)[0] ] = np.nan
    # localise ~nan pressure
    pres_idx=np.where(~np.isnan(pres))[0]
    
    # 10-20 shallow points have two identical pressure measurements, potentially separated by nans. 
    # Remove latest measurement, otherwise, inflection points are missed.
    idx_rm=np.where(np.diff(pres[pres_idx])==0)[0]
    pres[pres_idx[idx_rm]]=np.nan
    # relocalise ~nan pressure without identical pressure measurements.
    pres_idx=np.where(~np.isnan(pres))[0]
    
    # Initial index based on 1st profile (398:dive or 409:climb).
    init_index=[1 if np.diff(pres[pres_idx[1:3]])>0 else 1.5][0]
    # Was throwing an error for the last valid pressure point on unit409
    for i in range(pres_idx[1],pres.shape[0]-1):
        # Localise inflection point for 3 valid depths: i.e. sign((p_{z+1}-p_{z})*(p_{z}-p_{z-1}))<0
        if ~np.isnan(pres[i]) :
            idxp = np.where(pres_idx==i)[0][0]
            p0 = pres_idx[idxp-1]

            if not idxp+2>len(pres_idx):
                # Typical process
                p1=pres_idx[idxp+1]
                if np.sign((pres[p1]-pres[i])*(pres[i]-pres[p0])) == -1:
                    init_index+=.5
            else:
                if np.sign(pres[i]-pres[p0]) == -1:
                    init_index += .5
            unit['profile_index'][i]=init_index

    idx_c = np.where( unit['profile_index'] % 1 == .5)[0] # climb indices
    idx_d = np.where( unit['profile_index'] % 1 == 0)[0] # dive indices
    
    return unit, idx_d, idx_c


def bin_dp(u1, unitname, dp):
    """ Bin profile data from glider into coarse bins for quick plotting:
    default of dp=10m
    After binning, may want to linearly interpolate over the gaps.
    Parameters
    ----------
    u1: xarray.Dataset for a glider
    dp: size of pressure bin, eg. gridding onto a 10-m bin
    """
    
    # NOTE: Hard-coded the variable names that will be used to create the 
    # x- and y-axis of the gridded data.  Not super happy about hard-coding. EFW
    presname = 'pressure_dbar'
    idxname = 'profile_index'
    
    # Create a temporary time variable
    u1 = u1.assign(timevar=u1.time.astype('float'))

    # Extract an array of pressure and profile_index values
    pres = u1[presname].values
    prof_idx = u1[idxname].values
    
    # Compute the bins, evenly spaced between the surface and 1000m
    # Since this is for glider data, it should be OK to hardcode the maximum pressure
    pmin = 0
    pmax = 1000
    nbins = int(round((pmax-pmin)/dp))
    bins10m = np.linspace(pmin, pmax, nbins+1)
    pres10m = np.linspace(pmin+dp/2, pmax-dp/2, nbins)

    # Get a unique list of all the dive numbers (whole *.0 - dive, and half *.5 - climb)
    divenum_vec = np.unique(u1[idxname].values)
    # Remove the nan values
    divenum_vec = divenum_vec[~np.isnan(divenum_vec)]

    # Updated to operate on all variables in the data array
    varlist = list(u1.keys())
    varlist.remove(idxname)
    varlist.remove(presname)
    
    # Initialise an empty dictionary of variables
    myvars = dict()
    
    counter = 0
    for i in divenum_vec:
        # Subselect the xarray dataset for the dive or climb indices only
        u11 = u1.where(u1[idxname]==i, drop=True)

        # Calculate the median of values within the bin
        ddive = u11.groupby_bins(presname,bins10m, squeeze=False).mean()
        
        for varname in varlist:
            # Assign dimensions
            data1 = ddive[varname].values

            # Reshape
            data1 = data1[:, np.newaxis]
        
            myvars[varname] = (['pressure', 'divenum'], data1)
        
        mycoords = dict(
            divenum = (["divenum"], [i],
                      dict(long_name = "Profile index")),
            pressure = (["pressure"], pres10m,
                       dict(long_name = "Pressure",
                           units = "dbar")),
        )
    
        myattrs = dict(
            unit = unitname,
        )

        blank_new = xr.Dataset(data_vars=myvars, coords=mycoords, attrs=myattrs)

        if counter==0:
            blank_full = blank_new
            counter += 1
        else:
            blank_full = blank_full.combine_first(blank_new)
            counter += 1
            

    # Convert the temporary time variable back into datetime64, then drop it.
    blank_full = blank_full.assign(time=blank_full.timevar.astype('datetime64[ns]'))
    blank_full = blank_full.drop('timevar')
                                
    return blank_full