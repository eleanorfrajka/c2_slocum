from scipy.io import loadmat # to load bathymetry
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from setdir import *
import seaborn as sns

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

# Set a mercator aspect ratio for maps
# This one is for figure sizes when there's only one axis
def compute_ysize(xsize,lonlim,latlim):
    scale_lon = np.cos(np.mean(latlim)*np.pi/180)
    dlon = lonlim[1]-lonlim[0]
    dlat = latlim[1]-latlim[0]
    ysize = xsize/(dlon*scale_lon)*dlat
    return ysize

# Adjust aspect ratios in subplots
def forceAspect(ax,ratio=1):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)


# Plot maps (needs two glider tracks)
def map_tracks_pos(bathylon,bathylat,bathy,unit409,unit398):

    # Choose axis limits
    latlim = [52, 67]
    lonlim = [-66, -45]

    axes = plt.subplots(nrows=1, ncols=1)
    ax1 = plt.subplot(1,1,1)
    ax1.contour(bathylon, bathylat, bathy)
    ax1.set_ylim(latlim)
    ax1.set_xlim(lonlim)
    fig = plt.gcf()

    lonname = 'longitude'
    latname = 'latitude'
    ax1.plot(unit398[lonname], unit398[latname], color='r')
    ax1.plot(unit409[lonname], unit409[latname], color='b')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel("Latitude")
    plt.legend(['unit_398','unit_409'])
#    sns.set(style='whitegrid')

    xsize = 5
    ysize = compute_ysize(xsize, lonlim, latlim)

    fig.set_size_inches(xsize, ysize)

    #fig.savefig('output.png')
    save_figure(fig,'map_units_pos')
    
# Plot maps (needs two glider tracks)
def map_tracks(bathylon, bathylat, bathy, unit409, unit398):

    # Choose axis limits
    latlim = [52, 67]
    lonlim = [-66, -45]

    axes = plt.subplots(nrows=1, ncols=1)
    ax1 = plt.subplot(1,1,1)
    ax1.contour(bathylon, bathylat, bathy)
    ax1.set_ylim(latlim)
    ax1.set_xlim(lonlim)
    fig = plt.gcf()

    lonname = 'm_gps_lon'
    latname = 'm_gps_lat'
    ax1.plot(unit398[lonname], unit398[latname], color='r')
    ax1.plot(unit409[lonname], unit409[latname], color='b')
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    plt.legend(['unit_398','unit_409'])
#    sns.set(style='white')



    xsize = 5
    ysize = compute_ysize(xsize, lonlim, latlim)

    fig.set_size_inches(xsize, ysize)

    #fig.savefig('output.png')
    save_figure(fig,'map_unit409')

# Plot full time series of pressure from a single glider
def plot_pressure(unit409,titlestr):
    # Plot pressure against time
    xsize=10
    ysize=3
    timename = 'time'
    presname = 'pressure_dbar'

    # Plot the time series of pressure
    ax1 = plt.subplot(1,1,1)
    ax1.plot(unit409[timename], unit409[presname])
    ax1.set_ylabel('Pressure [dbar]')
    ax1.set_title(titlestr)
    plt.gca().invert_yaxis()
    fig = plt.gcf()
    fig.set_size_inches(xsize,ysize)

    # Text in the x axis will be displayed in 'YYYY-mm' format.
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
    # Rotates and right-aligns the x labels so they don't crowd each other.
    for label in ax1.get_xticklabels(which='major'):
        label.set(rotation=30, horizontalalignment='right')

    plt.show()
    

    # Save the figure
    fname = titlestr+'_pseries'
    save_figure(fig, fname)


    
# Plots the most recent 'ndays' of data (pressure, temperature, salinity')
def plot_tseries(unit409,ndays,titlestr):
    axes = plt.subplots(nrows=3, ncols=1,figsize=(10,10))

    timename = 'time'
    presname = 'pressure_dbar'
    salname = 'derived_salinity'
    tempname = 'sci_water_temp'

    max_time = unit409.time.max().values
    dt1 = np.timedelta64(-ndays, 'D')
    ds1 = unit409.where(unit409[timename]>=(max_time+dt1))
    
    
    # Pressure
    ax1 = plt.subplot(3,1,1)
    ax1.plot(ds1[timename], ds1[presname])
    ax1.set_ylabel('Pressure [dbar]')
    ax1.set_xlabel('')
    ax1.set_title(titlestr)
    ax1.xaxis.set_ticklabels([])
    ax1.invert_yaxis()

    # Temperature
    ax2 = plt.subplot(3,1,2)
    ax2.plot(ds1[timename], ds1[tempname])
    ax2.set_ylabel('Temperature [deg C]')
    ax2.set_xlabel('')
    ax2.xaxis.set_ticklabels([])

    
    # Salnity
    # Were some values very small
    #data_lastdays['derived_salinity'] = data_lastdays['derived_salinity'].replace(0,np.nan)
    #    data_lastdays['derived_salinity'].where(data_lastdays['derived_salinity'] < 0)
    ax3 = plt.subplot(3,1,3)
    ax3.plot(ds1[timename], ds1[salname])
    ax3.set_ylabel('Salinity')
    
    # Save
    fig = plt.gcf()
    fname = titlestr+'_tseries'
    save_figure(fig, fname)
    
def plot_profiles(unit409,ndays,titlestr):
    # Variable names (could be passed as a dictionary)
    timename = 'time'
    presname = 'pressure_dbar'
    salname = 'derived_salinity'
    tempname = 'sci_water_temp'
    
    # Most recent profiles
    max_time = unit409.time.max().values
    dt1 = np.timedelta64(-ndays, 'D')
    ds1 = unit409.where(unit409[timename]>=(max_time+dt1))

    
    # Maximum pressure for plot limits
    maxp = ds1[presname].max()
                  
    # Profile plot of recent data
    axes = plt.subplots(nrows=1, ncols=2,figsize=(10,10))
    ax1 = plt.subplot(1,2,1)
    ax1.plot(ds1[salname], ds1[presname])
    ax1.set_ylabel('Pressure [dbar]')
    ax1.set_xlabel('Salinity')
    ax1.set_title(titlestr)
    ax1.invert_yaxis()
    ax1.set_ylim([maxp,0])
    forceAspect(ax1,ratio=2)
    
    ax2 = plt.subplot(1,2,2)
    ax2.plot(ds1[tempname], ds1[presname])
    ax2.set_xlabel('Temperature [deg C]')
    ax2.set_ylabel('')
    
    ax2.invert_yaxis()
    ax2.set_ylim([maxp,0])
    forceAspect(ax2,ratio=2)

    # Save
    fig = plt.gcf()

    fname = titlestr+'_profiles'
    save_figure(fig, fname)
    

def dec2deg(dec1):
    import math
    if dec1==0:
        deg1 = 0
        mindec = 0
        degstr = str(deg1)+u"\N{DEGREE SIGN}"
    else:
        if dec1<0:
            dec1 = abs(dec1)
            dirstr = 'W/S'
        elif dec1>0:
            dirstr = 'E/N'

        deg1 = math.floor(dec1)
            
        mindec = round(100*(dec1-deg1)*60)/100
        
        if mindec==60:
            deg1 += deg1
            mindec = 0
            
        if mindec==0:
            degstr = str(deg1)+u"\N{DEGREE SIGN}"+dirstr
        else:
            degstr = str(deg1)+u"\N{DEGREE SIGN}"+str(mindec)+dirstr
            
    return degstr, deg1, mindec


def deg2dec(deg1, mindec):
    dec1 = deg1 + mindec/60
    
    return dec1

