from scipy.io import loadmat # to load bathymetry
import numpy as np
import matplotlib.pyplot as plt
from setdir import *

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

    unit398.plot(x='longitude',y='latitude', color='r', ax=ax1)
    unit409.plot(x='longitude',y='latitude', color='b', ax=ax1, xlabel='Longitude', ylabel='Latitude')
    plt.legend(['unit_398','unit_409'])


    xsize = 5
    ysize = compute_ysize(xsize, lonlim, latlim)

    fig.set_size_inches(xsize, ysize)

    #fig.savefig('output.png')
    save_figure(fig,'map_units_pos')
    
# Plot maps (needs two glider tracks)
def map_tracks(bathylon,bathylat,bathy,unit409,unit398):

    # Choose axis limits
    latlim = [52, 67]
    lonlim = [-66, -45]

    axes = plt.subplots(nrows=1, ncols=1)
    ax1 = plt.subplot(1,1,1)
    ax1.contour(bathylon, bathylat, bathy)
    ax1.set_ylim(latlim)
    ax1.set_xlim(lonlim)
    fig = plt.gcf()

    unit398.plot(x='m_gps_lon',y='m_gps_lat', color='r', ax=ax1)
    unit409.plot(x='m_gps_lon',y='m_gps_lat', color='b', ax=ax1, xlabel='Longitude', ylabel='Latitude')
    plt.legend(['unit_398','unit_409'])


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
    unit409.plot(x='time', y='pressure_dbar', ylabel='pressure(dbar)', title=titlestr)
    plt.legend().remove()
    plt.gca().invert_yaxis()
    fig = plt.gcf()
    fig.set_size_inches(xsize,ysize)

    plt.show()
    fname = titlestr+'_pseries'
    save_figure(fig, fname)


    
# Plots the most recent 'ndays' of data (pressure, temperature, salinity')
def plot_tseries(unit409,ndays,titlestr):
    axes = plt.subplots(nrows=3, ncols=1,figsize=(10,10))

    maxtime = max_time = unit409.time.max()
    data_lastdays = unit409[unit409.time>=(max_time+datetime.timedelta(days=-ndays))].copy()
    
    # Pressure
    ax1 = plt.subplot(3,1,1)
    data_lastdays.plot(y='pressure_dbar',x='time', 
                          ax=ax1, ylabel='Pressure (dbar)',xlabel='',
                          title=titlestr)
    ax1.xaxis.set_ticklabels([])
    ax1.invert_yaxis()
    ax1.get_legend().remove()

    # Temperature
    ax2 = plt.subplot(3,1,2)
    data_lastdays.plot(x='time',y='sci_water_temp',
                  ax=ax2,
                  xlabel='',ylabel='Temperature')
    ax2.xaxis.set_ticklabels([])
    ax2.get_legend().remove()

    
    # Salnity
    # Were some values very small
    #data_lastdays['derived_salinity'] = data_lastdays['derived_salinity'].replace(0,np.nan)
    #    data_lastdays['derived_salinity'].where(data_lastdays['derived_salinity'] < 0)
    ax3 = plt.subplot(3,1,3)
    data_lastdays.plot(x='time',y='derived_salinity',
                  ax=ax3,
                  ylabel='Salinity')
    ax3.get_legend().remove()
    
    # Save
    fig = plt.gcf()
    fname = titlestr+'_tseries'
    save_figure(fig, fname)
    
def plot_profiles(unit409,ndays,titlestr):
    # Most recent profiles
    maxtime = max_time = unit409.time.max()
    data_lastdays = unit409[unit409.time>=(max_time+datetime.timedelta(days=-ndays))].copy()
    # Maximum pressure for plot limits
    maxp = data_lastdays.pressure_dbar.max()
                  
    # Profile plot of recent data
    axes = plt.subplots(nrows=1, ncols=2,figsize=(10,10))
    ax1 = plt.subplot(1,2,1)
    data_lastdays.plot(y='pressure_dbar',x='derived_salinity', 
                          ax=ax1,
                         ylabel='Pressure (dbar)',xlabel='Salinity',
                         title=titlestr)
    ax1.legend().remove()
    ax1.invert_yaxis()
    ax1.set_ylim([maxp,0])
    forceAspect(ax1,ratio=2)
    ax2 = plt.subplot(1,2,2)
    data_lastdays.plot(y='pressure_dbar',x='sci_water_temp',
                         ax=ax2,
                         xlabel='Temperature',ylabel='')
    ax2.legend().remove()

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

