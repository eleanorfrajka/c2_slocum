"""
Set of functions for plotting glider data from time series data - either unit409 or position data.

Does not handle gridded data.
"""
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from setdir import *
from datetime import datetime, timedelta
import cmocean
from niceplotting import *


#------------------------------------------------------------
# PLOT MAPS FROM GLIDER TSERIES
#------------------------------------------------------------
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

    # Set the aspect ratio to Mercator-like
    xsize = 5
    ysize = compute_ysize(xsize, lonlim, latlim)

    fig.set_size_inches(xsize, ysize)

    #fig.savefig('output.png')
    save_figure(fig,'map_units_pos')

#------------------------------------------------------------
# PLOT MAPS FROM GLIDER POSITION
#------------------------------------------------------------
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

    # Set the aspect ratio to Mercator-like
    xsize = 5
    ysize = compute_ysize(xsize, lonlim, latlim)

    fig.set_size_inches(xsize, ysize)

    #fig.savefig('output.png')
    save_figure(fig,'map_unit409')

#------------------------------------------------------------
# PLOT FULL TIME SERIES OF GLIDER PRESSURE
#------------------------------------------------------------
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



#------------------------------------------------------------
# PLOT RECENT 'ndays' of glider T, S and P time series
#------------------------------------------------------------
# Allows a quick check of whether data are being sent
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

#------------------------------------------------------------
# PLOT RECENT 'ndays' of glider T, S and P profiles
#------------------------------------------------------------
# Quick and dirty check of stratification
def plot_profiles(unit409,ndays,titlestr):
    # Variable names (could be passed as a dictionary)
    timename = 'time'
    presname = 'pressure_dbar'
    salname = 'derived_salinity'
    tempname = 'sci_water_temp'
    
    # Most recent profiles
    max_time = unit409.time.max().values
    dt1 = np.timedelta64(-ndays, 'D')
    ds1 = unit409.where(unit409[timename]>=(max_time+dt1), drop=True)
    min_time = ds1.time.min().values
    
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
    timestr1 = pd.to_datetime(min_time).strftime('%b %d')# b for month MMM
    timestr2 = pd.to_datetime(max_time).strftime('%b %d')
    timestr3 = pd.to_datetime(max_time).strftime('%Y')

    ax2.set_title(timestr1+' - '+timestr2+', '+timestr3)
    ax2.set_ylabel('')
    
    ax2.invert_yaxis()
    ax2.set_ylim([maxp,0])
    forceAspect(ax2,ratio=2)

    # Save
    fig = plt.gcf()

    fname = titlestr+'_profiles'
    save_figure(fig, fname)



    
def map_tracks_mld(bathylon,bathylat,bathy,grid409,grid398,ag):

    # Choose axis limits
    latlim = [54, 64]
    lonlim = [-60, -45]

    axes = plt.subplots(nrows=1, ncols=2)
    ax1 = plt.subplot(1,2,1)
    ax1.contour(bathylon, bathylat, bathy, levels=[-4000,-3000,-2000,-1000,0], colors='k', linestyles='solid')
    ax1.set_ylim(latlim)
    ax1.set_xlim(lonlim)
    fig = plt.gcf()

    lonname = 'longitude'
    latname = 'latitude'
    sc=ax1.scatter(grid409.m_lon[0,:].values,grid409.m_lat[0,:].values,c=grid409.MLD.values,s=60,vmin=-800,vmax=0,cmap='viridis_r')
    sc=ax1.scatter(grid398.m_lon[0,:].values,grid398.m_lat[0,:].values,c=grid398.MLD.values,s=60,vmin=-800,vmax=0,cmap='viridis_r')
    sc=ax1.scatter(ag['LONGITUDE'].values,ag['LATITUDE'].values,c=ag['MLD'].values,s=60,vmin=-800,vmax=0,cmap='viridis_r')
    ax1.set_title('Argo last 2 months')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel("Latitude")

    # Set the aspect ratio to Mercator-like
    xsize = 7
    ysize = compute_ysize(xsize, lonlim, latlim)
    fig.set_size_inches(2*xsize, ysize)
    
    #------------------------------------------
    # Create ds1 with ndays most recent profiles
    #------------------------------------------
    max_time = grid409["timevec"].max().values
    dt1 = np.timedelta64(-15, 'D')
    ds1 = grid409.where(grid409["timevec"]>=(max_time+dt1), drop=True).copy()
    ds2 = grid398.where(grid398["timevec"]>=(max_time+dt1), drop=True).copy()
    min_time = max_time+dt1
    # Create a time string
    timestr1 = pd.to_datetime(min_time).strftime('%b %d')# b for month MMM
    timestr2 = pd.to_datetime(max_time).strftime('%b %d')
    timestr3 = pd.to_datetime(max_time).strftime('%Y')
    timestr = timestr1+' - '+timestr2+', '+timestr3
    last_days_ag=np.where(pd.to_datetime(ag['TIME'].values)>pd.to_datetime((max_time+dt1)))[0]
    
    ax2 = plt.subplot(1,2,2)
    ax2.contour(bathylon, bathylat, bathy, levels=[-4000,-3000,-2000,-1000,0], colors='k', linestyles='solid')
    ax2.set_ylim(latlim)
    ax2.set_xlim(lonlim)
    ax2.set_title(timestr)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel("Latitude")
    
    sc=ax2.scatter(ds1.m_lon[0,:].values,ds1.m_lat[0,:].values,c=ds1.MLD.values,s=60,vmin=-800,vmax=0,cmap='viridis_r')
    sc=ax2.scatter(ds2.m_lon[0,:].values,ds2.m_lat[0,:].values,c=ds2.MLD.values,s=60,vmin=-800,vmax=0,cmap='viridis_r')
    sc=ax2.scatter(ag['LONGITUDE'].values[last_days_ag],ag['LATITUDE'].values[last_days_ag],c=ag['MLD'].values[last_days_ag],s=60,vmin=-800,vmax=0,cmap='viridis_r')
    
    plt.tight_layout()
    pltci=fig.colorbar(sc, orientation="horizontal", ax=[ax1,ax2], fraction=.03, pad=0.1, shrink=.8)
    pltci.ax.set_xlabel(r'MLD [m]',fontsize=14,fontweight='bold')
    
    save_figure(fig,'map_tracks_mld')
    
    