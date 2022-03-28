"""
Set of functions for plotting glider data - where the glider data input are 
the gridded files e.g. grid409 and grid398.

NOT timeseries or position data.
"""
from scipy.io import loadmat # to load bathymetry
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from setdir import *
from datetime import datetime, timedelta
import cmocean
import gsw
from niceplotting import *




#--------------------------------------------------------------------------------------------------------------
# Plots a basic section against profile number.
#
# Varlist corresponds to the variable names to be plotted.
# Note that as written, this hard-codes the colour limits on variables.  
# Additionally functionality is needed to allow limits to be passed (for consistent limits between two figures)
# or to allow them to be determined from the data.
#--------------------------------------------------------------------------------------------------------------
def plot_sxn(ds1, varlist, xcoord):
    
    # How many variables to plot
    nn = len(varlist)


    slevels = [32, 34, 34.6, 34.8, 34.85, 34.9, 35]
    smapstr = 'viridis'
    tlevels = [0, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7]

        

    # Make some simple section plots
    fig,  axes = plt.subplots(nrows=nn, figsize=(10,3*nn))
    
    counter = 0
    for varname in varlist:
        data1 = ds1[varname]

        if varname=='derived_salinity':
            levels = slevels
            cmapstr = smapstr
        elif varname=='sci_water_temp':
            levels=[0, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7]
            cmapstr = 'RdYlBu_r'
        elif varname=='sci_oxy4_oxygen':
            levels = [350, 360, 370, 380, 390, 400, 410, 420]
            cmapstr = 'BrBG'
            cmapstr = cmocean.cm.oxy
        elif varname=='o2conc_cal':
            levels = [280, 290, 300, 310, 320, 330, 340]
            cmapstr = 'BrBG'
            cmapstr = cmocean.cm.oxy

        elif (varname=='sci_flbbcd_chlor_units') | (varname=='sci_bb2flsv9_chl_scaled'):
            levels = [0, .05, .1, .15, .2]
            cmapstr='YlGn'
        elif varname=='sci_flbbcd_cdom_units':
            levels = [-.2, 0, .1, .2, .5]
            levels = [x /1000 for x in levels]
            cmapstr = 'YlOrRd'
        elif varname=='sci_flbbcd_bb_units':
            levels = [1, 1.5, 2, 2.5, 3 , 5]
            levels = [x / 10000 for x in levels]
            cmapstr ='Reds'
        elif (varname=='sci_bb2flsv9_b532_scaled'):
            levels = [0, .2, .5, .7, 1]
            levels = [x / 1000 for x in levels]
            cmapstr = 'Blues'

        elif (varname=='sci_bb2flsv9_b700_scaled'):
            levels = [0, .5, 1, 1.5, 2, 2.5, 3]
            levels = [x/1000 for x in levels]
            cmapstr='Reds'
        else:
            levels = []

        if len(levels)>2: 
            data1.plot.pcolormesh(ax=axes[counter], x=xcoord, y='pressure',
                               ylim=[1000, 0], yincrease=False,
                               add_labels=True, levels=levels, cmap=cmapstr)
        else:
            data1.plot.pcolormesh(ax=axes[counter], x=xcoord, y='pressure',
                               ylim=[1000, 0], yincrease=False,
                               add_labels=True)
    
        if (varname=='derived_salinity') | (varname=='sci_water_temp'):
            axes[counter].plot(ds1[xcoord].values,
                               gsw.p_from_z(ds1['MLD'].values,np.nanmean(ds1['m_lat'].values)),'k')
    
        tstr = ds1.attrs['Serial number']+': '+ds1.attrs['Platform name']
        axes[0].set_title(tstr)

        # Trying to fix xaxis labels
        if xcoord=='timevec':
            axes[counter].xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(axes[-1].xaxis.get_major_locator()))

        counter += 1

    xticks1   = axes[-1].get_xticks()
    plt.tight_layout()
    # ISSUE: Need closer spacing between individual subplots
    
    # Save
    fname = ds1.attrs['Serial number']+'_sxn2'
    save_figure(fig, fname)
    
    
def plot_waterfall(ds_grid1, varlist):
    # Expects the gridded data as input
    divenum = ds_grid1.divenum
    
    presname = 'pressure'
    pres1 = ds_grid1[presname].values.copy()
    mp = len(divenum)

    # Choose some colors
    colors = plt.cm.rainbow(np.linspace(0, 1, mp))


    # How many variables to plot
    nn = len(varlist)

    # Record the max-min values in ddiff and the avg value in dmean
    ddiff = dict()
    dmean = dict()

    for dataname in varlist:
        # Choose some parameters for the waterfall plot
        dsal=np.zeros(len(divenum))
        msal=np.zeros(len(divenum))
        for ddo in range(len(divenum)):
            sal1 = ds_grid1[dataname][:,ddo]
            dsal1 = sal1.max()-sal1.min()
            dsal[ddo]=dsal1
            msal[ddo]=sal1.mean()

        # Typical difference between the maximum and minimum salinty in a profile
        dS = np.nanmedian(dsal) 
        dM = np.nanmedian(msal)
        ddiff[dataname] = dS
        dmean[dataname] = dM
    
    # Make some simple section plots
    fig,  axes = plt.subplots(nrows=nn, figsize=(10,3*nn))
    
    counter = 0
    for dataname in varlist:
        if nn==1:
            ax1 = axes
        else:
            ax1 = axes[counter]

        # scale factor on a profile
        s1 = ddiff[dataname]/4

        for ddo in range(len(divenum)):
            sal1 = ds_grid1[dataname][:,ddo] - dmean[dataname]
            ax1.plot(sal1/s1+ddo, pres1, color=colors[ddo])
    
        # Pressure increases with depth
        ax1.invert_yaxis()
        ax1.set_title(dataname)
        counter += 1
        
    ax1.set_xlabel('Profile index')
    
    
    
    
def plot_gridprof(grid409,ndays,varlist,titlestr,
                  bathylon,bathylat, bathy):
    #------------------------------------------
    # Variable names (could be passed as a dictionary)
    #------------------------------------------
    timename = 'time'
    presname = 'pressure'
    lonname = 'm_gps_lon'
    latname = 'm_gps_lat'
    


    #------------------------------------------
    # Get a mean time per profile - call it timevec
    #------------------------------------------
    time1 = grid409[timename].values
    timenan = pd.isnull(time1)
    time2 = time1.astype('float')
    time2[timenan] = np.nan
    time3 = np.nanmean(time2,axis=0)
    timevec = time3.astype('datetime64[ns]')

    grid409["timevec"] = ('divenum', timevec)

    
    #------------------------------------------
    # Create ds1 with ndays most recent profiles
    #------------------------------------------
    max_time = grid409["timevec"].max().values
    dt1 = np.timedelta64(-ndays, 'D')
    ds1 = grid409.where(grid409["timevec"]>=(max_time+dt1), drop=True).copy()
    min_time = max_time+dt1
    # Create a time string
    timestr1 = pd.to_datetime(min_time).strftime('%b %d')# b for month MMM
    timestr2 = pd.to_datetime(max_time).strftime('%b %d')
    timestr3 = pd.to_datetime(max_time).strftime('%Y')
    timestr = timestr1+' - '+timestr2+', '+timestr3

    
    # Maximum pressure for plot limits
    maxp = ds1[presname].max()
    pres1 = ds1[presname].values.copy()

    # Expects the gridded data as input
    divenum = ds1.divenum
    mp = len(divenum)



    #------------------------------------------
    # Choose some colors
    #------------------------------------------
#    colors = plt.cm.rainbow(np.linspace(0, 1, mp))
#    colors = plt.cm.coolwarm(np.linspace(0, 1, mp))
    colors = cmocean.cm.balance(np.linspace(0,1,mp))
#    cmap = plt.get_cmap('coolwarm',mp)
#    cmap = cmocean.tools.get_dict(cmocean.cm.balance, N=9)
    cmap = cmocean.cm.balance._resample(mp)
    # How many variables to plot
    nn = len(varlist)

    #------------------------------------------
    # Set up a wide, multi-panel plot
    #------------------------------------------
    fig = plt.figure(constrained_layout=True, figsize=(3*(nn+1),7))
    # Specify widths for plots (based on 3 profiles + one map)
    ax_widths = [2.5, 2.5, 2.5, 3]
    ax_height = [3.5]
    # Create the subplots
    gs = fig.add_gridspec(ncols=nn+1, nrows=1, width_ratios=ax_widths,
                          height_ratios=ax_height)
    
    #    axes[0,nn] = gs.subplots(sharey='row')
    counter = 0
    for dataname in varlist:
        # Pick the subplot axes
        if counter==0:
            ax1 = fig.add_subplot(gs[0, counter])
            ax0 = ax1
        else:
            ax1 = fig.add_subplot(gs[0, counter], sharey=ax0)
            plt.setp(ax1.get_yticklabels(), visible=False)
            
        #------------------------------------------
        # Plot the individual profiles
        #------------------------------------------
        for ddo in range(len(divenum)):
            sal1 = ds1[dataname][:,ddo]
            ax1.plot(sal1, pres1, color=colors[ddo])

        # Try setting axis limits
        lb,ub = middle_percent(ds1[dataname], 99)
        lb,ub = expandx(lb, ub) # Make the limits slightly wider
#        print('lb is '+str(lb)+', ub is '+str(ub))
        ax1.set_xlim([lb, ub])

        # Pressure increases with depth
        ax1.invert_yaxis()
        ax1.set_xlabel(dataname)
        if counter==0:
            ax1.set_ylabel('Pressure [dbar]')
            ax1.set_title(titlestr)
        elif counter==nn-1:
            ax1.set_title(timestr)
            
        counter += 1
        ax1.set_ylim([maxp,0])

    #------------------------------------------
    # Add a map to the right of the profiles
    #------------------------------------------
    lonvec = ds1[lonname].values
    xmin = lonvec[~np.isnan(lonvec)].min()-1.75
    xmax = lonvec[~np.isnan(lonvec)].max()+2
    latvec = ds1[latname].values
    ymin = latvec[~np.isnan(latvec)].min()-1.25
    ymax = latvec[~np.isnan(latvec)].max()+1.25

    # Use grid_spec with custom widths.  Not sure how to set the aspect ratio to Mercator.
    ax4 = fig.add_subplot(gs[0,counter])
    import matplotlib 
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs = ax4.contour(bathylon, bathylat, bathy, [-3000, -2000, -1000, 0], 
                     colors=[.75, .75, .75], linewidths=3)
    ax4.set_ylim([ymin, ymax])
    ax4.set_xlim([xmin, xmax])
    ax4.clabel(cs, inline=1, fontsize=10)
    ax4.set_title("Profile locations")
    
    # Plot old positions (before ndays) as gray
    ax4.plot(grid409[lonname],grid409[latname],marker='o', 
             markerfacecolor='w', color=[.25, .25, .25])

    # Plot recent positions (in ndays) with colors
    for ddo in range(len(divenum)):
        lon1 = ds1[lonname][:,ddo]
        lat1 = ds1[latname][:,ddo]
        ax4.plot(lon1,lat1,marker='o',markerfacecolor=colors[ddo],color=[.25, .25, .25])
    
    # Add a colorbar for the time of each profile
    timevec = ds1.time.values
    dmin = timevec[~np.isnan(timevec)].min()
    dmax = timevec[~np.isnan(timevec)].max()
    norm = matplotlib.colors.Normalize(vmin=dmin, vmax=dmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    tvec = np.arange(dmin, dmax, timedelta(days=1))

    # Save
    fig = plt.gcf()

    fname = titlestr+'_gridprof'
    save_figure(fig, fname)
    
    
def plot_gridscat(grid409,ndays,varlist,titlestr,
                  bathylon,bathylat, bathy):
    #------------------------------------------
    # Variable names (could be passed as a dictionary)
    #------------------------------------------
    timename = 'time'
    presname = 'pressure'
    lonname = 'm_gps_lon'
    latname = 'm_gps_lat'
    


    #------------------------------------------
    # Get a mean time per profile - call it timevec
    #------------------------------------------
    time1 = grid409[timename].values
    timenan = pd.isnull(time1)
    time2 = time1.astype('float')
    time2[timenan] = np.nan
    time3 = np.nanmean(time2,axis=0)
    timevec = time3.astype('datetime64[ns]')

    grid409["timevec"] = ('divenum', timevec)

    #------------------------------------------
    # Create ds1 with ndays most recent profiles
    #------------------------------------------
    max_time = grid409["timevec"].max().values
    dt1 = np.timedelta64(-ndays, 'D')
    ds1 = grid409.where(grid409["timevec"]>=(max_time+dt1), drop=True).copy()
    min_time = max_time+dt1
    # Create a time string
    timestr1 = pd.to_datetime(min_time).strftime('%b %d')# b for month MMM
    timestr2 = pd.to_datetime(max_time).strftime('%b %d')
    timestr3 = pd.to_datetime(max_time).strftime('%Y')
    timestr = timestr1+' - '+timestr2+', '+timestr3

    
    # Maximum pressure for plot limits
    maxp = ds1[presname].max()
    pres1 = ds1[presname].values.copy()

    # Expects the gridded data as input
    divenum = ds1.divenum
    mp = len(divenum)


    #------------------------------------------
    # Choose some colors
    #------------------------------------------
    cmap = plt.get_cmap('coolwarm')
    indices = np.linspace(0, cmap.N, len(divenum))
    my_colors = [cmap(int(i)) for i in indices]

    # How many variables to plot
    nn = len(varlist)

    #------------------------------------------
    # Set up a wide, multi-panel plot
    #------------------------------------------
    fig = plt.figure(constrained_layout=True, figsize=(3*(nn+1),4))
    
    ax_widths = [2.5, 2.5, 2.5, 3]
    ax_height = [3.5]
    gs = fig.add_gridspec(ncols=nn+1, nrows=1, width_ratios=ax_widths,
                          height_ratios=ax_height)
    
    #    axes[0,nn] = gs.subplots(sharey='row')
    counter = 0
    for dataname in varlist:
        # Pick the subplot axes
        if counter==0:
            ax1 = fig.add_subplot(gs[0, counter])
            ax0 = ax1
        else:
            ax1 = fig.add_subplot(gs[0, counter], sharey=ax0)
            plt.setp(ax1.get_yticklabels(), visible=False)
            
        presmat = np.tile(pres1, (len(divenum),1))
        presmat = presmat.transpose()
        print(presmat.shape)
        sal1 = ds1[dataname]
        print(sal1.shape)
        
        smap = ax1.scatter(sal1,presmat,s=5,c=ds1['time'],
                           edgecolors='none', marker='o', cmap=cmap)
        counter += 1
        
    cb = fig.colorbar(smap, orientation='vertical', use_gridspec=True,
                     ticks = ds1['timevec'].astype(int))
    
    tvec = pd.to_datetime(ds1['timevec'].values)
    cb.ax.set_yticklabels([tvec.strftime('%b %Y') for index in indices])
    print(ds1)


    fig = plt.gcf()

    fname = titlestr+'_gridscat'
    save_figure(fig, fname)