"""
Set of functions for tweaking plots for niceness.

NOT specific to glider data
"""
import matplotlib.pyplot as plt
import math
import numpy as np
import cmocean
import seaborn as sns
import matplotlib.dates as mdates

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

#------------------------------------------------------------
# ASPECT RATIO OF MAPS
#------------------------------------------------------------
def compute_ysize(xsize,lonlim,latlim):
    # Set a mercator aspect ratio for maps
    # This one is for figure sizes when there's only one axis
    scale_lon = np.cos(np.mean(latlim)*np.pi/180)
    dlon = lonlim[1]-lonlim[0]
    dlat = latlim[1]-latlim[0]
    ysize = xsize/(dlon*scale_lon)*dlat
    return ysize

def forceAspect(ax,ratio=1):
    # Adjust aspect ratios in subplots
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)


#--------------------------------------------------------------
# CONVERT DECIMAL DEGREES TO DECIMAL MINUTES
# For waypoint planning but could also be for axis ticks
#--------------------------------------------------------------
def dec2deg(dec1):
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


#--------------------------------------------------------------
# Middle_percent - to remove outliers from plot limits
#--------------------------------------------------------------
def middle_percent(data1, perc_keep):
    # Works on a vector/dataarray (not a matrix, yet)
    sal_sort = np.sort(data1)
    sal_sort_nonnan = sal_sort[~np.isnan(sal_sort)]
    NN = len(sal_sort_nonnan)

    # Lower and upper limits (in decimal, so that keep 98% means .01 and .99)
    alpha1 = (100-perc_keep)/2/100
    # Index cutoff
    low_index_limit = round(alpha1*NN)
    high_index_limit = round((1-alpha1)*NN)
    # Value limits
    low_limit = sal_sort_nonnan[low_index_limit]
    high_limit = sal_sort_nonnan[high_index_limit]
    
    return low_limit, high_limit


