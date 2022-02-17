"""
Calculate o2 concentration from a glider unit409
"""
from scipy.io import loadmat # to load bathymetry
import numpy as np
from setdir import *
import xarray as xr
import gsw
from seawater import eos80 as sw

#function [o2_mmolL,o2_mmolkg,o2_mlL,o2_solmmolkg] = #geto2(o2phase,T,PSS,P,lon,lat,varargin)
# O2 concentration function according to procedure detailed under:
# http://www.aanderaa.com/productsdetail.php?Oxygen-Optodes-2
#
# Procedure for calculating the oxygen concentration externally according
# to procedure in spreadsheet (TD-280-Oxygen-Optode-Calculation.xls) given
# the foil 7 part foil coefficient(SVUFoilCoef). Calculations are based on
# the modified Stern-Volmer formula proposed by Uchida et al. (2008, J.
# Atm. Oceanic Tech.)
#
# Hardcoded foil ceofficients for VITALS2016, TRINITY2018 Missions
# Foils are SN-124, SN-333
#
# Returns [mmol/L, mmol/kg, ml/L and solubility mmol/kg] 
#
# Translated from geto2.m, originally by
# Nicolai Bronikowski
# Memorial University of Newfoundland
# nbronikowski@mun.ca

def geto2(o2phase,T,PSS,P,lon,lat,snstr):
    SVUFoilCoef = []; # Oxygen Optode Stern Volmer Coefficients
    Mission = []; # Mission
    SN = [];  # Oxygen Sensor Serial Number
    sal_set = 0;

    # Serial number 232 (on glider 398 for TERIFIC2)
    SVUFoilCoef232 = [ 2.72495E-03, 
                    1.14452E-04, 
                    2.32433E-06, 
                    2.17482E+02,  
                   -2.94188E-01, 
                   -5.43610E+01,
                    4.27084E+00]
    
    # Serial number 268 (on glider 409 for TERIFIC2)
    SVUFoilCoef268 = [ 2.71282E-03, 
                    1.14279E-04, 
                    2.26400E-06, 
                    1.59244E+02,  
                   -2.27435E-01, 
                   -3.99901E+01,
                    3.13635E+00]
    
    o2coeff = {
      "232": SVUFoilCoef232,
      "268": SVUFoilCoef268,
    }

    

    ##Calculation routine:
    PT0 = sw.ptmp(PSS,T,P,0)
    C1 = o2coeff[snstr]
    Pr = o2phase
    
    ##Response time correction (Bittig)
    
    ## Apply Stern Volmer Eq. Uchida (2008) to calculate molar oxy
    Ksv = C1[0]+T*C1[1]+np.square(T)*C1[2]
    Po = T*C1[4]+C1[3]
    Pc = (Pr*C1[6])+C1[5]
    
    o2 = np.divide(np.divide(Po,Pc)-1,Ksv)
    
    ## Next compute pres. and salt compensation (Uchida)
    B0 = -6.24097e-3
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    C0 = -3.11680e-7
    
    coeff_p = 0.032 # Updated Uchida 2008, Aanderaa used to be 0.04
    
    Ts = np.log(np.divide(298.15 - T,273.15+T))
    
    comp_factor = np.exp((PSS-sal_set)*(B0+(Ts*B1)+(np.square(Ts)*B2)+
                                       ((Ts*Ts*Ts)*B3)) +
                       (np.square(PSS)-np.square(sal_set))*C0)
    
    press_comp = ((abs(P)/1000)*coeff_p)+1
    
    ## 5) Compute molar, molal and ml/L oxy conc
    dens = sw.dens(PSS,T,P)
    o2_mmolL = o2*comp_factor*press_comp
    o2_mmolkg = (o2_mmolL/(dens/1000))
    o2_mlL = (dens-1000-T + 1000)*o2_mmolkg/44660
    o2_solmmolkg = gsw.O2sol_SP_pt(PSS,PT0)
                       
    

    return o2_mmolL



def calc_o2conc_cal(unit409):
    temp = unit409.sci_water_temp
    salinity = unit409.derived_salinity
    pres = unit409.pressure_dbar
    o2calp = unit409.sci_oxy4_calphase
    o2conc = unit409.sci_oxy4_oxygen
    Sset = 0
    OxyLim = [290, 315]
    id = ~np.isnan(o2conc)

    N = len(unit409.time)
    longi = np.nanmean(unit409.m_gps_lon)*np.ones(N)
    lati = np.nanmean(unit409.m_gps_lat)*np.ones(N)
    snstr = unit409.attrs["optode S/N"]
    
    o2conc_cal = geto2(o2calp,temp,salinity,pres,
                       longi,lati,snstr)
    
    unit409['o2conc_cal'] = o2conc_cal
    
    return unit409