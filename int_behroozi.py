#!/bin/python

import yt.units as units
import yt.utilities.physical_constants as phys_const
from yt import YTArray,YTQuantity
from numpy import pi
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from math import *
import matplotlib.pyplot as plt
import pickle

#get stellar mass fraction from Behroozi+2010
fname_mass_frac='mass_frac_tot.txt'
Mhost_beh,mass_frac,mass_frac_z5,mass_frac_z1=np.loadtxt(fname_mass_frac,unpack=True)
interp_mass_frac=interp1d(Mhost_beh,mass_frac)

#fname_mass_frac_z1='mass_frac_tot.txt'
#Mhost_beh_z1,mass_frac_z1=np.loadtxt(fname_mass_frac,unpack=True,usecols=(0,2))
interp_mass_frac_z1=interp1d(Mhost_beh,mass_frac_z1)
interp_mass_frac_z5=interp1d(Mhost_beh,mass_frac_z5)

#fname = 'linear_interpolation_SFR_of_zp1_logMh.pkl'
fname = 'logsfr_of_zp1_logMhz0.pkl'
inf = open(fname,'r')
logSFR = pickle.load(inf)

lMhost_Min=11.5
lMhost_Max=14
nhostbins=int((lMhost_Max-lMhost_Min)/.25 +1)
lMhostbins=np.linspace(lMhost_Min,lMhost_Max,num=nhostbins)
dt=1e8


def Stellar_mass(lMhalo,z):
    #gives total stellar mass at z=.1
    if z==.1 :
        frac_mass=(interp_mass_frac(lMhalo))
    elif z==.5:
        frac_mass=(interp_mass_frac_z5(lMhalo))
    elif z==1:
        frac_mass=(interp_mass_frac_z1(lMhalo))
    else :
        print "ERROR can't compute stellar mass frac"
#    print frac_mass,lMhalo,'halo'
    return lMhalo+frac_mass



def z_to_Time(z,H0=YTQuantity(70.2,'km/s/Mpc')):
    #fitting function from http://arxiv.org/pdf/gr-qc/0506079v2.pdf
    if type(H0) != type(YTQuantity(123.)):
        H0 = YTQuantity(H0,'km/s/Mpc')  #assume it's passed in in standard units
    return ((2.0/H0)/(1+ (1+z)**2)).in_units('yr')

def time_to_z(t,H0=YTQuantity(70.2,'km/s/Mpc')):
    from numpy import sqrt
    #inverting the fitting function form http://arxiv.org/pdf/gr-qc/0506079v2.pdf
    if type(H0) != type(YTQuantity(123.,'AU')):
        H0 = YTQuantity(H0,'km/s/Mpc')  #assume it's passed in in standard units
    if type(t) == type(np.array([1,2,3])):
        if (t<20).all():
            t = YTArray(t*1e9,'yr')
        else:
            t = YTArray(t,'yr')
    elif type(t) == type(1.2):
        if t < 20:
            t = YTQuantity(t*1e9,'yr')
        else:
            t = YTQuantity(t,'yr')
    #otherwise, assume it was passed in as a YTArray or YTQuantity with the right units
    #
    # if type(t) == type(YTArray([1,2,3],'yr')):
    #     if (t<20).all():
    #         t *= 1e9    #assume t was passed in in Gyr
    # elif type(t) == type(YTQuantity(123.,'AU')):
    #     if t < 20:  #assume Gyr:
    #         t = YTQuantity(1e9*t,'yr')
    #     else:
    #         t = YTQuantity(t,'yr')
    H0 = H0.in_units('1/yr')
    return (-H0*t + sqrt(2.0*H0*t - H0*H0*t*t))/(H0*t)


def Stellar_mass_formed(lMhalo,zform):

    logSFR_thishost = logSFR(zform+1,lMhalo)
    total_stellar_mass= 10**logSFR_thishost[0]*dt
    #    print total_stellar_mass,'total stellar mass'
#    print 'a lMhalo=', lMhalo,' forms', logSFR_thishost[0], 'solar masses in stars at z=',zform
    return log10(total_stellar_mass)

    




def Total_Stellar_mass(lMhost,zmin):
    #check for negative contributions
    t_min=z_to_Time(8)
    t_max=z_to_Time(zmin)
#    print t_min,t_max,'tmin'
    time_bins=np.linspace(t_min,t_max,t_max/dt)
#    print t_min,t_max,<time_bins
    all_mass=0
    mass_per_host=np.ones(nhostbins)
    for i in time_bins:
        z=time_to_z(i)
        #get the amount of stars that form during that time and then merge
#        print 'looking at z=',z
        contrib= 10**Stellar_mass_formed(lMhost,z)
#            print contrib,z,'mergers for z'
        if contrib >0:
            all_mass=all_mass+contrib
    print all_mass,lMhost,'merger number ,host'
    return all_mass



def Total_mass(zmin):
    mass_density=np.ones(nhostbins)
    mass_gal=np.ones(nhostbins)
    print nhostbins,'nhostbins'
    for j in range(0,nhostbins):
        mass_density[j]=Total_Stellar_mass(lMhostbins[j],zmin)
        mass_gal[j]=Stellar_mass(lMhostbins[j],zmin)
#    plt.semilogy(mass_gal,mass_density)#, label="zform={0}".format(zform)) 
    plt.plot(mass_gal,np.log10(mass_density), label="integral")
    plt.plot(mass_gal,mass_gal, label="y=x")
    plt.grid()
#    plt.semilogy(lMhostbins,mass_density)#, label="zform={0}".format(zform)) 
    plt.legend(loc="best")
    plt.xlabel(r"$M_{gal}$")
    plt.ylabel(r"$inetgrated$ $mass$")
    plt.savefig("integral_behroozi_dt8_z"+str(zmin)+".pdf",dpi=400)
    plt.show()

    return    



