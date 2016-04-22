
import yt.units as units
import yt.utilities.physical_constants as phys_const
from yt import YTArray,YTQuantity
from numpy import pi
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from math import *
import matplotlib.pyplot as plt
from cosmological_utils import *
import pickle

#get stellar mass fraction from Behroozi+2010
fname_mass_frac='mass_frac_tot.txt'
Mhost_beh,mass_frac,mass_frac_z5,mass_frac_z1=np.loadtxt(fname_mass_frac,unpack=True)
interp_mass_frac=interp1d(Mhost_beh,mass_frac)
interp_mass_frac_z1=interp1d(Mhost_beh,mass_frac_z1)
interp_mass_frac_z5=interp1d(Mhost_beh,mass_frac_z5)

fname = 'logsfr_of_zp1_logMhz0.pkl'
fname_gal = 'behroozi_am.pkl'

inf = open(fname,'r')
logSFR = pickle.load(inf)
logSFR_gal=pickle.load(open(fname_gal,'r'))

lMhost_Min=11.5
lMhost_Max=14
nhostbins=int((lMhost_Max-lMhost_Min)/.25 +1)
lMhostbins=np.linspace(lMhost_Min,lMhost_Max,num=nhostbins)

lMgal_Min=7
lMgal_Max=11.25
ngalbins=int((lMgal_Max-lMgal_Min)/.25 +1)
lMgalbins=np.linspace(lMgal_Min,lMgal_Max,num=ngalbins)

dt=1e8


def Stellar_mass(lMhalo,z):
    #gives total stellar mass at z
    if z==.1 :
        frac_mass=(interp_mass_frac(lMhalo))
    elif z==.5:
        frac_mass=(interp_mass_frac_z5(lMhalo))
    elif z==1:
        frac_mass=(interp_mass_frac_z1(lMhalo))
    else :
        print "ERROR can't compute stellar mass frac"

    return lMhalo+frac_mass


def Stellar_mass_formed(lMhalo,zform):

    logSFR_thishost = logSFR(zform+1,lMhalo)
    total_stellar_mass= 10**logSFR_thishost[0]*dt
    #    print total_stellar_mass,'total stellar mass'
#    print 'a lMhalo=', lMhalo,' forms', logSFR_thishost[0], 'solar masses in stars at z=',zform
    return log10(total_stellar_mass)


def Stellar_mass_formed_gal(lMgal,zform):

    logSFR_thishost = logSFR_gal((zform+1),lMgal)
    total_stellar_mass= 10**logSFR_thishost[0]*dt
    #    print total_stellar_mass,'total stellar mass'
#    print 'a lMhalo=', lMhalo,' forms', logSFR_thishost[0], 'solar masses in stars at z=',zform
    return log10(total_stellar_mass)



def Total_Stellar_mass(lMhost,zmin):
    t_min=z_to_Time(8.5)
    t_max=z_to_Time(zmin)
    time_bins=np.linspace(t_min,t_max,(t_max-t_min)/dt)
    all_mass=0
    mass_per_host=np.ones(nhostbins)
    for i in time_bins:
        z=time_to_z(i)
        contrib= 10**Stellar_mass_formed(lMhost,z)
        all_mass=all_mass+contrib
    print all_mass,lMhost,'merger number ,host'
    return all_mass


def Total_Stellar_mass_gal(lMgal,zmin):
    t_min=z_to_Time(8.5)
    t_max=z_to_Time(zmin)
    time_bins=np.linspace(t_min,t_max,(t_max-t_min)/dt)
    all_mass=0
    mass_per_host=np.ones(nhostbins)
    for i in time_bins:
        z=time_to_z(i)

        contrib= 10**Stellar_mass_formed_gal(lMgal,z)
        all_mass=all_mass+contrib
        print 'looking at z=',z,i/1e9,contrib,all_mass
    return all_mass



def Total_mass(zmin):
    mass_density=np.ones(nhostbins)
    mass_gal=np.ones(nhostbins)
    print nhostbins,'nhostbins'
    for j in range(0,nhostbins):
        mass_density[j]=Total_Stellar_mass(lMhostbins[j],zmin)
        mass_gal[j]=Stellar_mass(lMhostbins[j],zmin)
    plt.plot(mass_gal,np.log10(mass_density), label="integral")
    plt.plot(mass_gal,mass_gal, label="y=x")
    plt.grid()
    plt.legend(loc="best")
    plt.xlabel(r"$M_{gal}$")
    plt.ylabel(r"$inetgrated$ $mass$")
    plt.savefig("integral_behroozi_dt8_z"+str(zmin)+".pdf",dpi=400)
    #    plt.show()

    return    


def Total_mass_gal(zmin=.001):
    mass_density=np.ones(ngalbins)
    mass_gal=np.ones(ngalbins)
    for j in range(0,ngalbins):
        mass_density[j]=Total_Stellar_mass_gal(lMgalbins[j],zmin)
    plt.plot(lMgalbins,np.log10(mass_density), label="integral")
    plt.plot(lMgalbins,lMgalbins, label="y=x")
    plt.grid()
    plt.legend(loc="best")
    plt.xlabel(r"$M_{gal}$")
    plt.ylabel(r"$inetgrated$ $mass$")
    plt.savefig("integral_behroozi_dt8_z"+str(zmin)+"_gal.pdf",dpi=400)
    #    plt.show()

    return    



