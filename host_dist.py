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
from cosmological_utils import *


#fname = 'linear_interpolation_SFR_of_zp1_logMh.pkl'
fname = 'logsfr_of_zp1_logMhz0.pkl'
inf = open(fname,'r')
logSFR = pickle.load(inf)

#get press-shechter data from http://hmf.icrar.org
fname_ps='PS.txt'
Mhost_ST,dn_dlogM_ST=np.loadtxt(fname_ps,unpack=True,skiprows=12,usecols=(0,7))
interp_ST=interp1d(Mhost_ST,dn_dlogM_ST)


Mmin=30
Mmax=150
dt=1e8
z_obs_min=0.05



def Period_to_Semimajor(P,M1,M2):
    Mtot = M1 + M2
    return (((phys_const.G*Mtot)/(4.*pi*pi)) * P**2)**(1./3.)


def Semimajor_to_Period(A,M1,M2):
    Mtot = M1 + M2
    return np.sqrt((4*np.pi*np.pi) * A**3 / (phys_const.G*Mtot))


def Stellar_mass(lMhalo,zform):
    #get star formation rate
    #has to take into account metallicity

    logSFR_thishost = logSFR(zform+1,lMhalo)
    
    #    print lMhalo,z,'getting stellar mass'
    #stellar mass as a function of halo mass and redshift
#    lMstar=lMhalo-2 

    total_stellar_mass= 10**logSFR_thishost[0]*dt
    #    print total_stellar_mass,'total stellar mass'
#    print 'a lMhalo=', lMhalo,' forms', logSFR_thishost[0], 'solar masses in stars at z=',zform
    return log10(total_stellar_mass)

def Coalescence_to_Period(t,M1,M2):
    #Initial period for given coalescence time
    Mtot=M1+M2
    mu=(M1*M2)/(M1+M2)
    a= (phys_const.G**3*Mtot**2*mu*t/phys_const.clight**5*256./5.)**(1./4)
    return Semimajor_to_Period(a,M1,M2)


def metallicity(lMgal):
    #oxygen abundance of Sun in Tremonti+2004
    log_O_sun=(8.69-12) 
    #from Tremonti +2004
    log_O_gal=-1.492+1.847*lMgal-0.08026*lMgal**2-12
    metal=log_O_gal-log_O_sun
    #    print metal,10**metal,log_O_gal+12
    return metal_frac=10**metal 

def cdf(x,mu,sigma):
    #cumulative probability of value below x for normal distribution of mean mu 
    # and standard deviation sigma, variance sigma**2
    return .5*(1+erf((x-mu)/sigma/sqrt(2))


def Number_Massive_Stars(lMstar):
    #change to Chabrier IMF
    #change "lower limit" for massive star definition
    #how to count the number of binaries with a massive companion without counting them twice?
    
#Salpeter IMF #ksi(m) dm = ksi0*(m/Msun)**(-2.35)*(dm)/Msun
#ksi(logM)=(m/Msun)**(-1.35)
    lMmin=.1 
    lMmax=2.5
    alpha=-1.35
    #normalisation is set by total stellar mass
    ksi0=lMstar*(-alpha)/(lMmin**(alpha+1)-lMmax**(alpha+1))
    
    #number of massive stars
    lMlim=log10(5.) #look above 5 solar mass stars
    lNstar=ksi0*(lMlim**(alpha)-lMmax**(alpha))/(1-alpha)
    # print 'total number of stars is',10**ksi0,' and number of massive stars is ',10**lNstar
    #    print "total number of stars is {0:e} and number of massive ones is {1:e}".format(10**ksi0,10**lNstar)
    return 10**lNstar

def Number_Binaries(Nstar):
    binary_frac=.5
    Nbinary=Nstar*binary_frac
    return Nbinary


def integrand(M,lMhost,zform):
    M=M*units.Msun
#    print M,'M'    
    # print lMhost,'lMhost'
    # print zform,'zform'
    number_massive=Number_Massive_Stars(Stellar_mass(lMhost,zform))
#    print number_massive,lMhost,zform,'in integrand'
    number_binaries=Number_Binaries(number_massive)
    # fraction of these binaries that can actually merge after with z<.5
    frac=Fraction_Merging(M,zform)
#    print 'out of ',number_binaries,' massive binaries, ', frac*number_binaries,' will actually merge a z<.1'
    
    return number_binaries*frac


def Fraction_Merging(M,zform):
    # this can return <0, check what to do
    #assumes binaris have same mass
    M1=M
    M2=M
    max_coal_time=z_to_Time(z_obs_min) -z_to_Time(zform)
    #    Mcut = 30.  #Msun
    Rstar1 = (M1/YTQuantity(1,'Msun'))**(15./19.) * YTQuantity(1,'Rsun').in_units('AU')     #http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf
    Rstar2 = (M2/YTQuantity(1,'Msun'))**(15./19.) * YTQuantity(1,'Rsun').in_units('AU')     #http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf
#    print "Minimum separation is coming from the constraint that the two stars not touch"
    minR = Rstar1+Rstar2
    minP = Semimajor_to_Period(minR,M1,M2)#.in_units('day')
#    print "Minimum separation comes out to {0}, which corresponds to an orbital period of {1}".format(minR.in_units('AU'),minP.in_units('day'))
    
    #maximum peroid is set by necessity to merge before z=.1
    maxP= Coalescence_to_Period(max_coal_time,M1,M2)
    maxa=Period_to_Semimajor(maxP,M1,M2)
#    print "Maximum separation to have merger is  {0}, which corresponds to an orbital period of {1}".format(maxa.in_units('AU'),maxP.in_units('day'))

    #number of binaries within that period, just integrate Moe+2013 powerlaw
    Pmax_Moe=20*units.day
    Pmin_Moe=2*units.day
    P0=1*units.day
    gam_Moe=-.4
    L_Moe=.22
    kappa_Moe=L_Moe*(gam_Moe-1)/((Pmax_Moe/P0)**(gam_Moe-2)-(Pmin_Moe/P0)**(gam_Moe-2))
    #    non_mormalized=L*(maxP**(gam_Moe-2)-minP**(gam_Moe-2))/(Pmax_Moe**(gam_Moe-2)-Pmin_Moe**(gam_Moe-2)
    #    nmerger=Nstars*((maxP/P0)**(gam_Moe-2)-(minP/P0)**(gam_Moe-2))/((Pmax_Moe/P0)**(gam_Moe-2)-(Pmin_Moe/P0)**(gam_Moe-2))
    
    frac_merger=max(L_Moe/(gam_Moe-1)*((maxP/P0)**(gam_Moe-2)-(minP/P0)**(gam_Moe-2)),0.)
    if isnan(frac_merger):
        frac_merger=0.0
#    print ' about ', 100*frac_merger,' % of the binaries of mass ',M, 'merge after z=.1'
    #non_mormalized=L*(maxP**(gam_Moe-2)-minP**(gam_Moe-2))/(Pmax_Moe**(gam_Moe-2)-Pmin_Moe**(gam_Moe-2)
    return frac_merger




    
def Number_Merging(M1,zform,Nstars):
    #assumes period distribution is indep of mass
    #NUMBER MERGING AT zMIN
    # we have to include mass loss at supernova
    # the maximal period is set by merger within time downto z=1. (let's say Hubble for now)
    # the minimal period is set by no touching of the stars
    # assume period distribution independent of mass, z and host, and metallicity
    # assumes all binaries have P between 2 days and 20. We should nclude all the binaries in order to renormalize properly
    
    M2=M1    
    zmin=.1
    max_coal_time=z_to_Time(zmin) -z_to_Time(zform)
    #    Mcut = 30.  #Msun
    Rstar1 = (M1/YTQuantity(1,'Msun'))**(15./19.) * YTQuantity(1,'Rsun').in_units('AU')     #http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf
    Rstar2 = (M2/YTQuantity(1,'Msun'))**(15./19.) * YTQuantity(1,'Rsun').in_units('AU')     #http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf
    print "Minimum separation is coming from the constraint that the two stars not touch"
    minR = Rstar1+Rstar2
    minP = Semimajor_to_Period(minR,M1,M2)#.in_units('day')
    print "Minimum separation comes out to {0}, which corresponds to an orbital period of {1}".format(minR.in_units('AU'),minP.in_units('day'))

    
    #maximum peroid is set by necessity to merge before z=.1
    maxP= Coalescence_to_Period(max_coal_time,M1,M2)
    maxa=Period_to_Semimajor(maxP,M1,M2)
    print "Maximum separation to have merger is  {0}, which corresponds to an orbital period of {1}".format(maxa.in_units('AU'),maxP.in_units('day'))

    #number of binaries within that period, just integrate Moe+2013 powerlaw
    Pmax_Moe=20*units.day
    Pmin_Moe=2*units.day
    P0=1*units.day
    gam_Moe=-.4
    L_Moe=.22
    kappa_Moe=L_Moe*(gam_Moe-1)/((Pmax_Moe/P0)**(gam_Moe-2)-(Pmin_Moe/P0)**(gam_Moe-2))
    #    non_mormalized=L*(maxP**(gam_Moe-2)-minP**(gam_Moe-2))/(Pmax_Moe**(gam_Moe-2)-Pmin_Moe**(gam_Moe-2)
    #    nmerger=Nstars*((maxP/P0)**(gam_Moe-2)-(minP/P0)**(gam_Moe-2))/((Pmax_Moe/P0)**(gam_Moe-2)-(Pmin_Moe/P0)**(gam_Moe-2))
    nmerger=Nstars *L_Moe/(gam_Moe-1)*((maxP/P0)**(gam_Moe-2)-(minP/P0)**(gam_Moe-2))
#    print 'out of ', Nstars,' stars ',nmerger,' merge after z=.1'
    #non_mormalized=L*(maxP**(gam_Moe-2)-minP**(gam_Moe-2))/(Pmax_Moe**(gam_Moe-2)-Pmin_Moe**(gam_Moe-2)

    
    
    
def All_Mergers(lMhost,zform,Mmin,Mmax):
   # print zform, lMhost, Mmax,Mmin,'test'
    #print Mmin,Mmax,'masses considered for merger'
    merging=quad(integrand,Mmin,Mmax,args=(lMhost,zform))[0]
#    print merging ,'stars between M1=',Mmin,'Msun and M2=',Mmax,' Msun merge if they have formed at zform=',zform ,'in a host of logM*=',lMhost
    return merging

    #return quad(integrand,Mnin,Mmax,args=(zform,LMstar))[0]

lMhost_Min=9
lMhost_Max=15
nhostbins=int((lMhost_Max-lMhost_Min)/.5 +1)
lMhostbins=np.linspace(lMhost_Min,lMhost_Max,num=nhostbins)


def Num_per_Host(zform,Mmin,Mmax):
    num_per_host=[]
    for i in lMhostbins:
#        print i,Mmin,Mmax,zform,'testing'
        num_per_host.append(All_Mergers(i,zform,Mmin,Mmax))
#        print i,  all_mergers(i,zform,Mmin,Mmax),'mergers'
    return num_per_host

def binaries_with_z():
    plt.clf()
    for zform in (1,1.5,2.5,3,3.5,4):
        plt.semilogy(lMhostbins,Num_per_Host(zform,Mmin,Mmax), label="zform={0}".format(zform)) 
    plt.legend(loc="best")
    plt.xlabel(r"$lMhost$")
    plt.ylabel(r"number mergers")
    #pl.savefig("blazar_lum_evolution.pdf",dpi=400)
    plt.show()



#total number of binaries, is the integral over formation times
def Total_mergers(Mmin,Mmax):
    #check for negative contributions
    # we will need to remove the binaries merging at z>.5
    #this is lookback time
    t_min=z_to_Time(10)
    t_max=z_to_Time(0.)
#    print t_min,t_max,'tmin'
    time_bins=np.linspace(t_min,t_max,t_max/dt)
#    print t_min,t_max,time_bins
    all_mergers=0
    lMhost=10
    mergers_per_host=np.ones(nhostbins)
    print nhostbins,'nhostbins'
    for j in range(0,nhostbins):
#        print j,lMhostbins[j],'testing'
        for i in time_bins:
            z=time_to_z(i)
        #get the amount of stars that form during that time and then merge
#            print 'looking at z=',z
            contrib= All_Mergers(lMhostbins[j],z,Mmin,Mmax)
#            print contrib,z,'mergers for z'
            if contrib >0:
                all_mergers=all_mergers+contrib
        print all_mergers,lMhostbins[j],'mergers,host',j                    
        mergers_per_host[j]=all_mergers

    plt.semilogy(lMhostbins,mergers_per_host)#, label="zform={0}".format(zform)) 
    plt.legend(loc="best")
    plt.xlabel(r"$lMhost$")
    plt.ylabel(r"number mergers")
    plt.savefig("mergers_host.pdf",dpi=400)
    plt.show()


    return    


def Total_Mergers_per_Host(Mmin,Mmax,lMhost):
    #check for negative contributions
    # we will need to remove the binaries merging at z>.5
    #this is lookback time
    t_min=z_to_Time(10)
    t_max=z_to_Time(0.)
#    print t_min,t_max,'tmin'
    time_bins=np.linspace(t_min,t_max,t_max/dt)
#    print t_min,t_max,time_bins
    all_mergers=0
    mergers_per_host=np.ones(nhostbins)
    for i in time_bins:
        z=time_to_z(i)
        #get the amount of stars that form during that time and then merge
        print 'looking at z=',z
        contrib= All_Mergers(lMhost,z,Mmin,Mmax)
#            print contrib,z,'mergers for z'
        if contrib >0:
            all_mergers=all_mergers+contrib
    print all_mergers,lMhost,'merger number ,host'
    return all_mergers



def Density_mergers(Mmin,Mmax):
    #number of mergers per Mpc-3, we multiply Number_mergers by halo mass function
    #check for negative contributions
    # we will need to remove the binaries merging at z>.5

    #this is lookback time
#     t_min=z_to_Time(10)
#     t_max=z_to_Time(0.)
# #    print t_min,t_max,'tmin'
#     time_bins=np.linspace(t_min,t_max,t_max/dt)
# #    print t_min,t_max,time_bins
#     all_mergers=0
#     lMhost=10
    merger_density=np.ones(nhostbins)
    print nhostbins,'nhostbins'
    for j in range(0,nhostbins):
#        print j,lMhostbins[j],'testing'

        mergers_per_host=Total_Mergers_per_Host(Mmin,Mmax,lMhostbins[j])
        merger_density[j]=mergers_per_host*interp_ST(10**lMhostbins[j])
        print merger_density[j],lMhostbins[j],'merger density ,host',j                    
    plt.semilogy(lMhostbins,merger_density)#, label="zform={0}".format(zform)) 
    plt.legend(loc="best")
    plt.xlabel(r"$lMhost$")
    plt.ylabel(r"merger density")
    plt.savefig("stars_density.pdf",dpi=400)
    plt.show()


    return    






    

