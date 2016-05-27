#!/bin/python

#import yt.units as units
#import yt.utilities.physical_constants as phys_const
#from yt import YTArray,YTQuantity
from numpy import pi
import numpy as np
#from scipy.integrate import quad
#from scipy.interpolate import interp1d
from math import *
import matplotlib.pyplot as plt
#import pickle
from cosmological_utils import *
from merger_tools import *
from os.path import expanduser
from yt.utilities.cosmology import Cosmology
import matplotlib.ticker as plticker
import matplotlib as mpl
import itertools
dt=1e8
if dt==1e8:
    DT="1e8"
elif dt==1e7:
    DT="1e7"
elif dt==5e8:
    DT="5e8"
elif dt==5e7:
    DT="5e7"
if dt==1e9:
    DT="1e9"


binary_frac=.69

z_min=.1
z_max=7.95
lMgal_Min=7
lMgal_Max=11.25
ldM=.2
ngalbins=int((lMgal_Max-lMgal_Min)/ldM +1)
lMgalbins=np.linspace(lMgal_Min,lMgal_Max,num=ngalbins)

# lMhalo_Min=10
# lMhalo_Max=15
# ldM=.5
# nhalobins=int((lMhalo_Max-lMhalo_Min)/ldM +1)
# lMhalobins=np.linspace(lMhalo_Min,lMhalo_Max,num=nhalobins)


t_min=z_to_Time(z_max)
t_max=z_to_Time(z_min)
time_bins=np.linspace(t_min,t_max,(t_max-t_min)/dt)
ntimebins=len(time_bins)
calib="PP04"
#Zcut=.1  #maximal metallicity we allow: Zcut=Z_gal/Z_sun

#merger rate from Abott+2016:2-53 Gpc-3/yr
 # 2 *10^9 Mpc^(-3)/yr
# dt=1e8 donc 2*10^9/10^8 Mpc(-3)/yr=20-530 /yr


def Total_lowZ_Stellar_mass_gal(lMgal,Zcut):
    #returns total stellar mass formed in galaxy of final mass lMgal
    #should be equal to lMgal but is not because it doesn't have mass loss
    t_min=z_to_Time(z_max)
    t_max=z_to_Time(z_min)
    time_bins=np.linspace(t_min,t_max,(t_max-t_min)/dt)
    all_mass=0
    all_formed=0
    for i in time_bins:
        z=time_to_z(i)
        contrib= 10**Stellar_mass_formed_gal(lMgal,z,dt)
        all_formed=all_formed+contrib
        all_mass=all_mass+contrib*Frac_LowZ_Stars_Avail(log10(all_formed),z,Zcut,calib)
       #        print 'looking at z=',z,i/1e9,contrib,all_mass
    return all_mass




def Frac_LowZ_Stars_Avail(lMgal,zform,Zcut,calib):
    #determines the fraction of lowZ stars for given environemnt
    #will have to include redshift evolution+ improvements Evan
    #account for fact that oxygen is not best indicater of [Fe/H] for certain galaxies
    #need a better idea of scatter

    
    #returns log10(O/H)+12
    #    print lMgal,'lMgal'
    mean_Z=metallicity_z_ma(lMgal,zform,calib)

    #    print mean_Z,10**mean_Z,lMgal,zform
#    mean_Z=metallicity_z_ma(lMgal,zform)
    sigma=.3 #.1 from Tremonti+radial (.2 in Berg+2013) +.2 from Sanders+2012
    prob_low_z=prob(Z_to_OH(Zcut),mean_Z,sigma)
    #    print zform,lMgal,mean_Z,OH_to_Z(mean_Z),Zcut,Z_to_OH(Zcut),prob_low_z
    #    print 'a lMsun=',lMgal,'galaxy has a mean Z of ',mean_Z,', and a P(Z<)',Zcut,' of ',prob_low_z
    return prob_low_z



def Density_mergers(lMgal,Zcut):
    #number density of available low Z stellar mass as a function of galaxy size
    stars_per_gal=Total_lowZ_Stellar_mass_gal(lMgal,Zcut)
    star_density=stars_per_gal*Stellar_Mass_Function(lMgal,ldM,z_min)
    return star_density

    

def make_plot():
    testing_SMF=np.arange(7,11.25,ldM)
    SMF=[]

    #print testing_SMF
    for i in testing_SMF :
        SMF.append(Stellar_Mass_Function(i,ldM,z_min))
        #    print SMF[i],i
        #print SMF    
    plt.semilogy(testing_SMF,SMF)#, label="zform={0}".format(zform)) 
    plt.legend(loc="best")
    plt.xlabel(r"$lMgal$")
    plt.ylabel(r"number density (Mpc$^{-3}$")
    plt.savefig("Galaxy_SMF.pdf",dpi=400)

def get_distrib(lMgal,Zcut):
    #provides merging binaries as a function of redhsift of formation
    z_bins=np.zeros(len(time_bins))
    all_mergers=np.zeros(len(time_bins))
    all_mass=0
    lMhalo=Gal_to_Halo_mass(lMgal,z_min)
    dlMh=1e-3
    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])

#        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        lowZ_stars_in_halo=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):
            #get corresponding galaxy mass through abundance matching 
            lMsubgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt*dlMh
            #multiply by fraction in mass bin
            stellar_mass_formed=stellar_mass_formed*mass_frac[m]
            #get amount of lowZ stuff
            lowZ_stars_in_halo[m]=stellar_mass_formed*Frac_LowZ_Stars_Avail(lMsubgal,z_bins[i],Zcut,calib)

        #now add up contributions of subhalos
        all_lowZ_stars=np.sum(lowZ_stars_in_halo)
        time_delay=t_max-time_bins[i]
        if time_delay != 0: 
            if Zcut==.3:
                contrib_merger=all_lowZ_stars*interp_merg_03(log10(time_delay))
            elif Zcut==.1:
                contrib_merger=all_lowZ_stars*interp_merg_01(log10(time_delay))
            elif Zcut==.01:
                contrib_merger=all_lowZ_stars*interp_merg_001(log10(time_delay))
            else:
                print "we do not model this metallicity"
        all_mergers[i]=contrib_merger
    return all_mergers,z_bins


def get_distrib_noEPS(lMgal,Zcut):
    #provides merging binaries as a function of redhsift of formation
    z_bins=np.zeros(len(time_bins))
    all_mergers=np.zeros(len(time_bins))
    all_mass=0
    lMhalo=Gal_to_Halo_mass(lMgal,z_min)
    dlMh=1e-3
    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])
        stellar_mass_formed= SFR_lMhalo(lMhalo,z_bins[i])*dt*dlMh
        all_lowZ_stars=stellar_mass_formed*Frac_LowZ_Stars_Avail(lMgal,z_bins[i],Zcut,calib)
        time_delay=t_max-time_bins[i]
        if time_delay != 0: 
            if Zcut==.3:
                contrib_merger=all_lowZ_stars*interp_merg_03(log10(time_delay))
            elif Zcut==.1:
                contrib_merger=all_lowZ_stars*interp_merg_01(log10(time_delay))
            elif Zcut==.01:
                contrib_merger=all_lowZ_stars*interp_merg_001(log10(time_delay))
            else:
                print "we do not model this metallicity"
        all_mergers[i]=contrib_merger*binary_frac
    return all_mergers,z_bins


def get_distrib_full(lMgal,Zcut):
    #provides merging binaries as a function of redhsift of formation for zmerg=1;.0.5;.1
    mass_z1=1e-10
    mass_z05=1e-10
    z_bins=np.zeros(len(time_bins))
    all_mergers=np.zeros(len(time_bins))
    all_mergers_z1=np.zeros(len(time_bins))
    all_mergers_z05=np.zeros(len(time_bins))
    all_stars=np.zeros(len(time_bins))
    all_mass=0
    lMhalo=Gal_to_Halo_mass(lMgal,z_min)
    dlMh=1e-3
    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])
        
#        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        local_stellar_mass_formed=np.zeros(len(mass_ratios))
        lowZ_stars_in_halo=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):
            #get corresponding galaxy mass through abundance matching 
            lMsubgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt*dlMh
            #multiply by fraction in mass bin
            stellar_mass_formed=stellar_mass_formed*mass_frac[m]
            #get amount of lowZ stuff
            lowZ_stars_in_halo[m]=stellar_mass_formed*Frac_LowZ_Stars_Avail(lMsubgal,z_bins[i],Zcut,calib)
            if np.isnan(mass_frac[m]):
                mass_frac[m]=0
                local_stellar_mass_formed[m]=0
        all_stars[i]=max(np.sum(stellar_mass_formed),1e-10)
        #now add up contributions of subhalos
        all_lowZ_stars=np.sum(lowZ_stars_in_halo)

#        if i>1:
        if (z_bins[i]>1) :
                #store z=1 case
            time_delay=z_to_Time(1.) -time_bins[i]
            if time_delay>10**6.7:
                #remove case not modeled by Drew
                if Zcut==.3:
                    contrib_merger_z1=all_lowZ_stars*interp_merg_03(log10(time_delay))*binary_frac
                elif Zcut==.1:
                    contrib_merger_z1=all_lowZ_stars*interp_merg_01(log10(time_delay))*binary_frac
                elif Zcut==.01:
                    contrib_merger_z1=all_lowZ_stars*interp_merg_001(log10(time_delay))*binary_frac
                else:
                    print "we do not model this metallicity"
            else:         
                fname_merger_t='bhbh_mergers_dNdt_SGK_m_150_ecc_FIXED.dat'
                tdelay,merger_frac_001,merger_frac_01,merger_frac_03=np.loadtxt(fname_merger_t,unpack=True,skiprows=1)
                if Zcut==.3:
                    contrib_merger_z1=all_lowZ_stars*merger_frac_03[0]*binary_frac
                elif Zcut==.1:
                    contrib_merger_z1=all_lowZ_stars*merger_frac_01[0]*binary_frac
                elif Zcut==.01:
                    contrib_merger_z1=all_lowZ_stars*merger_frac_001[0]*binary_frac
                else:
                    print "we do not model this metallicity"
                
            all_mergers_z1[i]=contrib_merger_z1
        mass_z1=all_stars[i]


        if (z_bins[i]>.5) :
                #store z=.5 case
            time_delay=z_to_Time(.5) -time_bins[i]
            if time_delay>10**6.7:
                #remove case not modeled by Drew
                if Zcut==.3:
                    contrib_merger_z05=all_lowZ_stars*interp_merg_03(log10(time_delay))*binary_frac
                elif Zcut==.1:
                    contrib_merger_z05=all_lowZ_stars*interp_merg_01(log10(time_delay))*binary_frac
                elif Zcut==.01:
                    contrib_merger_z05=all_lowZ_stars*interp_merg_001(log10(time_delay))*binary_frac
                else:
                    print "we do not model this metallicity"
            else:         
                fname_merger_t='bhbh_mergers_dNdt_SGK_m_150_ecc_FIXED.dat'
                tdelay,merger_frac_001,merger_frac_01,merger_frac_03=np.loadtxt(fname_merger_t,unpack=True,skiprows=1)
                if Zcut==.3:
                    contrib_merger_z05=all_lowZ_stars*merger_frac_03[0]*binary_frac
                elif Zcut==.1:
                    contrib_merger_z05=all_lowZ_stars*merger_frac_01[0]*binary_frac
                elif Zcut==.01:
                    contrib_merger_z05=all_lowZ_stars*merger_frac_001[0]*binary_frac
                else:
                    print "we do not model this metallicity"
            all_mergers_z05[i]=contrib_merger_z05
        mass_z05=all_stars[i]




        time_delay=t_max-time_bins[i]
        if time_delay != 0: 
            if Zcut==.3:
                contrib_merger=all_lowZ_stars*interp_merg_03(log10(time_delay))*binary_frac
            elif Zcut==.1:
                contrib_merger=all_lowZ_stars*interp_merg_01(log10(time_delay))*binary_frac
            elif Zcut==.01:
                contrib_merger=all_lowZ_stars*interp_merg_001(log10(time_delay))*binary_frac
            else:
                print "we do not model this metallicity"
        all_mergers[i]=contrib_merger
    return all_mergers,z_bins,all_mergers_z1,all_mergers_z05,np.log10(mass_z1),np.log10(mass_z05)





def get_distrib_with_cut(lMgal,Zcut):
    #provides merging binaries as a function of redhsift of formation
    z_bins=np.zeros(len(time_bins))
    all_mergers=np.zeros(len(time_bins))
    all_mass=0
    lMhalo=Gal_to_Halo_mass(lMgal,z_min)

    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])
        dlMh=1e-3
        #loop over all subhalo masses
        lsubhM_min=-4.95
        lsubhM_max=-.05
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        lowZ_stars_in_halo=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):

            #get corresponding galaxy mass through abundance matching 
            lMsubgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt*dlMh
            #multiply by fraction in mass bin
            stellar_mass_formed=stellar_mass_formed*mass_frac[m]            
            #get amount of lowZ stuff
            lowZ_stars_in_halo[m]=stellar_mass_formed*Frac_LowZ_Stars_Avail(lMsubgal,z_bins[i],Zcut,calib)
            if lowZ_stars_in_halo[m]<5e3:
                lowZ_stars_in_halo[m]=0.


        #now add up contributions of subhalos
        all_lowZ_stars=np.sum(lowZ_stars_in_halo)
        time_delay=t_max-time_bins[i]
        if time_delay != 0: 
            if Zcut==.3:
                contrib_merger=all_lowZ_stars*interp_merg_03(log10(time_delay))
            elif Zcut==.1:
                contrib_merger=all_lowZ_stars*interp_merg_01(log10(time_delay))
            elif Zcut==.01:
                contrib_merger=all_lowZ_stars*interp_merg_001(log10(time_delay))
            else:
                print "we do not model this metallicity"
        all_mergers[i]=contrib_merger
    return all_mergers,z_bins


def get_distrib_lowZ_stars(lMgal,Zcut):
    #provides low Z stellar mass as a function of redhsift of formation
    z_bins=np.zeros(len(time_bins))
    lowZ_stars=np.zeros(len(time_bins))
    lMhalo=Gal_to_Halo_mass(lMgal,z_min)
    dlMh=1e-3
    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])
#        print z_bins[i],time_bins[i]/1e9,'t in distrib'
        #loop over all subhalo masses
#        print z_bins[i],time_bins[i]/1e9,dt/1e9,'z,t,dt'
        lsubhM_min=-4.95
        lsubhM_max=-.05
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        lowZ_stars_in_halo=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):
            #get corresponding galaxy mass through abundance matching 
            lMgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt*dlMh
            #multiply by fraction in mass bin
            stellar_mass_formed=stellar_mass_formed*mass_frac[m]
            #get amount of lowZ stuff
            lowZ_stars_in_halo[m]=stellar_mass_formed*Frac_LowZ_Stars_Avail(lMgal,z_bins[i],Zcut,calib)
        lowZ_stars[i]=np.sum(lowZ_stars_in_halo)
        
#        print lowZ_stars
    return lowZ_stars,z_bins


def get_distrib_stars(lMgal):
    z_bins=np.zeros(len(time_bins))
    all_stars=np.zeros(len(time_bins))
    all_mass=0

    lMhalo=Gal_to_Halo_mass(lMgal,z_min)
    print z_min,lMhalo
    dlMh=1e-3
    for i in range (0,len(time_bins)):
        
        z_bins[i]=time_to_z(time_bins[i])
        print z_bins[i]
#        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)\
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i],dlMh,z_min)
        local_stellar_mass_formed=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):

            #get corresponding galaxy mass through abundance matching 
            #lMgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
#            print SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i]),'SFR'

            local_stellar_mass_formed[m]= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt*dlMh
            if np.isnan(local_stellar_mass_formed[m]):
                local_stellar_mass_formed[m]=0.
            # if np.isnan(mass_frac[m]):
            #     print 'pb'
            #multiply by fraction in mass bin
#            print local_stellar_mass_formed[m],mass_frac[m]
#            print local_stellar_mass_formed[m],mass_frac[m]
            local_stellar_mass_formed[m]=local_stellar_mass_formed[m] *mass_frac[m]
 #           print local_stellar_mass_formed[m],mass_frac[m]
#            print np.log10(np.sum(stellar_mass_formed)),'stellar mass formed at z=',z_bins[i]
            if np.isnan(mass_frac[m]):
                mass_frac[m]=0
                local_stellar_mass_formed[m]=0
        all_stars[i]=max(np.sum(local_stellar_mass_formed),1e-10)
#        print np.log10(all_stars[i]),'stellar mass formed at z=',z_bins[i],local_stellar_mass_formed
    print np.log10(np.sum(all_stars)),'all mass formed'    
    return #all_stars,z_bins



def get_all_mergers_allz(Zcut):

    #make data for EPS
    size=ngalbins*ntimebins
    all_mgal=np.zeros(size)
    all_mgal_z1=np.zeros(size)
    all_mgal_z05=np.zeros(size)
    all_zform=np.zeros(size)
    all_rates=np.zeros(size)
    all_rates_z1=np.zeros(size)
    all_rates_z05=np.zeros(size)
    lmass_z1=np.zeros(ngalbins)
    lmass_z05=np.zeros(ngalbins)
    all_mergers=np.empty((ngalbins,ntimebins))
    all_mergers_z1=np.empty((ngalbins,ntimebins))
    all_mergers_z05=np.empty((ngalbins,ntimebins))
    for m in range(0,ngalbins):
        mergers,z_bins,mergers_z1,mergers_z05,lmass_z1[m],lmass_z05[m]=get_distrib_full(lMgalbins[m],Zcut)
        
        all_mergers[m,:]=mergers*Stellar_Mass_Function(lMgalbins[m],ldM,z_min)
        all_mergers_z1[m,:]=mergers_z1*Stellar_Mass_Function(lmass_z1[m],ldM,1)
        all_mergers_z05[m,:]=mergers_z05*Stellar_Mass_Function(lmass_z05[m],ldM,0.5)
    print lmass_z05,lmass_z1,mergers_z05,mergers_z1
    for m in range(0,ngalbins):
        for t in range(0,ntimebins):
            r=m*ntimebins+t
            all_mgal[r]=lMgalbins[m]
            all_mgal_z1[r]=lmass_z1[m]
            all_mgal_z05[r]=lmass_z05[m]

            all_zform[r]=z_bins[t]
            all_rates[r]=all_mergers[m,t]
            all_rates_z1[r]=all_mergers_z1[m,t]
            all_rates_z05[r]=all_mergers_z05[m,t]
    np.save("lMgal_zform_z01"+str(Zcut)+str(calib)+"_dt"+DT,(all_mgal,all_zform,all_rates))
    np.save("lMgal_zform_z05"+str(Zcut)+str(calib)+"_dt"+DT,(all_mgal_z05,all_zform,all_rates_z05))
    np.save("lMgal_zform_z1"+str(Zcut)+str(calib)+"_dt"+DT,(all_mgal_z1,all_zform,all_rates_z1))
    return 


def get_all_mergers(Zcut):

    #make data for Shea
    size=ngalbins*ntimebins
#    print size,'size'
    all_mgal=np.zeros(size)
    all_zform=np.zeros(size)
    all_rates=np.zeros(size)
    all_mergers=np.empty((ngalbins,ntimebins))
#    print lMgalbins,'galbins'
#    print len(all_mgal)
#    all_stars=np.zeros(ngalbins,ntimebins)
    for m in range(0,ngalbins):
        # for t in ntimebins:
#        z[m]=time_to_z(time_bins[t])
#        lMhalo=Gal_to_Halo_mass(lMgalbins[m])
        mergers,z_bins=get_distrib(lMgalbins[m],Zcut)
        all_mergers[m,:]=mergers*Stellar_Mass_Function(lMgalbins[m],ldM,z_min)#/dt*1e9 # convert to Gpc^-3
#        print Stellar_Mass_Function(lMgalbins[m],ldM)
#        print mergers

    for m in range(0,ngalbins):
        for t in range(0,ntimebins):
            r=m*ntimebins+t
            all_mgal[r]=lMgalbins[m]
            all_zform[r]=z_bins[t]
            all_rates[r]=all_mergers[m,t]
#            np.save("lMgal_zform_rates_EPSz_1_Z"+str(Zcut)+str(calib)+"_dt"+DT,(all_mgal,all_zform,all_rates))
            np.save("lMgal_zform_rates_EPS_test_Z"+str(Zcut)+str(calib)+"_dt"+DT,(all_mgal,all_zform,all_rates))
#    print np.sum(all_mergers)*ldM*dt,np.sum(all_mergers,axis=0)*M,np.sum(all_mergers,axis=1)*dt,'summing'
#    print np.sum(all_mergers)*ldM*dt,'summing',np.mean(all_mergers)/ldM*dt            

    return #all_mergers,z_bins



def get_all_mergers_noEPS(Zcut):


    size=ngalbins*ntimebins
    all_mgal=np.zeros(size)
    all_zform=np.zeros(size)
    all_rates=np.zeros(size)
    all_mergers=np.empty((ngalbins,ntimebins))
    for m in range(0,ngalbins):
        mergers,z_bins=get_distrib_noEPS(lMgalbins[m],Zcut)
        all_mergers[m,:]=mergers*Stellar_Mass_Function(lMgalbins[m],ldM,z_min)#/dt*1e9 # convert to Gpc^-3

    for m in range(0,ngalbins):
        for t in range(0,ntimebins):
            r=m*ntimebins+t
            all_mgal[r]=lMgalbins[m]
            all_zform[r]=z_bins[t]
            all_rates[r]=all_mergers[m,t]
            np.save("lMgal_zform_noEPS_Z"+str(Zcut)+str(calib)+"_dt"+DT,(all_mgal,all_zform,all_rates))

    return #all_mergers,z_bins


def check_mass():
#    size=ngalbins*ntimebins
#    print size,'size'
#    all_mgal=np.zeros(size)
#    all_zform=np.zeros(size)
#    all_rates=np.zeros(size)
#    all_mgal=np.empty(ngalbins,ntimebins)
    print lMgalbins,'galbins'
#    print len(all_mgal)
#    all_stars=np.zeros(ngalbins,ntimebins)
    for m in range(0,ngalbins):
        # for t in ntimebins:
#        z[m]=time_to_z(time_bins[t])
#        lMhalo=Gal_to_Halo_mass(lMgalbins[m])
        stars,z_bins=get_distrib_stars(lMgalbins[m])
#        all_gal[m,:]=mergers*Stellar_Mass_Function(lMgalbins[m],ldM)
#        print Stellar_Mass_Function(lMgalbins[m],ldM)
#        print mergers
    # for m in range(0,ngalbins):
    #     for t in range(0,ntimebins):
    #         r=m*ntimebins+t
    #         all_mgal[r]=lMgalbins[m]
    #         all_zform[r]=z_bins[t]
    #         all_rates[r]=all_mergers[m,t]
    # np.save("lMgal_zform_rates_Z_dt1e7_dlm02"+str(Zcut)+str(calib),(all_mgal,all_zform,all_rates))
        print "a ", lMgalbins[m],"galaxy has" ,np.log10(np.sum(stars)),"formed"
#    return all_gal,z_bins
    



def make_merger_formation_plot():
#    stellar_dens=np.arange(7,11.25,1)
    #merger_dens_01=[]
    #merger_dens_05=[]
    #merger_dens_001=[]


    merg_7_03,z=get_distrib_with_cut(7.,.3) 
    merg_8_03,z=get_distrib_with_cut(8.,.3)
    merg_9_03,z=get_distrib_with_cut(9.,.3)
    merg_10_03,z=get_distrib_with_cut(10.,.3) 
    merg_11_03,z=get_distrib_with_cut(11.,.3)
    print merg_7_03,merg_8_03,merg_9_03,merg_10_03,merg_10_03

    merg_7_01,z=get_distrib_with_cut(7.,.1)
    merg_8_01,z=get_distrib_with_cut(8.,.1)
    merg_9_01,z=get_distrib_with_cut(9.,.1)
    merg_10_01,z=get_distrib_with_cut(10.,.1)
    merg_11_01,z=get_distrib_with_cut(11.,.1)
    print merg_7_01,merg_8_01,merg_9_01,merg_10_01,merg_10_01

    merg_7_001,z=get_distrib_with_cut(7.,.01)
    merg_8_001,z=get_distrib_with_cut(8.,.01)
    merg_9_001,z=get_distrib_with_cut(9.,.01)
    merg_10_001,z=get_distrib_with_cut(10.,.01)
    merg_11_001,z=get_distrib_with_cut(11.,.01)
#    print merg_7_001,merg_8_001,merg_9_001,merg_10_001,merg_10_001
#     merg_7_03,z=get_distrib(7.,.3) 
# #   print np.sum(merg_7_03),'all mergers'
# #   print merg_7_03,'per time'
#     merg_8_03,z=get_distrib(8.,.3)
#     merg_9_03,z=get_distrib(9.,.3)
#     merg_10_03,z=get_distrib(10.,.3)
#     merg_11_03,z=get_distrib(11.,.3)
    

#     merg_7_01,z=get_distrib(7.,.1)
#     merg_8_01,z=get_distrib(8.,.1)
#     merg_9_01,z=get_distrib(9.,.1)
#     merg_10_01,z=get_distrib(10.,.1)
#     merg_11_01,z=get_distrib(11.,.1)


#     merg_7_001,z=get_distrib(7.,.01)
#     merg_8_001,z=get_distrib(8.,.01)
#     merg_9_001,z=get_distrib(9.,.01)
#     merg_10_001,z=get_distrib(10.,.01)
#     merg_11_001,z=get_distrib(11.,.01)
    

    plt.clf()
    #divide by ten to get rates in events/yr/Gpc

    fig=plt.figure(figsize=(8.5,7.5))
    ax1 = fig.add_subplot(111)
    lookback_time=[(z_to_Time(0)-z_to_Time(i))/1e9 for i in z ]

    ax1.semilogy(lookback_time,(merg_7_01*Stellar_Mass_Function(7,1)*1e9,z_min) ,label='10$^{7}M_{\odot}$',color='k',linewidth=2)
    ax1.semilogy(lookback_time,(merg_7_001*Stellar_Mass_Function(7,1,z_min)*1e9),':',color='k',linewidth=2)
    ax1.semilogy(lookback_time,(merg_7_03*Stellar_Mass_Function(7,1,z_min)*1e9),'--',color='k',linewidth=2)

    ax1.semilogy(lookback_time,(merg_8_01*Stellar_Mass_Function(8,1,z_min)*1e9),label='10$^{8}M_{\odot}$',color='b',linewidth=2)
    ax1.semilogy(lookback_time,(merg_8_001*Stellar_Mass_Function(8,1,z_min)*1e9),':',color='b',linewidth=2)
    ax1.semilogy(lookback_time,(merg_8_03*Stellar_Mass_Function(8,1,z_min)*1e9),'--',color='b',linewidth=2)

    ax1.semilogy(lookback_time,(merg_9_01*Stellar_Mass_Function(9,1,z_min)*1e9),label='10$^{9}M_{\odot}$',color='c',linewidth=2)
    ax1.semilogy(lookback_time,(merg_9_001*Stellar_Mass_Function(9,1,z_min)*1e9),':',color='c',linewidth=2)
    ax1.semilogy(lookback_time,(merg_9_03*Stellar_Mass_Function(9,1,z_min)*1e9),'--',color='c',linewidth=2)
    
    ax1.semilogy(lookback_time,(merg_10_01*Stellar_Mass_Function(10,1,z_min)*1e9),label='10$^{10}M_{\odot}$',color='g',linewidth=2)
    ax1.semilogy(lookback_time,(merg_10_001*Stellar_Mass_Function(10,1,z_min)*1e9),':',color='g',linewidth=2)
    ax1.semilogy(lookback_time,(merg_10_03*Stellar_Mass_Function(10,1,z_min)*1e9),'--',color='g',linewidth=2)
    

    ax1.semilogy(lookback_time,(merg_11_01*Stellar_Mass_Function(11,1,z_min)*1e9),label='10$^{11}M_{\odot}$',color='r',linewidth=2)
    ax1.semilogy(lookback_time,(merg_11_001*Stellar_Mass_Function(11,1,z_min)*1e9),':',color='r',linewidth=2)
    ax1.semilogy(lookback_time,(merg_11_03*Stellar_Mass_Function(11,1,z_min)*1e9),'--',color='r',linewidth=2)


    ax2 = ax1.twiny()
    zloc = [.1,.2,.5,1,2,3,4,5,7]

    tloc = [(z_to_Time(0)-z_to_Time(zl))/1e9 for zl in zloc]

    zstr = ['$'+str(i)+'$' for i in zloc]


#    print lookback_time,z
#    print tloc, zloc

#    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xlim(tloc[0],tloc[-1])
    ax2.set_xticks(tloc)
    ax2.set_xticklabels(zloc)
    fontsize=18
    ax2.set_xlabel(r'$z_{form}$',fontsize=fontsize)

    ax1.tick_params(labelsize=18)
    ax1.set_xlim(1.3,13)
#    plt.xlim(0,7)
    ax1.legend(bbox_to_anchor=(10, 1e8),bbox_transform=plt.gcf().transFigure)
#    ax1.legend(loc="best")
    ax1.legend(loc=4)

    ax1.set_ylim(1e-2,1e2)
    ax1.set_xlabel(r"$t_{lookback}^{form}$ $(Gyr)$",fontsize=18)
    ax1.set_ylabel(r"$BH$ $merger$ $rate$ $(yr^{-1} Gpc^{-3})$",fontsize=18)
#    plt.savefig("test.pdf")
    plt.savefig("merger_formation_density_withsubhalos_withcut5e3"+calib+"_dt"+DT+".pdf",dpi=400)
#    plt.savefig("merger_formation_density_withsubhalos_"+calib+"_dt"+DT+".pdf",dpi=400)
    plt.clf()



def make_lowZstar_formation_plot():
    #convert into Mpc-3 Gyr-1    
    merg_7_03,z=get_distrib_lowZ_stars(7.,.3)
    merg_8_03,z=get_distrib_lowZ_stars(8.,.3)
    merg_9_03,z=get_distrib_lowZ_stars(9.,.3)
    merg_10_03,z=get_distrib_lowZ_stars(10.,.3)
    merg_11_03,z=get_distrib_lowZ_stars(11.,.3)

    merg_7_01,z=get_distrib_lowZ_stars(7.,.1)
    merg_8_01,z=get_distrib_lowZ_stars(8.,.1)
    merg_9_01,z=get_distrib_lowZ_stars(9.,.1)
    merg_10_01,z=get_distrib_lowZ_stars(10.,.1)
    merg_11_01,z=get_distrib_lowZ_stars(11.,.1)

    merg_7_001,z=get_distrib_lowZ_stars(7.,.01)
    merg_8_001,z=get_distrib_lowZ_stars(8.,.01)
    merg_9_001,z=get_distrib_lowZ_stars(9.,.01)
    merg_10_001,z=get_distrib_lowZ_stars(10.,.01)
    merg_11_001,z=get_distrib_lowZ_stars(11.,.01)

    lookback_time=[(z_to_Time(0)-z_to_Time(i))/1e9 for i in z ]

    #store data so that Shea can play with it
    d7_03=merg_7_03*Stellar_Mass_Function(7,1,z_min)*1e9/dt
    d8_03=merg_8_03*Stellar_Mass_Function(8,1,z_min)*1e9/dt
    d9_03=merg_9_03*Stellar_Mass_Function(9,1,z_min)*1e9/dt
    d10_03=merg_10_03*Stellar_Mass_Function(10,1,z_min)*1e9/dt
    d11_03=merg_11_03*Stellar_Mass_Function(11,1,z_min)*1e9/dt

    d7_01=merg_7_01*Stellar_Mass_Function(7,1,z_min)*1e9/dt
    d8_01=merg_8_01*Stellar_Mass_Function(8,1,z_min)*1e9/dt
    d9_01=merg_9_01*Stellar_Mass_Function(9,1,z_min)*1e9/dt
    d10_01=merg_10_01*Stellar_Mass_Function(10,1,z_min)*1e9/dt
    d11_01=merg_11_01*Stellar_Mass_Function(11,1,z_min)*1e9/dt

    d7_001=merg_7_001*Stellar_Mass_Function(7,1,z_min)*1e9/dt
    d8_001=merg_8_001*Stellar_Mass_Function(8,1,z_min)*1e9/dt
    d9_001=merg_9_001*Stellar_Mass_Function(9,1,z_min)*1e9/dt
    d10_001=merg_10_001*Stellar_Mass_Function(10,1,z_min)*1e9/dt
    d11_001=merg_11_001*Stellar_Mass_Function(11,1,z_min)*1e9/dt


    np.save("lowZstars_zform_tform",(z,lookback_time,d7_03,d8_03,d9_03,d10_03,d11_03,d7_01,d8_01,d9_01,d10_01,d11_01,d7_001,d8_001,d9_001,d10_001,d11_001))

    colors = ['k','b', 'c', 'g', 'r']
    cc=itertools.cycle(colors)

    print "Making plot..."
    plot_lines=[]
    plt.clf()
    fig=plt.figure(figsize=(9.5,8))
    ax1 = fig.add_subplot(111)

    masslines = []



    c=next(cc)
    l, =ax1.semilogy(lookback_time,merg_7_01*Stellar_Mass_Function(7,1,z_min),label='10$^{7} M_{\odot}$',color=c,linewidth=2)
    l2, =ax1.semilogy(lookback_time,merg_7_001*Stellar_Mass_Function(7,1,z_min),':',color=c,linewidth=2)
    l3, =ax1.semilogy(lookback_time,merg_7_03*Stellar_Mass_Function(7,1,z_min),'--',color=c,linewidth=2)

    masslines.append(l)
    metallines = [l,l2,l3]

    c=next(cc)
    ln, =ax1.semilogy(lookback_time,merg_8_01*Stellar_Mass_Function(8,1,z_min),label='10$^{8}M_{\odot}$',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_8_001*Stellar_Mass_Function(8,1,z_min),':',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_8_03*Stellar_Mass_Function(8,1,z_min),'--',color=c,linewidth=2)
    
    masslines.append(ln)

    c=next(cc)
    ln, =ax1.semilogy(lookback_time,merg_9_01*Stellar_Mass_Function(9,1,z_min),label='10$^{9}M_{\odot}$',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_9_001*Stellar_Mass_Function(9,1,z_min),':',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_9_03*Stellar_Mass_Function(9,1,z_min),'--',color=c,linewidth=2)
    
    masslines.append(ln)

    c=next(cc)
    ln, =ax1.semilogy(lookback_time,merg_10_01*Stellar_Mass_Function(10,1,z_min),label='10$^{10}M_{\odot}$',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_10_001*Stellar_Mass_Function(10,1,z_min),':',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_10_03*Stellar_Mass_Function(10,1,z_min),'--',color=c,linewidth=2)
    
    masslines.append(ln)

    c=next(cc)
    ln, =ax1.semilogy(lookback_time,merg_11_01*Stellar_Mass_Function(11,1,z_min),label='10$^{11}M_{\odot}$',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_11_001*Stellar_Mass_Function(11,1,z_min),':',color=c,linewidth=2)
    ax1.semilogy(lookback_time,merg_11_03*Stellar_Mass_Function(11,1,z_min),'--',color=c,linewidth=2)
    
    masslines.append(ln)



    # ax1.semilogy(lookback_time,(merg_7_01*Stellar_Mass_Function(7,1,z_min)),label='10$^{7} M_{\odot}$',color='k',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_7_001*Stellar_Mass_Function(7,1,z_min)),':',color='k',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_7_03*Stellar_Mass_Function(7,1,z_min)),'--',color='k',linewidth=2)

    # ax1.semilogy(lookback_time,(merg_8_01*Stellar_Mass_Function(8,1,z_min)),label='10$^{8}M_{\odot}$',color='b',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_8_001*Stellar_Mass_Function(8,1,z_min)),':',color='b',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_8_03*Stellar_Mass_Function(8,1,z_min)),'--',color='b',linewidth=2)

    # ax1.semilogy(lookback_time,(merg_9_01*Stellar_Mass_Function(9,1,z_min)),label='10$^{9}M_{\odot}$',color='c',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_9_001*Stellar_Mass_Function(9,1,z_min)),':',color='c',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_9_03*Stellar_Mass_Function(9,1,z_min)),'--',color='c',linewidth=2)

    # ax1.semilogy(lookback_time,(merg_10_01*Stellar_Mass_Function(10,1,z_min)),label='10$^{10}M_{\odot}$',color='g',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_10_001*Stellar_Mass_Function(10,1,z_min)),':',color='g',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_10_03*Stellar_Mass_Function(10,1,z_min)),'--',color='g',linewidth=2)

    # ax1.semilogy(lookback_time,(merg_11_01*Stellar_Mass_Function(11,1,z_min)),label='10$^{11}M_{\odot}$',color='r',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_11_001*Stellar_Mass_Function(11,1,z_min)),':',color='r',linewidth=2)
    # ax1.semilogy(lookback_time,(merg_11_03*Stellar_Mass_Function(11,1,z_min)),'--',color='r',linewidth=2)

    print "Making plot look nice!"

    legend1=plt.legend(masslines,[l.get_label() for l in masslines],loc=4)
    #plt.legend([l[0] for l in plot_lines],lookback_time,loc=4)
#    legend2 = plt.legend([metallines[1],metallines[0],metallines[2]],[r'$Z=0.01Z_\odot$',r'$Z=0.1Z_\odot$',r'$Z=0.3Z_\odot$'],loc=4)
    legend2 = plt.legend([metallines[1],metallines[0],metallines[2]],[r'$0.01Z_\odot$',r'$0.1Z_\odot$',r'$0.3Z_\odot$'],loc=8 )

    ax1.add_artist(legend1)
    ax2 = ax1.twiny()
    zloc = [.1,.2,.5,1,2,3,4,7]
    ztloc = [(z_to_Time(0)-z_to_Time(zl))/1e9 for zl in zloc]
    zstr = ['$'+str(i)+'$' for i in zloc]

    tloc = [2,4,6,8,10,12]
    tstr = ['$'+str(i)+'$' for i in tloc]

#    axwidth=3
#    axlength=12    
#    yloc = [1,10,1e2,1e3,1e4,1e5,1e6]
    #    yiloc = [0,1,2,3,4,5,6]
#    ystr = ['$'+str(i)+' $' for i in yloc]
#    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,6))
#     ax1.set_xticks(tloc)
#     ax1.set_xticklabels(tstr)
# #    ax1.set_yticklabels(ystr)
#     ax1.tick_params('x', length=14, width=2, which='major')
#     ax1.tick_params('x', length=8, width=1, which='minor')
#     ax1.tick_params('y', length=14, width=2, which='major')
#     ax1.tick_params('y', length=8, width=1, which='minor')
#     ax1.tick_params(labelsize=24)
    ax1.set_ylim(1e1,5e6)
    ax1.set_xlim(1.3,13)
#   
    
    #ax2.set_xlim(ax1.get_xlim())
    ax2.set_xlim(ztloc[0],ztloc[-1])
    ax2.set_xticks(ztloc)
#    ax2.set_xticklabels(zstr)
    ax2.set_xticklabels(zloc)

    #    ax2.tick_params('x', length=8, width=1, which='minor')
    
    #ax2.tick_params('both', length=14, width=2, which='major')
    #ax2.tick_params(labelsize=24)
        
    #plt.xlim(0,7)
#    ax1.legend(bbox_to_anchor=(10, 1e8),bbox_transform=plt.gcf().transFigure)
#    ax1.legend(loc="best")
#    ax1.legend(loc=4)#,fontsize=24)
    ax2.set_xlabel(r'$z_{\mathrm{form}}$')#,fontsize=28)
    ax1.set_xlabel(r"$t_{\mathrm{form}}$ $\mathrm{(Gyr)}$")#,fontsize=22)
    ax1.set_ylabel(r"$\frac{d\log M_* }{d\log M_{gal} dt_{form}}(M_{\odot}$  $\mathrm{Mpc^{-3}} Gyr^{-1})$")#,fontsize=22)
    plt.savefig("lowZ_star_formation_evol_density_lookback"+calib+"_dt"+DT+".pdf",dpi=400,bbox_inches='tight')

    plt.clf()




def make_star_formation_plot():

    merg_7_01,z=get_distrib_stars(7.,.1)
    merg_8_01,z=get_distrib_stars(8.,.1)
    merg_9_01,z=get_distrib_stars(9.,.1)
    merg_10_01,z=get_distrib_stars(10.,.1)
    merg_11_01,z=get_distrib_stars(11.,.1)


    plt.clf()
    fig=plt.figure(figsize=(8.5,7.5))
    ax1 = fig.add_subplot(111)


    lookback_time=[(z_to_Time(0)-z_to_Time(i))/1e9 for i in z ]

    ax1.semilogy(lookback_time,(merg_7_01*Stellar_Mass_Function(7,1)),label='10$^{7} M_{\odot}$',color='k',linewidth=2)

    ax1.semilogy(lookback_time,(merg_8_01*Stellar_Mass_Function(8,1)),label='10$^{8}M_{\odot}$',color='b',linewidth=2)

    ax1.semilogy(lookback_time,(merg_9_01*Stellar_Mass_Function(9,1)),label='10$^{9}M_{\odot}$',color='c',linewidth=2)

    ax1.semilogy(lookback_time,(merg_10_01*Stellar_Mass_Function(10,1)),label='10$^{10}M_{\odot}$',color='g',linewidth=2)

    ax1.semilogy(lookback_time,(merg_11_01*Stellar_Mass_Function(11,1)),label='10$^{11}M_{\odot}$',color='r',linewidth=2)



    ax2 = ax1.twiny()
    zloc = [.1,.2,.5,1,2,3,4,5,7]

    tloc = [(z_to_Time(0)-z_to_Time(zl))/1e9 for zl in zloc]

    zstr = ['$'+str(i)+'$' for i in zloc]


    print lookback_time,z
    print tloc, zloc

#    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xlim(tloc[0],tloc[-1])
    ax2.set_xticks(tloc)
    ax2.set_xticklabels(zloc)
    fontsize=18
    ax2.set_xlabel(r'$z_{form}$',fontsize=fontsize)

    ax1.tick_params(labelsize=18)
    ax1.set_ylim(1e5,5e9)
    ax1.set_xlim(1.3,13)
#    plt.xlim(0,7)
    ax1.legend(bbox_to_anchor=(10, 1e8),bbox_transform=plt.gcf().transFigure)
#    ax1.legend(loc="best")
    ax1.legend(loc=8)

#    plt.ylim(1e7,1e12)
#    plt.xlim(0,7)
#    plt.legend(loc="best")
#    plt.xlabel(r"z$_{form}$",fontsize=18)
    ax1.set_xlabel(r"$t_{lookback}^{form}$ $(Gyr)$",fontsize=18)
    plt.ylabel(r"$log(M/M_{\odot}) Mpc^{-3}$",fontsize=18)
    plt.savefig("star_formation_evol_density_"+calib+"_lookback.pdf",dpi=400)

    plt.clf()

