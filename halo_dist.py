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
import corner

#TO DO:
#-create star formation rate at given halo mass and redshift
#-include contributions from subhalos to low Z gas
#-substract non merging systems
#-make fancy plot

dt=1e8
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
calib="KK04"
#number subhalo mass bins
nmassbins=20
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
    star_density=stars_per_gal*Stellar_Mass_Function(lMgal,ldM)
    return star_density

    

def make_plot():
    testing_SMF=np.arange(7,11.25,ldM)
    SMF=[]

    #print testing_SMF
    for i in testing_SMF :
        SMF.append(Stellar_Mass_Function(i,ldM))
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
    lMhalo=Gal_to_Halo_mass(lMgal)

    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])

        #loop over all subhalo masses
        lsubhM_min=-4.95
        lsubhM_max=-.05

#        mass_ratios=np.linspace(lsubhM_min,lsubhM_max,nmassbins)
        size_bins=(lsubhM_min-lsubhM_max)/nmassbins
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i])
        lowZ_stars_in_halo=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):

            #get corresponding galaxy mass through abundance matching 
            lMsubgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt
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


def get_distrib_lowZ_stars(lMgal,Zcut):
    #provides low Z stellar mass as a function of redhsift of formation
    t_min=z_to_Time(z_max)
    t_max=z_to_Time(z_min)
    time_bins=np.linspace(t_min,t_max,(t_max-t_min)/dt)
    z_bins=np.zeros(len(time_bins))
    lowZ_stars=np.zeros(len(time_bins))
    lMhalo=Gal_to_Halo_mass(lMgal)
    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])
        #loop over all subhalo masses
        lsubhM_min=-4.95
        lsubhM_max=-.05
        size_bins=(lsubhM_min-lsubhM_max)/nmassbins
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i])
        lowZ_stars_in_halo=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):
            #get corresponding galaxy mass through abundance matching 
            lMgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt
            #multiply by fraction in mass bin
            stellar_mass_formed=stellar_mass_formed*mass_frac[m]
            #get amount of lowZ stuff
            lowZ_stars_in_halo[m]=stellar_mass_formed*Frac_LowZ_Stars_Avail(lMgal,z_bins[i],Zcut,calib)
        lowZ_stars[i]=np.sum(lowZ_stars_in_halo)
#        print lowZ_stars
    return lowZ_stars,z_bins









def get_distrib_stars(lMhalo,Zcut):
    z_bins=np.zeros(len(time_bins))
    all_stars=np.zeros(len(time_bins))
    all_mass=0

    for i in range (0,len(time_bins)):
        z_bins[i]=time_to_z(time_bins[i])
        #loop over all subhalo masses
        lsubhM_min=-4.95
        lsubhM_max=-.05

        #mass_ratios=np.linspace(lsubhM_min,lsubhM_max,nmassbins)
        size_bins=(lsubhM_min-lsubhM_max)/nmassbins
        mass_ratios,mass_frac=EPS_is_awesome(lMhalo,z_bins[i])
        stellar_mass_formed=np.zeros(len(mass_ratios))
        for m in range(0,len(mass_ratios)):
            #get corresponding galaxy mass through abundance matching 
            lMgal=np.log10(abundance_match_behroozi_2012(10**(lMhalo+mass_ratios[m]),z_bins[i]))
            #get star formation at that redshift,for each subhalo
            stellar_mass_formed[m]= SFR_lMhalo((lMhalo+mass_ratios[m]),z_bins[i])*dt
            #multiply by fraction in mass bin
            stellar_mass_formed[m]=stellar_mass_formed[m]*mass_frac[m]
        all_stars[i]=np.sum(stellar_mass_formed)
         
    return all_stars,z_bins

def get_all_mergers(Zcut):

    #make data for Shea
    size=ngalbins*ntimebins
    print size,'size'
    all_mgal=np.zeros(size)
    all_zform=np.zeros(size)
    all_rates=np.zeros(size)
    all_stars=np.empty((ngalbins,ntimebins))
    print lMgalbins,'galbins'
    print len(all_mgal)
#    all_stars=np.zeros(ngalbins,ntimebins)
    for m in range(0,ngalbins):
        # for t in ntimebins:
#        z[m]=time_to_z(time_bins[t])
        lMhalo=Gal_to_Halo_mass(lMgalbins[m])
        stars,z_bins=get_distrib_stars(lMhalo,Zcut)
        all_stars[m,:]=stars
        
    for m in range(0,ngalbins):
        for t in range(0,ntimebins):
            r=m*ntimebins+t
            all_mgal[r]=lMgalbins[m]
            all_zform[r]=z_bins[t]
            all_rates[r]=all_stars[m,t]
    np.save("lMgal_zform_rates_Z"+str(Zcut),(all_mgal,all_zform,all_rates))
            
    return all_stars,z_bins

    



def make_merger_formation_plot():
#    stellar_dens=np.arange(7,11.25,1)
    #merger_dens_01=[]
    #merger_dens_05=[]
    #merger_dens_001=[]


    merg_7_03,z=get_distrib(7.,.3) 
    merg_8_03,z=get_distrib(8.,.3)
    merg_9_03,z=get_distrib(9.,.3)
    merg_10_03,z=get_distrib(10.,.3)
    merg_11_03,z=get_distrib(11.,.3)


    merg_7_01,z=get_distrib(7.,.1)
    merg_8_01,z=get_distrib(8.,.1)
    merg_9_01,z=get_distrib(9.,.1)
    merg_10_01,z=get_distrib(10.,.1)
    merg_11_01,z=get_distrib(11.,.1)


    merg_7_001,z=get_distrib(7.,.01)
    merg_8_001,z=get_distrib(8.,.01)
    merg_9_001,z=get_distrib(9.,.01)
    merg_10_001,z=get_distrib(10.,.01)
    merg_11_001,z=get_distrib(11.,.01)


    plt.clf()
    #divide by ten to get rates in events/yr/Gpc


    plt.figure(figsize=(7.5,6))
    plt.semilogy(z,merg_7_01*Stellar_Mass_Function(7,1)/10,label='10$^{7}M_{\odot}$',color='k',linewidth=2)
    plt.semilogy(z,merg_7_001*Stellar_Mass_Function(7,1)/10,':',color='k',linewidth=2)
    plt.semilogy(z,merg_7_03*Stellar_Mass_Function(7,1)/10,'--',color='k',linewidth=2)

    plt.semilogy(z,merg_8_01*Stellar_Mass_Function(8,1)/10,label='10$^{8}M_{\odot}$',color='b',linewidth=2)
    plt.semilogy(z,merg_8_001*Stellar_Mass_Function(8,1)/10,':',color='b',linewidth=2)
    plt.semilogy(z,merg_8_03*Stellar_Mass_Function(8,1)/10,'--',color='b',linewidth=2)

    plt.semilogy(z,merg_9_01*Stellar_Mass_Function(9,1)/10,label='10$^{9}M_{\odot}$',color='c',linewidth=2)
    plt.semilogy(z,merg_9_001*Stellar_Mass_Function(9,1)/10,':',color='c',linewidth=2)
    plt.semilogy(z,merg_9_03*Stellar_Mass_Function(9,1)/10,'--',color='c',linewidth=2)
    
    plt.semilogy(z,merg_10_01*Stellar_Mass_Function(10,1)/10,label='10$^{10}M_{\odot}$',color='g',linewidth=2)
    plt.semilogy(z,merg_10_001*Stellar_Mass_Function(10,1)/10,':',color='g',linewidth=2)
    plt.semilogy(z,merg_10_03*Stellar_Mass_Function(10,1)/10,'--',color='g',linewidth=2)
    

    plt.semilogy(z,merg_11_01*Stellar_Mass_Function(11,1)/10,label='10$^{11}M_{\odot}$',color='r',linewidth=2)
    plt.semilogy(z,merg_11_001*Stellar_Mass_Function(11,1)/10,':',color='r',linewidth=2)
    plt.semilogy(z,merg_11_03*Stellar_Mass_Function(11,1)/10,'--',color='r',linewidth=2)



    plt.tick_params(labelsize=18)
    plt.ylim(1e0,1e4)
    plt.xlim(0,7)
    plt.legend(loc="best")
    plt.xlabel(r"z$_{form}$",fontsize=18)
    plt.ylabel(r"$BH$ $merger$ $rate$ $(yr^{-1} Gpc^{-3})$",fontsize=18)
    plt.savefig("merger_formation_density_withsubhalos"+calib+".pdf",dpi=400)




def make_lowZstar_formation_plot():
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


    
#    print merg_7_01
#    print z

    plt.clf()

    plt.figure(figsize=(7.5,6))

    # plt.semilogy(z,merg_7_01,label='10$^{15} M_{\odot}$',color='k',linewidth=2)
    # plt.semilogy(z,merg_7_001,':',color='k',linewidth=2)
    # plt.semilogy(z,merg_7_03,'--',color='k',linewidth=2)

    # plt.semilogy(z,merg_8_01,label='10$^{11}M_{\odot}$',color='b',linewidth=2)
    # plt.semilogy(z,merg_8_001,':',color='b',linewidth=2)
    # plt.semilogy(z,merg_8_03,'--',color='b',linewidth=2)

    # plt.semilogy(z,merg_9_01,label='10$^{12}M_{\odot}$',color='c',linewidth=2)
    # plt.semilogy(z,merg_9_001,':',color='c',linewidth=2)
    # plt.semilogy(z,merg_9_03,'--',color='c',linewidth=2)

    # plt.semilogy(z,merg_10_01,label='10$^{13}M_{\odot}$',color='g',linewidth=2)
    # plt.semilogy(z,merg_10_001,':',color='g',linewidth=2)
    # plt.semilogy(z,merg_10_03,'--',color='g',linewidth=2)

    # plt.semilogy(z,merg_11_01,label='10$^{14}M_{\odot}$',color='r',linewidth=2)
    # plt.semilogy(z,merg_11_001,':',color='r',linewidth=2)
    # plt.semilogy(z,merg_11_03,'--',color='r',linewidth=2)


    plt.semilogy(z,merg_7_01*Stellar_Mass_Function(7,1),label='10$^{7} M_{\odot}$',color='k',linewidth=2)
    plt.semilogy(z,merg_7_001*Stellar_Mass_Function(7,1),':',color='k',linewidth=2)
    plt.semilogy(z,merg_7_03*Stellar_Mass_Function(7,1),'--',color='k',linewidth=2)

    plt.semilogy(z,merg_8_01*Stellar_Mass_Function(8,1),label='10$^{8}M_{\odot}$',color='b',linewidth=2)
    plt.semilogy(z,merg_8_001*Stellar_Mass_Function(8,1),':',color='b',linewidth=2)
    plt.semilogy(z,merg_8_03*Stellar_Mass_Function(8,1),'--',color='b',linewidth=2)

    plt.semilogy(z,merg_9_01*Stellar_Mass_Function(9,1),label='10$^{9}M_{\odot}$',color='c',linewidth=2)
    plt.semilogy(z,merg_9_001*Stellar_Mass_Function(9,1),':',color='c',linewidth=2)
    plt.semilogy(z,merg_9_03*Stellar_Mass_Function(9,1),'--',color='c',linewidth=2)

    plt.semilogy(z,merg_10_01*Stellar_Mass_Function(10,1),label='10$^{10}M_{\odot}$',color='g',linewidth=2)
    plt.semilogy(z,merg_10_001*Stellar_Mass_Function(10,1),':',color='g',linewidth=2)
    plt.semilogy(z,merg_10_03*Stellar_Mass_Function(10,1),'--',color='g',linewidth=2)

    plt.semilogy(z,merg_11_01*Stellar_Mass_Function(11,1),label='10$^{11}M_{\odot}$',color='r',linewidth=2)
    plt.semilogy(z,merg_11_001*Stellar_Mass_Function(11,1),':',color='r',linewidth=2)
    plt.semilogy(z,merg_11_03*Stellar_Mass_Function(11,1),'--',color='r',linewidth=2)




    
    plt.tick_params(labelsize=18)
    plt.ylim(1e5,5e9)
    plt.xlim(0,7)
    plt.legend(loc="best")
    plt.xlabel(r"z$_{form}$",fontsize=18)
    plt.ylabel(r"$\log(M_{\odot})$ $ Mpc^{-3}$",fontsize=18)
    plt.savefig("lowZ_star_formation_evol_density"+calib+".pdf",dpi=400)

    plt.clf()


def make_star_formation_plot():
#    stellar_dens=np.arange(7,11.25,1)
    #merger_dens_01=[]
    #merger_dens_05=[]
    #merger_dens_001=[]

    merg_7_03,z=get_distrib_stars(15.,.3)
    merg_8_03,z=get_distrib_stars(11.,.3)
    merg_9_03,z=get_distrib_stars(12.,.3)
    merg_10_03,z=get_distrib_stars(13.,.3)
    merg_11_03,z=get_distrib_stars(14.,.3)

    merg_7_01,z=get_distrib_stars(15.,.1)
    merg_8_01,z=get_distrib_stars(11.,.1)
    merg_9_01,z=get_distrib_stars(12.,.1)
    merg_10_01,z=get_distrib_stars(13.,.1)
    merg_11_01,z=get_distrib_stars(14.,.1)

    merg_7_001,z=get_distrib_stars(15.,.01)
    merg_8_001,z=get_distrib_stars(11.,.01)
    merg_9_001,z=get_distrib_stars(12.,.01)
    merg_10_001,z=get_distrib_stars(13.,.01)
    merg_11_001,z=get_distrib_stars(14.,.01)


    
#    print merg_7_01
#    print z

    plt.clf()

    plt.figure(figsize=(7,6))

    plt.semilogy(z,merg_7_01,label='10$^{15} M_{\odot}$',color='k',linewidth=2)
    plt.semilogy(z,merg_7_001,':',color='k',linewidth=2)
    plt.semilogy(z,merg_7_03,'--',color='k',linewidth=2)

    plt.semilogy(z,merg_8_01,label='10$^{11}M_{\odot}$',color='b',linewidth=2)
    plt.semilogy(z,merg_8_001,':',color='b',linewidth=2)
    plt.semilogy(z,merg_8_03,'--',color='b',linewidth=2)

    plt.semilogy(z,merg_9_01,label='10$^{12}M_{\odot}$',color='c',linewidth=2)
    plt.semilogy(z,merg_9_001,':',color='c',linewidth=2)
    plt.semilogy(z,merg_9_03,'--',color='c',linewidth=2)

    plt.semilogy(z,merg_10_01,label='10$^{13}M_{\odot}$',color='g',linewidth=2)
    plt.semilogy(z,merg_10_001,':',color='g',linewidth=2)
    plt.semilogy(z,merg_10_03,'--',color='g',linewidth=2)

    plt.semilogy(z,merg_11_01,label='10$^{14}M_{\odot}$',color='r',linewidth=2)
#    plt.semilogy(z,merg_11_001,':',color='r',linewidth=2)
#    plt.semilogy(z,merg_11_03,'--',color='r',linewidth=2)


    # plt.semilogy(z,merg_7_01*Halo_Mass_Function(10,1),label='10$^{10} M_{\odot}$',color='k',linewidth=2)
    # plt.semilogy(z,merg_7_001*Halo_Mass_Function(10,1),':',color='k',linewidth=2)
    # plt.semilogy(z,merg_7_03*Halo_Mass_Function(10,1),'--',color='k',linewidth=2)

    # plt.semilogy(z,merg_8_01*Halo_Mass_Function(11,1),label='10$^{11}M_{\odot}$',color='b',linewidth=2)
    # plt.semilogy(z,merg_8_001*Halo_Mass_Function(11,1),':',color='b',linewidth=2)
    # plt.semilogy(z,merg_8_03*Halo_Mass_Function(11,1),'--',color='b',linewidth=2)

    # plt.semilogy(z,merg_9_01*Halo_Mass_Function(12,1),label='10$^{12}M_{\odot}$',color='c',linewidth=2)
    # plt.semilogy(z,merg_9_001*Halo_Mass_Function(12,1),':',color='c',linewidth=2)
    # plt.semilogy(z,merg_9_03*Halo_Mass_Function(12,1),'--',color='c',linewidth=2)

    # plt.semilogy(z,merg_10_01*Halo_Mass_Function(13,1),label='10$^{13}M_{\odot}$',color='g',linewidth=2)
    # plt.semilogy(z,merg_10_001*Halo_Mass_Function(13,1),':',color='g',linewidth=2)
    # plt.semilogy(z,merg_10_03*Halo_Mass_Function(13,1),'--',color='g',linewidth=2)

    # plt.semilogy(z,merg_11_01*Halo_Mass_Function(14,1),label='10$^{14}M_{\odot}$',color='r',linewidth=2)
    # plt.semilogy(z,merg_11_001*Halo_Mass_Function(14,1),':',color='r',linewidth=2)
    # plt.semilogy(z,merg_11_03*Halo_Mass_Function(14,1),'--',color='r',linewidth=2)




    
    plt.tick_params(labelsize=18)
    plt.ylim(1e7,1e12)
    plt.xlim(0,7)
    plt.legend(loc="best")
    plt.xlabel(r"z$_{form}$",fontsize=18)
    plt.ylabel(r"$log(M/M_{\odot}) Mpc^{-3}$",fontsize=18)
    plt.savefig("star_formation_evol_number"+calib+".pdf",dpi=400)

    plt.clf()
