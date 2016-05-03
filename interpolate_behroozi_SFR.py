#!/bin/python

import sys
from numpy import loadtxt,log10,arange,ones,exp,log,empty,unique,savetxt,column_stack
from pickle import dump
from mytools import griddata
from scipy.interpolate import interp2d
from mytools import abundance_match_behroozi_2012 as bAM

try:
    sfr_file,outbase = sys.argv[1:]
except Exception:
    print "usage:  python {0} <sfr_release.dat path> <outbase>".format(sys.argv[0].split('/')[-1])
    sys.exit(1)

#complication -- need to worry about the gray area of the plot
zp1,logMh,logSFR,logmedMstar = loadtxt(sfr_file,unpack=True)

Mhalo = 10**logMh
medMstar = 10**logmedMstar
# SFR = 10**logSFR
z = zp1 - 1

#to go from Median Mstar -> Avg Mstar, "simply" multiply by exp(0.5*(ln 10 * sigma)^2), where sigma is given by the following function:
def lognormalscatter(z):
    a = 1./(1+z)
    return 0.218 + (a-1)*(-0.023)

# avgMstar = empty(medMstar.shape[0])
# for ii in range(medMstar.shape[0]):
#     thisz = z[ii]
#     thissigma = lognormalscatter(thisz)
#     correction = exp(0.5*((log(10)*thissigma)**2))
#     avgMstar[ii] = medMstar[ii]*correction

#a-doy...don't need a loop
sigma = lognormalscatter(z)
correction = exp(0.5*((log(10)*sigma)**2))
avgMstar = medMstar*correction

logavgMstar = log10(avgMstar)

halomass,zmax = [],[]

for mh in unique(Mhalo):
    msk = Mhalo == mh
    thisz = z[msk]
    # thisSFR = SFR[msk]
    thislogSFR = logSFR[msk]

    maxz = max(thisz[thislogSFR!=-1000])

    halomass.append(mh)
    zmax.append(maxz)

savetxt(outbase+'Mhalo_MaxRedshift.txt',column_stack((halomass,zmax)),header="Mhalo\tMaxRedshift")

print "Making interpolation as a function of log Mhalo"
f_of_Mhalo = interp2d(z,logMh,logSFR,kind='linear')
print "Making interpolation as a function of log Mstar"
f_of_Mstar = interp2d(z,logavgMstar,logSFR,kind='linear')

outf = open(outbase+'logSFR_of_z_logMhalo.pkl','w')
dump(f_of_Mhalo,outf)
outf.close()
print "Saved interpolation as a function of logMhalo to {0}".format(outbase+'logSFR_of_z_logMhalo.pkl')


outf = open(outbase+'logSFR_of_z_logMstar.pkl','w')
dump(f_of_Mstar,outf)
outf.close()
print "Saved interpolation as a function of logMstar to {0}".format(outbase+'logSFR_of_z_logMstar.pkl')


#
#
#
#
#
#
# # zp1,logmh,logsfr = loadtxt(sfh_file,unpack=True)
# mh = 10**logmh
# # ms = bAM(mh,0)
# ms = bAM(mh,0,alpha=-1.92)      #use the SGK AM at low stellar masses
# logms = log10(ms)
#
# zp1_grid,logms_grid,logsfr_grid = griddata(zp1,logms,logsfr)
#
# zp1_grid_mh,logmh_grid,logsfr_grid_mh = griddata(zp1,logmh,logsfr)
#
# f = interp2d(zp1,logms,logsfr,kind='linear')
#
# outf = open(outfile,'w')
# dump(f,outf)
# outf.close()
# print "Saved interpolation to {0}".format(outfile)

import matplotlib.pyplot as plt

testlogMs = arange(4.6,11.4,.25)
testlogMh = arange(9,15.6,0.25)
testz = arange(0,8,.5)

testlogSFR_mh = f_of_Mhalo(testz,testlogMh)
testlogSFR_ms = f_of_Mstar(testz,testlogMs)

zgrid,mhgrid,logsfr_mhgrid = griddata(z,logMh,logSFR)

plt.pcolormesh(zgrid,mhgrid,logsfr_mhgrid,vmin=-3,vmax=3)
for ii in range(len(testlogMh)):
    plt.scatter(testz,ones(testz.shape[0])*testlogMh[ii],c=testlogSFR_mh[ii],s=50,vmin=-3,vmax=3)

plt.xlabel(r'$z$')
plt.ylabel(r'$\log M_\mathrm{halo} (M_\odot)$')
plt.xlim(0,8)
plt.ylim(9,15.5)
cbar = plt.colorbar()
cbar.set_label(r'$\log \mathrm{SFR}$')
plt.savefig(outbase+'interp_check_Mh.png')

plt.close('all')
#
# zgrid,msgrid,logsfr_msgrid = griddata(z,logavgMstar,logSFR)
#
# plt.pcolormesh(zgrid,msgrid,logsfr_msgrid,vmin=-3,vmax=3)
# for ii in range(len(testlogMs)):
#     plt.scatter(testz,ones(testz.shape[0])*testlogMs[ii],c=testlogSFR_ms[ii],s=50,vmin=-3,vmax=3)
#
# plt.xlabel(r'$z$')
# plt.ylabel(r'$\log M_\star (M_\odot)$')
# plt.xlim(0,8)
# plt.ylim(4.5,11.5)
# cbar = plt.colorbar()
# cbar.set_label(r'$\log \mathrm{SFR}$')
# plt.savefig(outbase+'interp_check_Ms.png')


print "Now filling in the area of the plot where halos don't exist..."
#now do i tall again, but fill in the -1000s with the same value as appears before
# logSFR_filled = array(logSFR,copy=True)
for mh in unique(Mhalo):
    msk = Mhalo == mh
    thisz = z[msk]
    thislogSFR = logSFR[msk]

    maxz = max(thisz[thislogSFR!=-1000])
    sfr_at_maxz = thislogSFR[thisz==maxz]
    logSFR[msk&(z>maxz)] = sfr_at_maxz

print "Making interpolation as a function of log Mhalo"
f_of_Mhalo = interp2d(z,logMh,logSFR,kind='linear')
print "Making interpolation as a function of log Mstar"
f_of_Mstar = interp2d(z,logavgMstar,logSFR,kind='linear')

outf = open(outbase+'logSFR_of_z_logMhalo-filled.pkl','w')
dump(f_of_Mhalo,outf)
outf.close()
print "Saved interpolation as a function of logMhalo to {0}".format(outbase+'logSFR_of_z_logMhalo-filled.pkl')


outf = open(outbase+'logSFR_of_z_logMstar-filled.pkl','w')
dump(f_of_Mstar,outf)
outf.close()
print "Saved interpolation as a function of logMstar to {0}".format(outbase+'logSFR_of_z_logMstar-filled.pkl')

testlogSFR_mh = f_of_Mhalo(testz,testlogMh)
testlogSFR_ms = f_of_Mstar(testz,testlogMs)

zgrid,mhgrid,logsfr_mhgrid = griddata(z,logMh,logSFR)

plt.pcolormesh(zgrid,mhgrid,logsfr_mhgrid,vmin=-3,vmax=3)
for ii in range(len(testlogMh)):
    plt.scatter(testz,ones(testz.shape[0])*testlogMh[ii],c=testlogSFR_mh[ii],s=50,vmin=-3,vmax=3)

plt.xlabel(r'$z$')
plt.ylabel(r'$\log M_\mathrm{halo} (M_\odot)$')
plt.xlim(0,8)
plt.ylim(9,15.5)
cbar = plt.colorbar()
cbar.set_label(r'$\log \mathrm{SFR}$')
plt.savefig(outbase+'interp_check_Mh-filled.png')

plt.close('all')

#
# zgrid,msgrid,logsfr_msgrid = griddata(z,logavgMstar,logSFR)
#
# plt.pcolormesh(zgrid,msgrid,logsfr_msgrid,vmin=-3,vmax=3)
# for ii in range(len(testlogMs)):
#     plt.scatter(testz,ones(testz.shape[0])*testlogMs[ii],c=testlogSFR_ms[ii],s=50,vmin=-3,vmax=3)
#
# plt.xlabel(r'$z$')
# plt.ylabel(r'$\log M_\star (M_\odot)$')
# plt.xlim(0,8)
# plt.ylim(4.5,11.5)
# cbar = plt.colorbar()
# cbar.set_label(r'$\log \mathrm{SFR}$')
# plt.savefig(outbase+'interp_check_Ms-filled.png')
