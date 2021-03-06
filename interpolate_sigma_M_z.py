#!/usr/bin/env python

import sys

if len(sys.argv) == 2:
    outname = sys.argv[1]
else:
    print "usage:  python {0} <outname>".format(sys.argv[0].split('/')[-1])
    sys.exit(1)

from pickle import dump
import numpy as np
from scipy.interpolate import interp2d,InterpolatedUnivariateSpline
from hmf import hmf
# from astropy.cosmology import FlatLambdaCDM


zbins = np.arange(0,8,.1)
mmin = 8
mmax = 15
dlog10m = 0.1   #steps of 0.1 in logM

siglist = []
mlist = []
zlist = []

# astropy_cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
# hmf_cosmo = cosmo.Cosmology(cosmo_model=astropy_cosmo)

print "Grabbing sigma(M,z)"
cosmodict = {'Om0':0.3,'H0':70.}

mf = hmf.MassFunction(z=0,Mmin=mmin,Mmax=mmax,dlog10m=dlog10m,cosmo_params=cosmodict)
for z in zbins:
    mf.update(z=z)
    siglist.append(mf.sigma)
    mlist.append(mf.M)
    zlist.append(np.ones(mf.M.shape[0])*z)

siglist = np.array(siglist).flatten()
zlist = np.array(zlist).flatten()
mlist = np.array(mlist).flatten()

logM = np.log10(mlist)
logsig = np.log10(siglist)

print "Making interpolation"
f = interp2d(logM,zlist,logsig,kind='linear')
print "Saving interpolation"
outf = open(outname,'w')
dump(f,outf)
outf.close()
print "Wrote interpolation to "+outname


import matplotlib.pyplot as plt
from mytools import griddata

testz = np.arange(0,8,.5)
testm = np.logspace(8,15,20)

# testlogsig = f(np.log10(testm),testz)

logmgrid,zgrid,logsgrid = griddata(logM,zlist,logsig)

fig = plt.figure()
ax = plt.gca()

# ax.set_xscale('log')


vmin = np.min(logsig)
vmax = np.max(logsig)
im = plt.pcolormesh(logmgrid,zgrid,logsgrid,vmin=vmin,vmax=vmax)

for ii in range(len(testz)):
    plt.scatter(np.log10(testm),np.ones(testm.shape[0])*testz[ii],c=f(np.log10(testm),testz[ii]),s=50,vmin=vmin,vmax=vmax)

# cbar = plt.colorbar(im)

plt.xlabel(r'$\mathrm{mass\,(M_\odot)}$')
plt.ylabel(r'$\mathrm{redshift}$')

plt.savefig(outname[:-4]+'2dinterp_test.png')

plt.close('all')

fig = plt.figure()
from mytools import gencolors
colors = gencolors(len([0,1,2,3,4,5,6,7,8]))

ii = 0
for z in [0,1,2,3,4,5,6,7,8]:
    myz = zbins[np.argmin(np.abs(z-zbins))]
    msk = zlist == myz

    plt.plot(logM[msk],logsig[msk],ls='-',color=colors[ii],label=r'$z = '+str(myz)+'$')
    plt.plot(np.log10(testm),f(np.log10(testm),myz),ls='--',color=colors[ii],label='_nolegend_')
    ii+=1

leg = plt.legend(loc=3,ncol=2,fontsize=24)
leg.get_frame().set_alpha(0)

plt.xlabel(r'$log Mhalo$')
plt.ylabel(r'$\sigma$')

plt.savefig(outname[:-4]+'1dinterp_test.png')
plt.close('all')

print "Now doing growth factor"
#now do growth factor:
growth = []
for z in zbins:
    mf.update(z=z)
    growth.append(mf.growth_factor)

growth = np.array(growth)
#now interpolate:
g_of_z = InterpolatedUnivariateSpline(zbins,growth)

outf = open(outname[:-4]+'growth_factor_of_z.pkl','w')
dump(g_of_z,outf)
outf.close()

fig = plt.figure()
plt.plot(zbins,growth,'k-')
plt.plot(zbins,g_of_z(zbins),'r--')
# plt.plot
plt.xlabel('redshift')
plt.ylabel('growth factor')
plt.axhline(y=1)

plt.savefig(outname[:-4]+'growth_plot.pdf')
