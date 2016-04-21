#!/bin/python

import sys
from numpy import loadtxt,log10,arange,ones
from pickle import dump
from mytools import griddata
from scipy.interpolate import interp2d
from mytools import abundance_match_behroozi_2012 as bAM

try:
    sfh_file,outfile = sys.argv[1:]
except Exception:
    print "usage:  python {0} <sfh_release.dat path> <outfile>".format(sys.argv[0].split('/')[-1])
    sys.exit(1)

zp1,logmh,logsfr = loadtxt(sfh_file,unpack=True)
mh = 10**logmh
ms = bAM(mh,0,alpha=-1.92)      #use the SGK AM at low stellar masses
logms = log10(ms)

zp1_grid,logms_grid,logsfr_grid = griddata(zp1,logms,logsfr)

zp1_grid_mh,logmh_grid,logsfr_grid_mh = griddata(zp1,logmh,logsfr)

f = interp2d(zp1,logms,logsfr,kind='linear')

outf = open(outfile,'w')
dump(f,outf)
outf.close()
print "Saved interpolation to {0}".format(outfile)

import matplotlib.pyplot as plt

testlogMs = arange(4.6,11.4,.25)
test1pz = arange(1,9.5,.5)
testlogSFR = f(test1pz,testlogMs)
plt.pcolormesh(zp1_grid,logms_grid,logsfr_grid,vmin=-3,vmax=3)
for ii in range(len(testlogMs)):
    plt.scatter(test1pz,ones(test1pz.shape[0])*testlogMs[ii],c=testlogSFR[ii],s=50,vmin=-3,vmax=3)

plt.xlabel(r'$1+z$')
plt.ylabel(r'$\log M_\star(z=0) (M_\odot)$')
plt.xlim(.9,9.1)
plt.ylim(4.5,11.5)
cbar = plt.colorbar()
cbar.set_label(r'$\log \mathrm{SFR}$')
plt.savefig(outfile[:-len('.pkl')]+'.png')
