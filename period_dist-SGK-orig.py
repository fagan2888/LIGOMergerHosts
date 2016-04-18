#!/bin/python


#ok, let's think this through:
#I have a distribution of periods -- let's think of it as f(<P)
#that's equivalent to a distribution of r0 from Kepler's Law:  a = ((G*(M1+M2)/4pi^2)*P^2)^(1/3)
#which means I have a distribution of infall times:   5./256.*(c**5 * a**4)/(G**3 * m**2 * mu)

#so, step 1 is to get f(<P) by intfrom astropy.cosmology import Planck13 as cosmo
egrating the Moe & di Stefano result:

#first though:  they are only sensitive to periods 2-20 days.  What coalescence times does that correspond to, for our system?
import yt.units as units
import yt.utilities.physical_constants as phys_const
from yt import YTArray,YTQuantity
import numpy as np
from numpy import pi,log10,logspace,array
from scipy.integrate import quad
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13 as cosmo


def Period_to_Semimajor(P,M1,M2):
    Mtot = M1 + M2
    return (((phys_const.G*Mtot)/(4.*pi*pi)) * P**2)**(1./3.)

def Semimajor_to_Period(A,M1,M2):
    Mtot = M1 + M2
    from numpy import sqrt,pi
    return sqrt((4*pi*pi) * A**3 / (phys_const.G*Mtot))

def Semimajor_to_Coalescence(a,M1,M2):
    mu = (M1*M2)/(M1+M2)
    Mtot = M1+M2
    return 5./256.*(phys_const.clight**5 * a**4)/(phys_const.G**3 * Mtot**2 * mu)

def Period_to_Coalescence(P,M1,M2):
    a = Period_to_Semimajor(P,M1,M2)
    return Semimajor_to_Coalescence(a,M1,M2)

th = cosmo.age(0)
thubble = YTQuantity(th*th.unit.in_units('yr'),'yr')

M1 = 30. * units.Msun
M2 = 30. * units.Msun
Mtot = M1 + M2
# mu = (M**2/Mtot)

Pmin = 2.*units.day
amin = Period_to_Semimajor(Pmin,M1,M2)
tmin = Semimajor_to_Coalescence(amin,M1,M2)
# amin = (((phys_const.G*Mtot)/(4.*pi*pi)) * Pmin**2)**(1./3.) #.in_units('kpc')
# tmin = 5./256.*(phys_const.clight**5 * amin**4)/(phys_const.G**3 * Mtot**2 * mu)

Pmax = 20.*units.day
amax = Period_to_Semimajor(Pmax,M1,M2)
tmax = Semimajor_to_Coalescence(amax,M1,M2)
# amax = (((phys_const.G*Mtot)/(4.*pi*pi)) * Pmax**2)**(1./3.) #.in_units('kpc')
# tmax = 5./256.*(phys_const.clight**5 * amax**4)/(phys_const.G**3 * Mtot**2 * mu)

print "Moe & Di Stefano are sensitive to coalescence times between {0} and {1}".format(tmin.in_units('yr'),tmax.in_units('yr'))
print "That is, log10 Tmin = {0}, log10 Tmax = {1}".format(log10(tmin.in_units('yr')),log10(tmax.in_units('yr')))
print "As a fraction of the age of the universe, that's {0} to {1}".format(tmin.in_units('yr')/thubble,tmax.in_units('yr')/thubble)



#now let's look at integrating the distribution:
def U_sub_P_dP(P,yp):
    """
    this assumes P is in days
    """
    return P**(yp-1)


#need to think about where I integrate from
#in principle, that lower limit should be set probably by the requirement that the typical surface density of a proto-stellar disk, times pi*r^2 is greater than the expecd mass
#from:  Dullemond thesis:
def sigma(r,runit,outunit="Msun/AU**2"):
    from yt import YTQuantity #,YTArray
    rinau = YTQuantity(r,runit).in_units('AU').item()
    sig_gcm2 = YTQuantity(1700*(rinau**-3./2.),'g/cm**2')
    return sig_gcm2.in_units(outunit).item()

def sigma_integrand(r,runit):
    from numpy import pi
    return 2*pi*r*sigma(r,runit)

constantsigma = True
#Rmin is set by the requirement that integral_0^Rmin (2*pi*R*dr*sigma) >= M
rbins = YTArray(logspace(-2,3,1e3),'AU')
if not constantsigma:
    print "Assuming surface density goes as 1700 g/cm^2 (R/AU)^(-3/2)"
    Mltr = [quad(sigma_integrand,0,r.in_units('AU').item(),args='AU')[0] for r in rbins]      #in Msun
else:
    print "Assuming a constant surface density of 1500 g/cm^2"
    Mltr = YTArray(np.pi*YTQuantity(1500,'g/cm**2')* (rbins**2),'Msun')

Mcut = 30.  #Msun
Rstar1 = (M1/YTQuantity(1,'Msun'))**(15./19.) * YTQuantity(1,'Rsun').in_units('AU')     #http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf
Rstar2 = (M2/YTQuantity(1,'Msun'))**(15./19.) * YTQuantity(1,'Rsun').in_units('AU')     #http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf


minR = max([rbins[Mltr>Mcut].min().in_units('AU'),(Rstar1+Rstar2).in_units('AU')])
if minR == Rstar1+Rstar2:
    print "Minimum separation is coming from the constraint that the two stars not touch"
else:
    print "Minimum separation is coming from the constraint that you have enough mass in a constant-density disk"
increase = 0.3
print "Now increasing minimum separation by {0}% due to mass loss during supernovae".format(increase*100.)
minR = minR*(1.+increase)

minP = Semimajor_to_Period(minR,M1,M2).in_units('day')
print "Minimum separation comes out to {0}, which corresponds to an orbital period of {1}".format(minR.in_units('AU'),minP.in_units('day'))
minCoal = Period_to_Coalescence(minP,M1,M2)
print "That corresponds to a coalesence time of {0}".format(minCoal.in_units('yr'))

#now I can integrate the period distribution from minR to R':
Pbins = YTArray(logspace(log10(minP.item()),log10(20),100),'day')        #days
Abins = Period_to_Semimajor(Pbins,M1,M2).in_units('AU')
Coalbins = Period_to_Coalescence(Pbins,M1,M2).in_units('yr')
# Abins = (((phys_const.G*Mtot)/(4.*pi*pi)) * Pbins**2)**(1./3.)

#TODO:  repeat this for the different powers that they get
NltP_pt9 = [quad(U_sub_P_dP,minP.in_units('day').item(),P.in_units('day').item(),args=-0.9)[0] for P in Pbins] #since I can convert P to coalescence, this is also N with coalescence less than some time
NltP_pt4 = [quad(U_sub_P_dP,minP.in_units('day').item(),P.in_units('day').item(),args=-0.4)[0] for P in Pbins]
NltP_pt3 = [quad(U_sub_P_dP,minP.in_units('day').item(),P.in_units('day').item(),args=-0.3)[0] for P in Pbins]
NltP_pt1 = [quad(U_sub_P_dP,minP.in_units('day').item(),P.in_units('day').item(),args=-0.1)[0] for P in Pbins]


# fig = plt.figure()
# plt.loglog(Coalbins,NltP,ls='-',color='k')
# plt.xlabel('$t_\mathrm{coalesce}\,(\mathrm{years})$')
# plt.ylabel('$f(\mathrm{coalescence\,time} < t_\mathrm{coalesce})$')
#
#
# plt.figure()
# plt.plot(Coalbins,NltP,ls='-',color='k')
# plt.xlabel('$t_\mathrm{coalesce}\,(\mathrm{years})$')
# plt.ylabel('$f(\mathrm{coalescence\,time} < t_\mathrm{coalesce})$')



def z_to_time(z,cosmology=None):
    """
    gives the time, in years, since the big bang given an astropy cosmology
    and a redshift.  if no cosmology passed in, assumes Planck13
    """
    if cosmology is None:
        from astropy.cosmology import Planck13 as cosmology
    from yt import YTArray
    t = cosmology.age(z)
    return YTArray(t*t.unit.in_units('yr'),'yr')


def time_to_z(age,cosmology=None,v=False):
    """
    returns the redshift of a given age using an astropy cosmology

    age is taken to be in Gyr if they're all less than 15, years otherwise
    (unless it's passed in as a YTQuantity/YTArray, then it's figured out)
    """
    from yt import YTArray,YTQuantity
    if cosmology is None:
        from astropy.cosmology import Planck13 as cosmology
    from astropy.cosmology import z_at_value
    import astropy.units as u
    import numpy as np

    gyrconv = False
    #numpy array?
    if type(age) == type(np.array([1,2,3])):
        if (age<15).all():
            gyrconv = True
            age = u.Quantity(age*1e9,u.yr)
        else:
            age = u.Quantity(age,u.yr)
    #single number?
    elif type(age) == type(1.2) or type(age) == type(1):
        if age < 15:
            gyrconv = True
            age = u.Quantity(age*1e9,u.yr)
        else:
            age = u.Quantity(age,u.yr)
    #yt quantity? convert it
    elif type(age) == type(YTQuantity(12e9,'yr')) or type(age) == type(YTArray([1.,2.])):
        age = u.Quantity(age.in_units('yr'),u.yr)
    #otherwise, gotta by an astropy quantity already
    else:
        assert type(age) == type(u.Quantity(13.6,u.yr))
    if v and gyrconv:
        print "Converted to Gyr"

    try:
        it = iter(age)
        z = []
        for ii in it:
            z.append(z_at_value(cosmology.age,ii))
        z = np.array(z)
    except TypeError, te:
        # age is not iterable
        z = z_at_value(cosmology.age,age)
    return z

"""
def z_to_time(z,H0=YTQuantity(70.2,'km/s/Mpc')):
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
"""





plt.figure()
ax = plt.gca()
ax.semilogx(Coalbins.in_units('yr'),array(NltP_pt4)/NltP_pt4[-1],ls='-',color='k',label=r'$\mathrm{MW:\,}y_\mathrm{p} = -0.4$')
ax.semilogx(Coalbins.in_units('yr'),array(NltP_pt3)/NltP_pt3[-1],ls='-',color='c',label=r'$\mathrm{OGLE-II,\,LMC:\,}y_\mathrm{p} = -0.3$')
ax.semilogx(Coalbins.in_units('yr'),array(NltP_pt1)/NltP_pt1[-1],ls='--',color='c',label=r'$\mathrm{OGLE-III,\,LMC:\,}y_\mathrm{p} = -0.1$')
ax.semilogx(Coalbins.in_units('yr'),array(NltP_pt9)/NltP_pt9[-1],ls='-',color='m',label=r'$\mathrm{OGLE-II,\,SMC:\,}y_\mathrm{p} = -0.9$')
ax.axvline(x=13.6e9,ls=':',color='grey')
ax.set_xlabel('$t_\mathrm{coalesce}\,(\mathrm{years})$')
ax.set_ylabel('$f(\mathrm{coalescence\,time} < t_\mathrm{coalesce})$')

plt.axvline(x=minCoal.in_units('yr').item(),ls='--',color='k')

leg = plt.legend(loc=2,fontsize=24)
# leg.get_frame().set_alpha(0)

ax2 = plt.twiny()
ax2.set_xscale('log')
ax2.set_xlim(ax.get_xlim())
ax2.set_xlabel(r'$\mathrm{Redshift\,at\,}t(z = 0.1) - t_\mathrm{coalesce}$')
t_zpt1 = z_to_time(0.1)
locs = YTArray([1e9,2.5e9,6e9,1.2e10],'yr')
formtime = t_zpt1 - locs        #time at which the system became a binary black hole
# msk = formtime > 0
redshifts = time_to_z(formtime)    #top labels
# bottomlocs = Coalbins[msk]
# toptickloc = bottomlocs[::20]
topticklab = ['{0:.2}'.format(z.item()) for z in redshifts]
ax2.set_xticks(locs)
ax2.set_xticklabels(topticklab)

plt.savefig('Coal_cumulative.png')


ax.set_xlim(1e9,thubble.in_units('yr').item())
ax.set_ylim(0,0.4)
ax2.set_xlim(ax.get_xlim())

locs = YTArray([1e9,2e9,3e9,5e9,7e9,1e10,1.2e10],'yr')
formtime = t_zpt1 - locs        #time at which the system became a binary black hole
# msk = formtime > 0
redshifts = time_to_z(formtime)    #top labels
# bottomlocs = Coalbins[msk]
# toptickloc = bottomlocs[::20]
topticklab = ['{0:.2}'.format(z.item()) for z in redshifts]
ax2.set_xticks(locs)
ax2.set_xticklabels(topticklab)

plt.savefig('Coal_cumulative2.png')

# plt.figure()
# plt.loglog(Pbins.in_units('day'),Coalbins.in_units('yr'),ls='-',color='k')
# plt.xlabel(r'$\mathrm{Orbital\,Period\,(days)}$')
# plt.ylabel(r'$\mathrm{coalescence\,time\,(years)}$')
#
# plt.figure()
# plt.loglog(Abins.in_units('AU'),Coalbins.in_units('yr'),ls='-',color='k')
# plt.xlabel(r'$\mathrm{Semimajor\,axis\,(AU)}$')
# plt.ylabel(r'$\mathrm{coalescence\,time\,(years)}$')

plt.figure()
ax = plt.gca()
plt.loglog(Pbins.in_units('day'),Coalbins.in_units('yr'),ls='-',color='k')
plt.xlabel(r'$\mathrm{Orbital\,Period\,(days)}$')
plt.ylabel(r'$\mathrm{coalescence\,time\,(years)}$')

ax.set_xlim(Pbins.in_units('day').min(),Pbins.in_units('day').max())
# ax.set_ylim(Coalbins.in_units('yr').min(),Coalbins.in_units('yr').max())
ax.set_ylim(1e9,1e12)

xtickloc = [2,3,4,5,8,10,20]
ax.set_xticks(xtickloc)
ax.set_xticklabels(xtickloc)

ax2 = plt.twiny()
ax2.set_xscale('log')
ax2.set_xlim(ax.get_xlim())
tickloc = Pbins[::len(Pbins)/5]
tickstr = ['{0:.5}'.format(Period_to_Semimajor(P,M1,M2).in_units('AU')) for P in tickloc]
ax2.set_xticks(tickloc)
ax2.set_xticklabels(tickstr)
ax2.set_xlabel(r'$\mathrm{Semimajor\,axis\,(AU)}$')

plt.savefig('Relations.png')

# # G = 4.301e-6 #kpc/Msun [km/s]^2
# t = 5./256.*(c**5 * r0_kpc**4)/(G**3 * m**2 * mu)
#
# plt.loglog(r0_kpc,t/3.15e7)
#
#
#
# yps = [-.4,-.3,-.1,-.9]
