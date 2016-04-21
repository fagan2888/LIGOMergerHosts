#!/bin/python


def z_to_Time(z,cosmo=None,verbose=False):
    if cosmo is None:
        #use Planck 2015 from last column (TT, TE, EE+lowP+lensing+ext) of Table 4 from http://arxiv.org/pdf/1502.01589v2.pdf
        from yt.utilities.cosmology import Cosmology
        h = 0.6774
        om = 0.3089
        ol = 0.6911
        if verbose: print "Assuming a Planck 2015 cosmology (H0 = {0}, Om0 = {1}, OL = {2})".format(h*100,om,ol)
        cosmo = Cosmology(hubble_constant=h,omega_matter=om,omega_lambda=ol)

    return cosmo.t_from_z(z).in_units('yr')
    #
    #
    # if type(H0) != type(YTQuantity(123.)):
    #     H0 = YTQuantity(H0,'km/s/Mpc')  #assume it's passed in in standard units
    # return ((2.0/H0)/(1+ (1+z)**2)).in_units('yr')

def time_to_z(t,cosmo=None,verbose=False):            # H0=YTQuantity(70.2,'km/s/Mpc')):
    from yt import YTQuantity,YTArray
    if cosmo is None:
        #use Planck 2015 from last column (TT, TE, EE+lowP+lensing+ext) of Table 4 from http://arxiv.org/pdf/1502.01589v2.pdf
        from yt.utilities.cosmology import Cosmology
        h = 0.6774
        om = 0.3089
        ol = 0.6911
        if verbose: print "Assuming a Planck 2015 cosmology (H0 = {0}, Om0 = {1}, OL = {2})".format(h*100,om,ol)
        cosmo = Cosmology(hubble_constant=h,omega_matter=om,omega_lambda=ol)

    if type(t) != type(YTQuantity(1,'Gyr')) and type(t) != type(YTArray([1,2,3],'Gyr')):
        #then I need to figure out units and wrap in a yt object
        if type(t) == type(1.23):   #single float
            if t < 15:  #assume Gyr
                t = YTArray(t,'Gyr')
                if verbose: print "Assuming time in Gyr"
            elif t < 1e11:  #assume yr
                t = YTArray(t,'yr')
                if verbose: print "Assuming time in yr"
            else:   #then it's probably in seconds
                t = YTArray(t,'s')
                if verbose: print "Assuming time in seconds"
        else:
            from numpy import array
            t = array(t)
            if (t < 15).all():
                t = YTArray(t,'Gyr')
                if verbose: print "Assuming time in Gyr"
            elif (t < 1e11).all():  #assume yr
                t = YTArray(t,'yr')
                if verbose: print "Assuming time in yr"
            else:   #then it's probably in seconds
                t = YTArray(t,'s')
                if verbose: print "Assuming time in seconds"

    return cosmo.z_from_t(t)

    # from numpy import sqrt
    # #inverting the fitting function form http://arxiv.org/pdf/gr-qc/0506079v2.pdf
    # if type(H0) != type(YTQuantity(123.,'AU')):
    #     H0 = YTQuantity(H0,'km/s/Mpc')  #assume it's passed in in standard units
    # if type(t) == type(np.array([1,2,3])):
    #     if (t<20).all():
    #         t = YTArray(t*1e9,'yr')
    #     else:
    #         t = YTArray(t,'yr')
    # elif type(t) == type(1.2):
    #     if t < 20:
    #         t = YTQuantity(t*1e9,'yr')
    #     else:
    #         t = YTQuantity(t,'yr')
    # #otherwise, assume it was passed in as a YTArray or YTQuantity with the right units
    # #
    # # if type(t) == type(YTArray([1,2,3],'yr')):
    # #     if (t<20).all():
    # #         t *= 1e9    #assume t was passed in in Gyr
    # # elif type(t) == type(YTQuantity(123.,'AU')):
    # #     if t < 20:  #assume Gyr:
    # #         t = YTQuantity(1e9*t,'yr')
    # #     else:
    # #         t = YTQuantity(t,'yr')
    # H0 = H0.in_units('1/yr')
    # return (-H0*t + sqrt(2.0*H0*t - H0*H0*t*t))/(H0*t)
