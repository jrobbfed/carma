"""Summary
"""
from shells import pv_slice
import numpy as np
import astropy.units as u

def ppv_model(outfile=None, dist=414*u.pc, pix_size=7.5*u.arcsec,
    vstep=0.099*u.km/u.s, acen=83.72707*u.deg, dcen=-5.07792*u.deg,
    thickness=0.0*u.pc, fwhm=0.0*u.km/u.s, beta=0.0, r=0.22*u.pc,
    dr=0.2*u.pc, vexp=2.2*u.km/u.s, depth_offset=0.*u.pc, 
    vel_offset=0.*u.km/u.s, v0=13.6*u.km/u.s, chan_pad=1.,
    ignore_cloud=True):
    """Summary
    
    Parameters
    ----------
    outfile : string, optional
        If set, write out fits ppv cube.
    dist : number or length-type astropy.Unit, optional
        Distance to shell. 
        Assumes pc if not specified. Orion A is at 414 pc (Menten + 2007)
    pix_size : number or angle-type astropy.Unit, optional
        Size of pixel in requested ppv cube. Assume arcseconds if not specified.
        Assumes 2*pix_size = gaussian beam FWHM
    vstep : number or speed-type astropy.Unit, optional
        Width of velocity channel in ppv cube. Assume km/s if not specified.
        NRO cube has 0.099 km/s channel width.
    acen : number or angle-type astropy.Unit, optional
        Right Ascension of shell. Assume degrees if not specified.
    dcen : number or angle-type astropy.Unit, optional
        Declination of shell. Assume degrees if not specified.
    thickness : number or length-type astropy.Unit, optional
        Thickness of cloud. Assume parsecs if not specified.
         Not implemented in this function, requires outside call to idl program.
    fwhm : number or speed-type astropy.Unit, optional
        Velocity fwhm of cloud. Assume km/s if not specified.
        Not implemented in this function, requires outside call to idl program.
    beta : number, optional
        Turbulent spectral index of cloud.
        Not implemented in this function, requires outside call to idl program.
    r : number or length-type astropy.Unit, optional
        Radius of shell. Assumes pc if not specified.
        Radius is measured from center to midpoint of shell rim.
    dr : number or length-type astropy.Unit, optional
        Thickness of shell rim. Assumes pc if not specified.
    vexp : number or speed-type astropy.Unit, optional
        Expansion velocity of shell. Assumes km/s if not specified.
    depth_offset : number or length-type astropy.Unit, optional
        Not implemented in this function, requires outside call to idl program.
    vel_offset : number or speed-type astropy.Unit, optional
        Not implemented in this function, requires outside call to idl program.
    v0 : number or speed-type astropy.Unit, optional
        Systemic velocity of shell. Assumes km/s if not specified.
    chan_pad : float, optional
        Factor to pad the ppv cube on either end in velocity channel space.
    ignore_cloud : bool, optional
        If True, ignore cloud completely. If not, will need to 
        call to idl program.
    """

   # Work in a pixel scale 2x finer than end result. 
   if ignore_cloud:
       scale = dist.to(u.pc) * pix_size.to(u.radian) / 2. # pc per pixel
       box_size = np.floor(4 * r / scale) # Number of pixels per side of ppp box.

       # den = np.zeros(3*[box_size])
       x, y, z = np.indices(3*[box_size]) - np.floor(box_size / 2) #Index w.r.t center of box.
       rr = np.sqrt(x ** 2. + y ** 2. + z ** 2.) * scale 

       #inside = rr < r - dr/2.
       #on = (rr >= r - dr/2.) & (rr < r + dr/2.) 
       #outside = rr >= r + dr/2.

       #ratio = np.sum(inside) / np.sum(on)
       den = (rr >= r - dr/2.) & (rr < r + dr/2.) #Density is 1 on rim of shell, 0 elsewhere.

       #Velocity field of bubble.
       #This is the z-component of velocity, assuming l.o.s. is in z direction. 
       #From similar triangles, z / r = v_z / v_r : v_r -> vexp, v_z -> vel, r*scale -> rr, z*scale -> z
       vel = on * vexp * z * scale / rr      
       vel[rr == 0] = 0 #Fix pole.

       vlo, vhi = -1. * chan_pad * vexp, chan_pad * vexp
       vcen = np.linspace(vlo, vhi, (vhi - vlo) / vstep)

       #Gridding PPV cube.
       ppv = ppp2ppv(den, vel, vcen)

def ppp2ppv(den, vel, vcen):
    """Summary
    
    Parameters
    ----------
    den : TYPE
        Description
    vel : TYPE
        Description
    vcen : array-like
        Array of velocity channel centers. Must be equally spaced.
    
    Returns
    -------
    TYPE
        Description
    """
    return ppv 

def analytic_profile(r, R, dR, radius_type="middle"):
    """Summary
    
    Parameters
    ----------
    r : TYPE
        Description
    R : TYPE
        Description
    dR : TYPE
        Description
    radius_type : str, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    
    Raises
    ------
    Exception
        Description
    """
    #From Beaumont & Williams (2010), assumes bare optically thin shell without cloud.
    #r: impact parameter from center of shell
    #R: radius of shell (measured at R_type of rim)
    #dR: thickness of shell
    #R_type: Where in the rim R is measured at. 
    
    if radius_type == "middle":
        R = R - dR / 2.
    elif radius_type == "inner":
        pass
    elif radius_type == "outer":
        R = R - dR
    else:
        raise Exception('radius_type must be one of "inner", "middle", or "outer".')

    r = np.asarray(r)
    profile = np.zeros_like(r)
    inside = (r < R)
    #print(inside)
    profile[inside] = \
                2*R * ((1 + dR/R)**2. - (r[inside]/R)**2.)**0.5 -\
                2*R * (1 - (r[inside]/R)**2.)**0.5

    on = (r >= R) & (r < R + dR)
    #print(on)
    profile[on] = \
                2*R * ((1 + dR/R)**2. - (r[on]/R)**2.)**0.5

    outside = (r >= R + dR)
    #print(outside)
    profile[outside] = 0

    return profile

def radial_profile(array, center=None, mode='average',
                   nbins=10, returnSEM=True, removeNaN=True):
    """Return a radial profile of `array`. The radial profile can
    be binned in `nbins` radial bins and then return a statistic given
    by `mode` such as the average in each radial bin. Or, if `mode` = 'alongx',
    returns the unbinned profile through the `center` of the array along the
    x-direction.
    
    Parameters
    ----------
    array : TYPE
        Description
    center : None, optional
        Description
    mode : str, optional
        Description
    nbins : int, optional
        Description
    returnSEM : bool, optional
        Description
    removeNaN : bool, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    
    Raises
    ------
    Exception
        Description
    """
    from scipy.stats import binned_statistic, sem
    import numpy as np
    y, x = np.indices(array.shape)
    if center is None:
        center = [(array.shape[1] - 1) / 2., (array.shape[0] - 1) / 2.]
    r = np.sqrt((center[0] - x) ** 2. + (center[1] - y) ** 2.)
    if removeNaN:
        r, array = r[np.isfinite(array)], array[np.isfinite(array)]
# Profile is the summed pixel values in each bin of radius.
# profile, rbins = np.histogram(r, weights=array, bins=nbins)
# print(r.max())
    if mode == 'average':
        # Divide by the number of pixels in each bin
        # profile = profile / np.histogram(r, bins=nbins)[0]
        profile, rbins, binnumber = binned_statistic(
            r, array, 'mean', bins=nbins)
        rbin_centers = (rbins[1:] + rbins[:-1]) / 2.
    elif mode == 'along x':
        print(np.shape(array), np.shape(x), np.shape(y))
        print(y, np.floor(center[1])) 
        profile = array[(y == np.floor(center[1])) & (x - np.floor(center[0]) >= 0)]
        rbin_centers = r[(y == np.floor(center[1])) & (x - np.floor(center[0]) >= 0)]
        returnSEM = False #wouldn't make sense in this context.
    else:
        raise Exception("mode {} not implemented".format(mode))

    if returnSEM:
        "Return the standard error on the mean of each bin."
        SEM = binned_statistic(
            r, array, lambda y: sem(y), bins=nbins)[0]
        return profile, SEM, rbin_centers
    else:
        return profile, rbin_centers
