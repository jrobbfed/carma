"""Summary
"""
#from shells import pv_slice
import numpy as np
import astropy.units as u
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
from astropy.io import fits

def ppv_model(outfile=None, dist=414*u.pc, pix_size=7.5*u.arcsec,
    vstep=0.099*u.km/u.s, acen=83.72707*u.deg, dcen=-5.07792*u.deg,
    thickness=0.0*u.pc, fwhm=0.0*u.km/u.s, beta=0.0, r=0.22*u.pc,
    dr=0.2*u.pc, vexp=2.2*u.km/u.s, depth_offset=0.*u.pc, 
    vel_offset=0.*u.km/u.s, v0=13.6*u.km/u.s, chan_pad=1.,
    ignore_cloud=True, write_fits=True, return_hdu=True,
    return_ppp=False, downsample=True, smooth=True):
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
    write_fits : bool, optional
        Description
    return_hdu : bool, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """

   # Work in a pixel scale 2x finer than end result. 
    if ignore_cloud:
        scale = dist.to(u.pc) * pix_size.to(u.radian).value / 2. # pc per pixel
        box_size = np.floor(4 * r / scale) # Number of pixels per side of ppp box.

        # den = np.zeros(3*[box_size])
        # x, y, z = np.indices(3*[box_size]) - np.floor(box_size / 2) #Index w.r.t center of box.
        x, y, z = np.indices(3*[box_size]) - np.floor(box_size / 2) #Index w.r.t center of box.

        rr = np.sqrt(x ** 2. + y ** 2. + z ** 2.) * scale 

        #inside = rr < r - dr/2.
        #on = (rr >= r - dr/2.) & (rr < r + dr/2.) 
        #outside = rr >= r + dr/2.

        #ratio = np.sum(inside) / np.sum(on)
        den = ((rr >= r - dr/2.) & (rr < r + dr/2.)).astype(int) #Density is 1 on rim of shell, 0 elsewhere.

        #Velocity field of bubble.
        #This is the z-component of velocity, assuming l.o.s. is in z direction. 
        #From similar triangles, z / r = v_z / v_r : v_r -> vexp, v_z -> vel, r*scale -> rr, z*scale -> z
        vel = den * vexp * z * scale / rr      
        vel[rr == 0] = 0 #Fix pole.

        vlo, vhi = -1. * chan_pad * vexp, chan_pad * vexp
        vcen = np.linspace(vlo, vhi, (vhi - vlo) / vstep)

        #Gridding PPV cube.
        ppv = ppp2ppv(den, vel.value, vcen.value)
        if downsample:
            ppv = congrid(ppv, (ppv.shape[0]/2, ppv.shape[1]/2, ppv.shape[2]))
        if smooth:
            gauss = Gaussian2DKernel(stddev=2) #2 pixels

            for i in range(ppv.shape[2]):
                 ppv[:,:,i] = convolve_fft(ppv[:,:,i], gauss, normalize_kernel=True)
        
        ppv = np.swapaxes(ppv, 0, 2)
 
        if return_hdu:
            h = fits.Header()
            h["CTYPE1"] = "RA---TAN"
            h["CRPIX1"] = ppv.shape[0]/2.
            h["CRVAL1"] = (u.Quantity(acen, u.deg).value, 'DEGREES')
            h["CDELT1"] = (u.Quantity(pix_size, u.arcsec).to(u.deg).value, 'DEGREES')

            h["CTYPE2"] = "DEC--TAN"
            h["CRPIX2"] = ppv.shape[1]/2.
            h["CRVAL2"] = (u.Quantity(dcen, u.deg).value, 'DEGREES')
            h["CDELT2"] = (u.Quantity(pix_size, u.arcsec).to(u.deg).value, 'DEGREES')

            h["CTYPE3"] = "VELO-LSR"
            h["CRPIX3"] = ppv.shape[2]/2.
            h["CRVAL3"] = (u.Quantity(v0, u.km/u.s).value, 'KM/S')
            h["CDELT3"] = (u.Quantity(vstep, u.km/u.s).value, 'KM/S')
            h["CUNIT3"] = "km/s" 

            h['THICK'] = (u.Quantity(thickness, u.pc).value, 'Cloud Thickness (pc)')
            h['DIST'] = (u.Quantity(dist, u.pc).value, 'Distance to cloud (pc)')
            h['V_FWHM'] = (u.Quantity(fwhm, u.km/u.s).value, 'Cloud velocity spread (km/s)')
            h['BETA'] = (beta, 'Cloud velocity power spectrum index')
            h['VEXP'] = (u.Quantity(vexp, u.km/u.s).value, 'Expansion vel (km/s - cap to midplane)')
            h['R'] = (u.Quantity(r, u.pc).value, 'Bubble size (pc)')
            h['DR'] = (u.Quantity(dr, u.pc).value, 'Bubble thickness (pc)')
            h['ZOFF'] = (u.Quantity(depth_offset, u.pc).value, 'Bubble depth offset (pc)')
            h['VOFF'] = (u.Quantity(vel_offset, u.km/u.s).value, 'Bubble-Cloud vel (km/s)')
           
            ppv = fits.PrimaryHDU(ppv, h)


        if write_fits:
            ppv.writeto(outfile, overwrite=True)

        return ppv, den, vel



def congrid(a, newdims, method='nearest'):
    """
    Regrid array `a` to new dimensions using nearest neigbor.
    
    Parameters
    ----------
    a : TYPE
        Description
    newdims : TYPE
        Description
    method : str, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    ndims = len(a.shape)
    index_list = []
    indices = np.indices(newdims)
    for i in range(ndims):
        base = indices[i]
        index_list.append(
            a.shape[i] / newdims[i] * base)
    x, y, z = np.array(index_list).round().astype(int)
    newa = a[x, y, z]
    return newa
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
    dv = vcen[1] - vcen[0]
    nchan = len(vcen)
    result = np.empty((den.shape[0], den.shape[1], nchan))

    vstart = vcen[0] - dv/2.
    voxel_channel = np.floor((vel - vstart) / dv).astype(int) #Which channel is each ppp voxel (3d pixel) in?
    voxel_valid = (voxel_channel >= 0) & (voxel_channel < nchan) #True for voxels within velocity range.
    voxel_channel_2d = voxel_channel.reshape(-1, voxel_channel.shape[-1]) # [Nx times Ny by Nz] array

    den *= voxel_valid #Remove voxels that fall outside of velocity range.

    #vbins = np.append(vcen - dv /2, vcen[-1] + dv/2) #Left edges + rightmost edge velocity channels.
    #The channel numbers in each xy pixel are scaled by a unique ID.
    voxel_channel_2d_scaled = nchan * np.arange(voxel_channel_2d.shape[0])[:,None] + voxel_channel_2d
    limit = nchan * voxel_channel_2d.shape[0]
    print(limit, voxel_channel_2d_scaled.shape, den.shape)
    ppv = np.bincount(voxel_channel_2d_scaled.ravel(), weights=den.ravel(), minlength=limit+1)[:-1]
    ppv.shape = den.shape[:-1] + (nchan,) #Reshape into a PPV cube.

    return ppv 

# def histogram_lastaxis(data, bins, weights=None):
#     """
#     Calculate histograms over the last axis of an ndarray.
#     For example, if 
#     Bins is a sequence of the left edges of each bin plus the
#     right edge of the final bin.
#     """

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
