"""Summary
"""
#from shells import pv_slice
import numpy as np
import astropy.units as u
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
from astropy.io import fits
from spectral_cube import SpectralCube
from aplpy import FITSFigure
import shells
import matplotlib.pyplot as plt

def main():
    """Summary
    """
    
    import numpy as np
    import os
    from astropy.io import fits
    import astropy.units as u
    import warnings
    import numpy as np
    from scipy.ndimage.filters import gaussian_filter
    import shell_model

    cubefile12co = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
    regionfile = "../shell_candidates/AllShells.reg"

    # For Shell 18
    n=17
    shell_list = shells.get_shells()
    shell = shell_list[n]

    outfile = '../turbulent_model/shell18_nocloud.fits'
    model_pars = {
        'outfile':"'{}'".format(outfile),
        'dist':414*u.pc, # pc
        'pix_size':7.5*u.arcsec, # arcsec
        'vstep':0.099*u.km/u.s, # km/s
        'acen':shell.ra.value, # deg
        'dcen':shell.dec.value, # deg
        'thickness':0.0, # pc
        'fwhm':0.0, # km/s
        'beta':0.0, # spectral index
        'R':0.22*u.pc, # pc
        'dr':0.2*u.pc, # pc
        'vexp':2.2*u.km/u.s, # km/s
        'depth_offset':0.0, # pc
        'vel_offset':0.0, # km/s
        'v0':13.6*u.km/u.s, # km/s
        'ignore_cloud':1, #Ignore cloud.
        'method':'sample',
        'write_fits':False,
        'samples_per_voxel':27
        }

    ppv = shell_model.ppv_model(dist=model_pars['dist']*u.pc, pix_size=model_pars['pix_size']*u.arcsec,\
                                             vstep=model_pars['vstep']*u.km/u.s, acen=shell.ra, dcen=shell.dec,\
                                             R=model_pars['R']*u.pc, dr=model_pars['dr']*u.pc,\
                                             vexp=model_pars['vexp']*u.km/u.s, v0=model_pars['v0']*u.km/u.s,\
                                             interpolate_ppv=True)


    pv = pv_average(
        cube=SpectralCube.read(ppv_interp),
        ra_center=shell.ra, dec_center=shell.dec,
        width=pv_width, length=pv_length, angle_step=10*u.deg)

    fig = plt.figure(figsize=(10,10))
    fig1 = FITSFigure(pv, figure=fig, subplot=(1,1,1))
    fig1.show_grayscale()
    fig1.savefig("pv_interp.png")

def compare(plot_file=None, cube_file="../nro_maps/../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits",
    regionfile="../shell_candidates/AllShells.reg", 
    pv_width_in_pixels=3., pv_length_in_radii=3., pv_angle=0*u.deg, pv_angle_step=90*u.deg,
    average_pv_obs=True, average_pv_model=False, model_pars=None,
    normalize=True, contour_levels=20, draw_radius=True):
    """Summary
    
    Parameters
    ----------
    cube_file : str, optional
        Observed cube file to read.
    regionfile : str, optional
        Description
    pv_width_pixels : float, optional
        Description
    pv_length_radii : float, optional
        Description
    pv_angle : TYPE, optional
        Description
    pv_angle_step : TYPE, optional
        Description
    average_pv_obs : bool, optional
        Description
    average_pv_model : bool, optional
        Description
    model_pars : dict, optional
        Parameters to pass to `ppv_model`. If None, 
        the parameters are the defaul values in `ppv_model`.
    
    Deleted Parameters
    ------------------
    model_pars : TYPE, optional
        Description
    """
    try:
        cube_model = SpectralCube.read(ppv_model(**model_pars))
    except:
        cube_model = SpectralCube.read(ppv_model())

    head_model = cube_model.header
    cube_obs = SpectralCube.read(cube_file).subcube(
                                cube_model.longitude_extrema[1],
                                cube_model.longitude_extrema[0],
                                cube_model.latitude_extrema[0],
                                cube_model.latitude_extrema[1],
                                cube_model.spectral_extrema[0],
                                cube_model.spectral_extrema[1])
    pv_width = pv_width_in_pixels * head_model['CDELT1'] * u.Unit(head_model['CUNIT1'])
    pv_length = pv_length_in_radii * (head_model['R'] / head_model['DIST']) * u.radian
    ra_center = head_model['CRVAL1'] * u.deg
    dec_center = head_model['CRVAL2'] * u.deg

    if average_pv_obs:
        pv_obs = shells.pv_average(cube=cube_obs,
         ra_center=ra_center, dec_center=dec_center,
         width=pv_width, length=pv_length, angle_step=pv_angle_step)
    else:
        pv_obs = shells.pv_slice(cube=cube_obs,
         ra_center=ra_center, dec_center=dec_center,
         width=pv_width, length=pv_length, angle=pv_angle)

    if average_pv_model:
        pv_model = shells.pv_average(cube=cube_model,
         ra_center=ra_center, dec_center=dec_center,
         width=pv_width, length=pv_length, angle_step=pv_angle_step)
    else:
        pv_model = shells.pv_slice(cube=cube_model,
         ra_center=ra_center, dec_center=dec_center,
         width=pv_width, length=pv_length, angle=pv_angle)

    if normalize:
        pv_obs.data /= np.nanmax(pv_obs.data)
        pv_model.data /= np.nanmax(pv_model.data)

    fig = plt.figure(figsize=(20,10))

    title = """angle = {}, length = {} , width = {},
        r = {} pc, dr = {} pc, v0 = {} km/s, vexp = {} km/s""".format(
        pv_angle.round(2), pv_length.to(u.arcmin).round(2), pv_width.to(u.arcsec).round(2),
        head_model['r'], head_model['dr'], head_model['v0'], head_model['vexp'])

    figpv = FITSFigure(pv_obs, figure=fig, subplot=((1,2,1)))
    figpv.show_grayscale(aspect='auto')
    figpv.show_contour(pv_model, levels=contour_levels)
    figpv.set_title(title)

    figmom0 = FITSFigure(cube_obs.moment0().hdu, figure=fig, subplot=((1,2,2)))
    figmom0.show_grayscale()
    #figmom0.show_contour(cube_model.moment0().hdu, levels=int(contour_levels/2))
    figmom0.set_title(title)

    if draw_radius:
        r_degrees = (head_model['R'] / head_model['DIST']) * 360. / (2 * np.pi) 
        dr_degrees = (head_model['DR'] / head_model['DIST']) * 360. / (2 * np.pi)

        figmom0.show_circles(ra_center.value, dec_center.value, r_degrees)
        figmom0.show_circles(ra_center.value, dec_center.value, r_degrees - dr_degrees / 2.,
         linestyle='--', edgecolor='red')
        figmom0.show_circles(ra_center.value, dec_center.value, r_degrees + dr_degrees / 2.,
         linestyle='--', edgecolor='red')

    if plot_file:
        fig.savefig(plot_file)

    return fig

def ppv_model(outfile=None, dist=414*u.pc, pix_size=7.5*u.arcsec,
    vstep=0.099*u.km/u.s, acen=83.72707*u.deg, dcen=-5.07792*u.deg,
    thickness=0.0*u.pc, fwhm=0.0*u.km/u.s, beta=0.0, R=0.22*u.pc,
    dr=0.2*u.pc, vexp=2.2*u.km/u.s, depth_offset=0.*u.pc, 
    vel_offset=0.*u.km/u.s, v0=13.6*u.km/u.s,
    ignore_cloud=True, write_fits=False, return_hdu=True,
    return_ppp=False, downsample=True, smooth=True, smooth_fwhm=2.
    working_grid_factor=2.,
    interpolate_ppv=False, method='sample', samples_per_voxel=27.,
    pad_pixels=5, pad_channels=5):
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
    R : TYPE, optional
        Description
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
    ignore_cloud : bool, optional
        If True, ignore cloud completely. If not, will need to 
        call to idl program.
    write_fits : bool, optional
        Description
    return_hdu : bool, optional
        Description
    return_ppp : bool, optional
        Description
    downsample : bool, optional
        Description
    smooth : bool, optional
        Description
    working_grid_factor : float, optional
        Description
    interpolate_ppv : bool, optional
        Description
    method : str, optional
        Description
    samples_per_voxel : float, optional
        Description
    pad_pixels : int, optional
        Description
    pad_channels : int, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    
    Deleted Parameters
    ------------------
    r : number or length-type astropy.Unit, optional
        Radius of shell. Assumes pc if not specified.
        Radius is measured from center to midpoint of shell rim.
    chan_pad : float, optional
        Factor to pad the ppv cube on either end in velocity channel space.
    """

    if ignore_cloud:
   # Work in a pixel scale 2x finer than end result. 
        if method == "sample":
            """
            So, we don’t need to create I(x,y,z) and v(x,y,z) grids, since I(x,y,vz) only depends on x, y, z and the constants v_exp and v_off.

    - Start off with a resolution greater than needed at the end, then can bin the resulting (x, y, vz) space to the requested resolution.

    We should not need any interpolation if working with an appropriate resolution.

    - This method will treat the expanding shell as many individual identical particles with (x,y,z,vz).

    Bin into voxels of size [dx,dy,channel_width]
            """
            
            pix_pc = u.Quantity(dist, u.pc) * u.Quantity(pix_size, u.radian).value # pc per pixel
            shell_volume = (4./3.) * np.pi * ((R+dr/2)**3. - (R-dr/2)**3.)
            n_points = int((shell_volume / pix_pc ** 3.) * samples_per_voxel)
            #Sample a Unit Sphere
            theta = np.random.uniform(0, 2*np.pi, n_points)
            z0 = np.random.uniform(-1., 1., n_points)
            x0 = np.sqrt(1 - z0**2.)*np.cos(theta)
            y0 = np.sqrt(1 - z0**2.)*np.sin(theta)

            #Sample radius from distribution that has uniform volume density
            r = np.random.uniform((R.value - dr.value/2)**3., (R.value + dr.value/2)**3., n_points) ** (1./3.)
            x, y, z = x0*r, y0*r, z0*r

            #Z-component of velocities
            vz = vexp * z / r + v0

            #Make spatial and velocity bins with midpoint in center of middle bin.
            #pad_pixels = 5
            pix_start = -R - dr/2 - pad_pixels*pix_pc
            pix_end = R + dr/2 + pad_pixels*pix_pc
            pix_bins = symmetric_bins(pix_start.value, pix_end.value, pix_pc.value)

            #pad_channels = 5
            vstart = v0 - vexp - pad_channels*vstep
            vend = v0 + vexp + pad_channels*vstep
            vz_bins = symmetric_bins(vstart.value, vend.value, vstep.value)

            ppv, edges = np.histogramdd((x, y, vz.value), bins=(pix_bins, pix_bins, vz_bins))

        else:
            scale = u.Quantity(dist, u.pc) * u.Quantity(pix_size, u.radian).value / working_grid_factor # pc per pixel
            box_size = np.floor(4 * u.Quantity(R, u.pc) / scale) # Number of pixels per side of ppp box.

            # den = np.zeros(3*[box_size])
            # x, y, z = np.indices(3*[box_size]) - np.floor(box_size / 2) #Index w.r.t center of box.
            x, y, z = np.indices(3*[box_size]) - np.floor(box_size / 2) #Index w.r.t center of box.

            rr = np.sqrt(x ** 2. + y ** 2. + z ** 2.) * scale 

            #inside = rr < r - dr/2.
            #on = (rr >= r - dr/2.) & (rr < r + dr/2.) 
            #outside = rr >= r + dr/2.

            #ratio = np.sum(inside) / np.sum(on)
            den = ((rr >= R - dr/2.) & (rr < R + dr/2.)).astype(int) #Density is 1 on rim of shell, 0 elsewhere.

            #Velocity field of bubble.
            #This is the z-component of velocity, assuming l.o.s. is in z direction. 
            #From similar triangles, z / r = v_z / v_r : v_r -> vexp, v_z -> vel, r*scale -> rr, z*scale -> z
            vel = den * vexp * z * scale / rr      
            vel[rr == 0] = 0 #Fix pole.

            vlo, vhi = -vexp - vstep * pad_channels, vexp + vstep * pad_channels
            vcen = np.linspace(vlo, vhi, (vhi - vlo) / vstep)

            #Gridding PPV cube.
            ppv = ppp2ppv(den, vel.value, vcen.value, interpolate=interpolate_ppv)
            if downsample:
                ppv = congrid(ppv, (ppv.shape[0]//working_grid_factor, ppv.shape[1]//working_grid_factor, ppv.shape[2]))

        if smooth:
            gauss = Gaussian2DKernel(stddev = smooth_fwhm / np.sqrt(8 * np.log(2))) #FWHM to std. dev. of gaussian.

            for i in range(ppv.shape[2]):
                 ppv[:,:,i] = convolve_fft(ppv[:,:,i], gauss, normalize_kernel=True)
        
        ppv = np.swapaxes(ppv, 0, 2)
 
        if return_hdu:
            h = fits.Header()
            h["CTYPE1"] = "RA---TAN"
            h["CRPIX1"] = ppv.shape[2]/2. + 0.5
            h["CRVAL1"] = (u.Quantity(acen, u.deg).value, 'DEGREES')
            h["CDELT1"] = (u.Quantity(pix_size, u.arcsec).to(u.deg).value, 'DEGREES')

            h["CTYPE2"] = "DEC--TAN"
            h["CRPIX2"] = ppv.shape[1]/2. + 0.5
            h["CRVAL2"] = (u.Quantity(dcen, u.deg).value, 'DEGREES')
            h["CDELT2"] = (u.Quantity(pix_size, u.arcsec).to(u.deg).value, 'DEGREES')

            h["CTYPE3"] = "VELO-LSR"
            h["CRPIX3"] = ppv.shape[0]/2. + 0.5
            h["CRVAL3"] = (u.Quantity(v0, u.km/u.s).value, 'KM/S')
            h["CDELT3"] = (u.Quantity(vstep, u.km/u.s).value, 'KM/S')
            h["CUNIT3"] = "km/s" 

            h['THICK'] = (u.Quantity(thickness, u.pc).value, 'Cloud Thickness (pc)')
            h['DIST'] = (u.Quantity(dist, u.pc).value, 'Distance to cloud (pc)')
            h['V_FWHM'] = (u.Quantity(fwhm, u.km/u.s).value, 'Cloud velocity spread (km/s)')
            h['BETA'] = (beta, 'Cloud velocity power spectrum index')
            h['VEXP'] = (u.Quantity(vexp, u.km/u.s).value, 'Expansion vel (km/s - cap to midplane)')
            h['R'] = (u.Quantity(R, u.pc).value, 'Bubble size (pc)')
            h['DR'] = (u.Quantity(dr, u.pc).value, 'Bubble thickness (pc)')
            h['ZOFF'] = (u.Quantity(depth_offset, u.pc).value, 'Bubble depth offset (pc)')
            h['VOFF'] = (u.Quantity(vel_offset, u.km/u.s).value, 'Bubble-Cloud velocity offset (km/s)')
            h['V0'] = (u.Quantity(v0, u.km/u.s).value, "Systemic Velocity (km/s)")
           
            ppv = fits.PrimaryHDU(ppv, h)


        if write_fits:
            ppv.writeto(outfile, overwrite=True)
        if return_ppp:
            return ppv, den, vel
        else:
            return ppv

def symmetric_bins(start, end, step):
    """
    Returns bin edges that are symmetric around the `mid` value,
    with the `mid` value falling into the midpoint of the 
    central bin. mid = (end - start) / 2
    
    Parameters
    ----------
    start : TYPE
        Description
    end : TYPE
        Description
    step : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    mid = (end + start) / 2.
    bins_upper = np.arange(mid + step/2.,
                           end + step/2.,
                           step)
    bins_lower = np.arange(mid - step/2.,
                           start - step/2.,
                           -step)[::-1]
    return np.append(bins_lower, bins_upper)

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
def ppp2ppv(den, vel, vcen, interpolate=False):
    """Convert density and velocity cubes into a
     Position-Position-Velocity cube. 
    
    Parameters
    ----------
    den : 3D array-like
        Density cube in xyz space.
    vel : 3D array-like
        Velocity cube in xyz space.
    vcen : array-like
        Array of velocity channel centers. Must be equally spaced.
    interpolate : bool, optional
        If True, implement an interpolation along the velocity axis
        cloning the implementation in the IDL program ppp2ppv.pro
    
    Returns
    -------
    3D array-like
        Position-Position-Velocity cube.
    """
    dv = vcen[1] - vcen[0]
    nchan = len(vcen)

    vstart = vcen[0] - dv/2.
    voxel_channel = np.floor((vel - vstart) / dv).astype(int) #Which channel is each ppp voxel (3d pixel) in?
    voxel_valid = (voxel_channel >= 0) & (voxel_channel < nchan) #True for voxels within velocity range.
    voxel_channel_2d = voxel_channel.reshape(-1, voxel_channel.shape[-1]) # [Nx times Ny by Nz] array

    den *= voxel_valid #Remove voxels that fall outside of velocity range.

    if interpolate:
        # x = np.arange(den.shape[0] * den.shape[1]) % den.shape[0]
        # y = np.arange(den.shape[0] * den.shape[1]) // den.shape[0]
         
        # for i in range(den.shape[2]):
        #     z = x * 0 + i
        #     if i == den.shape[2] - 1:                
        #         z_next = z
        #     else:
        #         z_next = z + 1 
        #     chan = voxel_channel[x, y, z]
        #     chan_next = voxel_channel[x, y, z_next]
        ppv = np.empty((den.shape[0], den.shape[1], nchan))
        nz = den.shape[2]
        for z in range(nz):
            print(z)
            z_next = min(z + 1, nz - 1)
            chan = voxel_channel[:, :, z]
            chan_next = voxel_channel[:, :, z_next]
            den_z = den[:, :, z]

            # Use maximum absolute difference between channels in z-adjacent voxels
            # to set how finely to resample and redistribute density.
            jump = 3 * int(np.max(abs(chan - chan_next)) + 1)

            for j in range(jump):
                weight = 1.0 * j / jump
                chan_intermediate = np.floor(chan * (1 - weight) + chan_next * weight).astype(int)
                ppv[:, :, chan_intermediate] += den_z / float(jump) 

    else:
        #vbins = np.append(vcen - dv /2, vcen[-1] + dv/2) #Left edges + rightmost edge velocity channels.
        #The channel numbers in each xy pixel are scaled by a unique ID.
        voxel_channel_2d_scaled = nchan * np.arange(voxel_channel_2d.shape[0])[:,None] + voxel_channel_2d
        limit = nchan * voxel_channel_2d.shape[0]
        
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

def pv_average(cube=None, ra_center=None, dec_center=None,
    width=22.5*u.arcsec, length=5*u.arcmin,
    angle_range=[0*u.deg, 360.*u.deg], angle_step=10*u.deg,
    mode='average'):
    """Returns a postion-velocity slice of `cube` from pvextractor,
    averaged over `angle`. If `angle_step` == None, step by `width`
    
    Parameters
    ----------
    cube : None, optional
        Description
    ra_center : None, optional
        Description
    dec_center : None, optional
        Description
    width : TYPE, optional
        Description
    length : TYPE, optional
        Description
    angle_range : TYPE, optional
        Description
    angle_step : TYPE, optional
        Description
    mode : str, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """


    import shells
    import numpy as np
    from astropy.io import fits
    
    angle_list = np.linspace(angle_range[0], angle_range[1],
                             (angle_range[1] + angle_step) / angle_step)
    pv_list = [shells.pv_slice(cube=cube, ra_center=ra_center, dec_center=dec_center,
                               angle=angle, width=width, length=length)
               for angle in angle_list]
    if mode == 'average':
        average_data = np.mean(np.array([hdu.data for hdu in pv_list]), axis=0)
        print(average_data.shape)
        average_HDU = fits.PrimaryHDU(average_data, header=pv_list[0].header)
        return average_HDU


if __name__ == "__main__":
    main()
