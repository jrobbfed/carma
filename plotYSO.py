import astropy.units as u
import pyfits
from astropy.coordinates import SkyCoord
import numpy as np

def cubeslice(cubefile='orion_13co.combine.fits', ralim=['5h37m30s', '5h34m30s'], declim=['-6d43m00s', '-5d54m00s'], 
vlim=[0*u.km/u.s, 20.*u.km/u.s], ra_axis=2, dec_axis=1, v_axis=0):
    """
    returns a chunk of the datacube within the specified ra, dec, velocity
    limits. uses standard fits header keywords, cdelt, crval, and naxis
    """
    f = pyfits.open(cubefile)
    head = f[0].header
    #first dimension of f[0].data represents different polarizations, this data
    #only has 1.
    data = f[0].data[0] 
    f.close()

    rastep, decstep, vstep = head['cdelt1']*u.deg, head['cdelt2']*u.deg, head['cdelt3']*u.m/u.s

    rarefpix, decrefpix, vrefpix = head['crpix1'], head['crpix2'], head['crpix3']

    raref, decref, vref = head['crval1']*u.deg, head['crval2']*u.deg, head['crval3']*u.m/u.s

    ran, decn, vn = head['naxis1'], head['naxis2'], head['naxis3']

    ra = np.linspace(raref.value - rastep.value*(rarefpix-1),
            raref.value + rastep.value*(ran-rarefpix), num=ran)*u.deg
    dec = np.linspace(decref.value - decstep.value*(decrefpix-1),
            decref.value + decstep.value*(decn-decrefpix), num=decn)*u.deg
    v = np.linspace(vref.value - vstep.value*(vrefpix-1),
            vref.value + vstep.value*(vn-vrefpix), num=vn)*u.m/u.s
    
    #find indices that correspond to the requested ra,dec,v ranges.
    clo = skycoord(ra=ralim[0], dec=declim[0])
    chi = skycoord(ra=ralim[1], dec=declim[1])

    iira = np.where((ra < clo.ra) & (ra > chi.ra))[0]  
    iidec = np.where((dec > clo.dec) & (dec < chi.dec))[0]
    iiv = np.where((v >= vlim[0]) & (v <= vlim[1]))[0]
    print iira[0], iidec, iiv
    print type(iira)     
    return data[iiv[0]:iiv[-1],iidec[0]:iidec[-1],iira[0]:iira[-1]]
   
def spectral_cubeslice(cubefile='orion_13co.combine.fits', ralim=['5h37m30s',
'5h34m30s'], declim=['-6d43m00s', '-5d54m00s'], 
vlim=[0*u.km/u.s, 20.*u.km/u.s]):
    """
    returns a chunk of the datacube within the specified ra, dec, velocity
    limits.
    this is the 2nd version, which uses the spectral_cube package to parse the
    fits file, allowing easy handling of wcs info. returns a cube object which
    retains the wcs info and the spectral-axis units.
    """
    from spectral_cube import spectralcube
    #if cubefile is a fits file, read it in with spectralcube. if cubefile is already
    #a cube object, bypass the read.

    try:
        cube = spectralcube.read(cubefile)
    except (typeerror, valueerror):
        cube = cubefile
        pass

    try:
        subcube = cube.spectral_slab(vlim[0], vlim[1])
        ra = subcube.world[0,0,:][2] 
        dec = subcube.world[0,:,0][1] 
    except (attributeerror):
        subcube = cube
        pass

    ra = subcube.world[0,0,:][2] 
    dec = subcube.world[0,:,0][1] 
    #find indices that correspond to the requested ra,dec,v ranges.
    clo = skycoord(ra=ralim[0], dec=declim[0])
    chi = skycoord(ra=ralim[1], dec=declim[1])

    iira = np.where((ra < clo.ra) & (ra > chi.ra))[0]  
    iidec = np.where((dec > clo.dec) & (dec < chi.dec))[0]
   
    return subcube[:,iidec[0]:iidec[-1],iira[0]:iira[-1]]

def ysoslice(ysofile='spitzer_orion.fits', ralim=['5h37m30s',
'5h34m30s'], declim=['-6d43m00s', '-5d54m00s'], ysotype='all'):
    f = pyfits.open(ysofile)
    data = f[1].data
    f.close()
  
    clo = skycoord(ra=ralim[0], dec=declim[0])
    chi = skycoord(ra=ralim[1], dec=declim[1])

    ra = data['ra']*u.deg
    dec = data['dec']*u.deg

    iira = (ra < clo.ra) & (ra > chi.ra)  
    iidec = (dec > clo.dec) & (dec < chi.dec)
    # if ysotype == 'all':
        
    return data[iira & iidec] 

def rmscalc(data=none, cubefile=none, ralim=none, declim=none, vlim=[15,20], 
        vunit='km/s'):
    """
    calculate the standard deviation in a section of a cube.

    data (spectral-cube.cube or str): cube object from spectral-cube package, which has wcs and spectral-cube
            axis info. if str, data should be a file which will be input into spectral_cubeslice
    #cubefile (str): fits file to be converted to a cube object.
    ralim (list[str]): limits on the right ascension, in 'hms'.
    declim (list[str]): limits on the declination, in 'dms'
    vlim (list[flt]): range to be considered for 
    """
    from spectral_cube import spectralcube

    if vunit == 'km/s':
        vlim *= u.km/u.s

    #assume data is either a cube object or a fits cube file name.    
    try:
       subcube = data.spectral_slab(vlim[0], vlim[1])
    except (attributeerror, typeerror):

        try:
            subcube = spectral_cubeslice(cubefile=data, ralim=ralim, declim=declim, vlim=vlim)
        except (attributeerror, typeerror):
            print ("rms expects either a cube object or a fits file name.")
            raise
    else:

        #slice on ra and dec ****rewrite spectral_cubeslice instead to accept both cube and filename
        if ralim != none:
            ra = subcube.world[0,0,:][2] 
            dec = subcube.world[0,:,0][1] 
            #find indices that correspond to the requested ra,dec,v ranges.
            clo = skycoord(ra=ralim[0], dec=declim[0])
            chi = skycoord(ra=ralim[1], dec=declim[1])
            iira = np.where((ra < clo.ra) & (ra > chi.ra))[0]  
            iidec = np.where((dec > clo.dec) & (dec < chi.dec))[0]
            subcube = subcube[:,iidec[0]:iidec[-1],iira[0]:iira[-1]]

    #calculate the standard deviation of all pixels in the subcube.
    return subcube.std()

def pv(data, ra=, dec=, width=0.):
    """
    data: ndarray, SpectralCube, str, HDU 
    """

######## PLOTTING FUNCTIONS #############

def plotyso(plotfile='yso.pdf', cubefile='orion_13co.combine.fits',
        ysofile='spitzer_orion.fits', ralim=['5h37m30s', '5h34m30s'],
        declim=['-6d43m00s', '-5d54m00s'], vlim=[0*u.km/u.s, 20.*u.km/u.s], 
        ysotype='all', sumtype='sum'):

    import matplotlib.pyplot as plt

    data = cubeslice(cubefile=cubefile, ralim=ralim, declim=declim, vlim=vlim)
    yso = ysoslice(ysofile=ysofile, ralim=ralim, declim=declim,
            ysotype=ysotype)
    if sumtype=='sum':
        datasum = np.sum(data, axis=0)
    if sumtype=='peak':
        datasum = np.max(data, axis=0)
    
    clo = skycoord(ra=ralim[0], dec=declim[0])
    chi = skycoord(ra=ralim[1], dec=declim[1])
    xlo, xhi = clo.ra.value, chi.ra.value
    ylo, yhi = clo.dec.value, chi.dec.value

    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(datasum, interpolation='none', extent=[xlo,xhi,ylo,yhi],
            aspect='auto')
    ax.plot(yso['ra'], yso['dec'], 'k+', markersize=8)
    
    plt.savefig(plotfile)

def plotredblue(plotfile='redblue.pdf', cubefile='l1461n.12co.fits',
        ralim=['5h37m30s', '5h34m30s'], declim=['-6d43m00s', '-5d54m00s'], bluev=[3.5, 6.5], redv=[9.5, 14.5],
        vunit='km/s', contour=true, cnorm='max',
        cmin=0.5, cmax=0.99, cstep=0.05,
        rmsralim=['5h36m50s', '5h35m50s'],
        rmsdeclim=['-6d24m','-6d2m'], rmsvlim=[15, 20], mommode='sum', 
        smooth=false, smoothwidth=20, showyso=true,
        ysofile='spitzer_orion.fits', ysotype='all', verbose=true):
    """
    plot the integrated 12co intensity in two velocity bins, corresponding
    default to Figure 6. in Nakamura et al. 2012.
    bluev: Blueshifted velocity range in km/s
    redv: Redshifted velocity range in km/s
    contour [T/F]: Make contour plot.
    cnorm ['max']: The value in each of red/blue moment image to normalize
        the contour levels to.
    cmin: Lowest contour level, in fraction of cnorm.
    cmax:
    cstep: Step in contour levels, in fraction of cnorm.

    smooth [T/F]: Do you want to gaussian smooth the moment images
    smoothwidth: Width in pixels for gaussian smoothing.
    showyso [T/F]: Over plot YSO positions on contours.
    ysofile: Filename of YSO table.
    ysotype ['all']: Type of YSO to select from YSO table.
    """

    from spectral_cube import SpectralCube, Projection
    import aplpy
    #from wcsaxes import WCS
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt 

    pass
    if vunit == 'km/s':
        if verbose:
            print "Assuming velocity inputs are in km/s..."
        bluev *= u.km/u.s
        redv *= u.km/u.s
    cube = SpectralCube.read(cubefile)

    deltav = cube.spectral_axis[1] - cube.spectral_axis[0]
    wcs_celest = WCS(cube.header).sub(['celestial'])
    bluecube = cube.spectral_slab(bluev[0], bluev[1])
    redcube = cube.spectral_slab(redv[0], redv[1])
    #Compute integrated intensity (zeroth moment) of 12CO in each red/blue
    #velocity bin.
    if mommode=='sum':
        if verbose:
            print "Calculating 0th moments using numpy.sum..."
        bluedata = bluecube.unmasked_data[:]
        reddata = redcube.unmasked_data[:]
        redmom0 = (reddata.sum(axis=0) * deltav).value #m/s
        bluemom0 = (bluedata.sum(axis=0) * deltav).value #m/s
    elif mommode == 'spectral_cube':
        if verbose:
            print "Calculating 0th moments using spectral_cube.moment (this is slow, try mommode=='sum')..."
        redmom0 = redcube.moment(order=0)
        bluemom0 = bluecube.moment(order=0)

    else:
        raise (ValueError('mommode should be either sum or spectral_cube'))
    
    if smooth:
        if verbose:
            print "Smoothing 0th moment images with a gaussian filter "+str(smoothwidth)+" pixels in width."
        from scipy import ndimage
        redmom0 = ndimage.gaussian_filter(redmom0, smoothwidth)
        bluemom0 = ndimage.gaussian_filter(bluemom0, smoothwidth)

    redhdu = Projection(redmom0, wcs=wcs_celest).hdu
    bluehdu = Projection(bluemom0, wcs=wcs_celest).hdu

    f = plt.figure()
    #Make contour level arrays in increments of percent of the maximum value. 
    if contour:
        if cnorm == 'max':
            if verbose:
                print "Normalizing the contours to the maximum value in the moment image..."
            rnorm = np.max(redmom0)
            bnorm = np.max(bluemom0)
            print rnorm, bnorm
        redlevels = np.arange(cmin*rnorm, cmax*rnorm, rnorm*cstep)
        bluelevels = np.arange(cmin*bnorm, cmax*bnorm, bnorm*cstep)
        fig = aplpy.FITSFigure(redhdu, figure=f, subplot=(1,1,1))
        redmap = fig.show_contour(data=redhdu, levels=redlevels, colors='red') 
        bluemap = fig.show_contour(data=bluehdu, levels=bluelevels, colors='blue')

    else:
        if verbose:
            print "Plotting red and blue moment maps as separate red and blue images, non-normalized..."
        figblue = aplpy.FITSFigure(bluehdu, figure=f, subplot=(1, 2, 1)) 
        figred = aplpy.FITSFigure(redhdu, figure=f, subplot=(1, 2, 2)) 
        redmap = figred.show_colorscale(cmap='Reds')
        bluemap = figblue.show_colorscale(cmap='Blues')
    if showyso:
        if verbose:
            print "Plotting YSOs from "+str(ysofile)+ "..."
        yso = ysoslice(ysofile=ysofile, ralim=ralim, declim=declim,
            ysotype=ysotype)
        ra = yso['RA'] * u.deg
        dec = yso['DEC'] * u.deg
        if contour:
            fig.show_markers(ra, dec, c='k', marker='+')
        else:
            figblue.show_markers(ra, dec, c='k', marker='+')
            figred.show_markers(ra, dec, c='k', marker='+')

    f.tight_layout()


    return redmom0, bluemom0


    

    
    
