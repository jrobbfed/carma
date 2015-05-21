import astropy.units as u
import pyfits
from astropy.coordinates import SkyCoord
import numpy as np

def cubeslice(cubefile='orion_13co.combine.fits', ralim=['5h37m30s', '5h34m30s'], declim=['-6d43m00s', '-5d54m00s'], 
vlim=[0*u.km/u.s, 20.*u.km/u.s], ra_axis=2, dec_axis=1, v_axis=0):
    """
    Returns a chunk of the datacube within the specified RA, Dec, velocity
    limits. Uses standard FITS header keywords, CDELT, CRVAL, and NAXIS
    """
    f = pyfits.open(cubefile)
    head = f[0].header
    #First dimension of f[0].data represents different polarizations, this data
    #only has 1.
    data = f[0].data[0] 
    f.close()

    rastep, decstep, vstep = head['CDELT1']*u.deg, head['CDELT2']*u.deg, head['CDELT3']*u.m/u.s

    rarefpix, decrefpix, vrefpix = head['CRPIX1'], head['CRPIX2'], head['CRPIX3']

    raref, decref, vref = head['CRVAL1']*u.deg, head['CRVAL2']*u.deg, head['CRVAL3']*u.m/u.s

    raN, decN, vN = head['NAXIS1'], head['NAXIS2'], head['NAXIS3']

    ra = np.linspace(raref.value - rastep.value*(rarefpix-1),
            raref.value + rastep.value*(raN-rarefpix), num=raN)*u.deg
    dec = np.linspace(decref.value - decstep.value*(decrefpix-1),
            decref.value + decstep.value*(decN-decrefpix), num=decN)*u.deg
    v = np.linspace(vref.value - vstep.value*(vrefpix-1),
            vref.value + vstep.value*(vN-vrefpix), num=vN)*u.m/u.s
    
    #Find indices that correspond to the requested ra,dec,v ranges.
    clo = SkyCoord(ra=ralim[0], dec=declim[0])
    chi = SkyCoord(ra=ralim[1], dec=declim[1])

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
    Returns a chunk of the datacube within the specified RA, Dec, velocity
    limits.
    This is the 2nd version, which uses the spectral_cube package to parse the
    fits file, allowing easy handling of WCS info. Returns a Cube object which
    retains the WCS info and the spectral-axis units.
    """
    from spectral_cube import SpectralCube
    cube = SpectralCube.read(cubefile)
    subcube = cube.spectral_slab(vlim[0], vlim[1])

    ra = subcube.world[0,0,:][2] 
    dec = subcube.world[0,:,0][1] 
    #Find indices that correspond to the requested ra,dec,v ranges.
    clo = SkyCoord(ra=ralim[0], dec=declim[0])
    chi = SkyCoord(ra=ralim[1], dec=declim[1])

    iira = np.where((ra < clo.ra) & (ra > chi.ra))[0]  
    iidec = np.where((dec > clo.dec) & (dec < chi.dec))[0]
   
    return subcube[:,iidec[0]:iidec[-1],iira[0]:iira[-1]]

def ysoslice(ysofile='spitzer_orion.fits', ralim=['5h37m30s',
'5h34m30s'], declim=['-6d43m00s', '-5d54m00s'], ysotype='all'):
    f = pyfits.open(ysofile)
    data = f[1].data
    f.close()
  
    clo = SkyCoord(ra=ralim[0], dec=declim[0])
    chi = SkyCoord(ra=ralim[1], dec=declim[1])

    ra = data['RA']*u.deg
    dec = data['DEC']*u.deg

    iira = (ra < clo.ra) & (ra > chi.ra)  
    iidec = (dec > clo.dec) & (dec < chi.dec)
    # if ysotype == 'all':
        
    return data[iira & iidec] 

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
    
    clo = SkyCoord(ra=ralim[0], dec=declim[0])
    chi = SkyCoord(ra=ralim[1], dec=declim[1])
    xlo, xhi = clo.ra.value, chi.ra.value
    ylo, yhi = clo.dec.value, chi.dec.value

    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(datasum, interpolation='none', extent=[xlo,xhi,ylo,yhi],
            aspect='auto')
    ax.plot(yso['RA'], yso['DEC'], 'k+', markersize=8)
    
    plt.savefig(plotfile)
    figure = plt.figure()
    df

def plotredblue(plotfile='redblue.pdf', cubefile='l1461n.12co.fits',
        ralim=None, declim=None, bluev=[3.5, 6.5], redv=[9.5, 14.5],
        vunit='km/s', contour=True):
    #Plot the integrated 12CO intensity in two velocity bins, corresponding
    #default to Figure 6. in Nakamura et al. 2012.
    #bluev: Blueshifted velocity range in km/s
    #redv: Redshifted velocity range in km/s
    from spectral_cube import SpectralCube
    import aplpy
    from astropy.wcs import WCS

    if vunit == 'km/s':
        bluev *= u.km/u.s
        redv *= u.km/u.s
    cube = SpectralCube.read(cubefile)
    deltav = cube.spectral_axis[1] - cube.spectral_axis[0]
    bluecube = cube.spectral_slab(bluev[0], bluev[1])
    redcube = cube.spectral_slab(redv[0], redv[1])
    # wcs = WCS(redcube.header) 
    bluedata = bluecube.unmasked_data[:]
    reddata = redcube.unmasked_data[:]

    #Compute integrated intensity (zeroth moment) of 12CO in each red/blue
    #velocity bin.
    redmom0 = reddata.sum(axis=0) * deltav
    bluemom0 = bluedata.sum(axis=0) * deltav
    # redmom0 = redcube.moment(order=0)
    # bluemom0 = bluecube.moment(order=0)
    print redmom0, bluemom0

    

    
    
