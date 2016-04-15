#Plots various maps of the HI distribution in M33. 
#Data adapted from GALFA-HI; Putman et al. 2009 
#
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.models import custom_model_1d
from astropy.modeling.fitting import LevMarLSQFitter

M33vcut = [-320., -65] #Liberal velocity cut for M33 in km/s
                       #following Putman09
M33xcut = [50, 250] 
M33ycut = [70, 250]
#Northern Arc/Warp
northvcut = [-324, -200]
northycut = [200, 240]
northxcut = [126, 195]


#Wright's cloud 
wrightvcut = [-400, -360]
wrightxcut = [210, 330]
wrightycut = [55, 150]

#Southern Cloud
southvcut = [-170, -125]
southxcut = [145, 177]
southycut = [89, 117]

#Mystery Galactic Tendril
tvcut = [10.6, 31]
txcut = [260, 330]
tycut = [50, 110]

#Sky background velocity cut; includes no HI line
#emission
bgvcut = [-400., -365]


cube = 'M33_galfa_mod.fits' #Data cube file
pixsize = 60. #Size of pixel in arcseconds.
D_M33 = 0.730 #Distance to M33 in Mpc


def cuberead(f):
    """
    Read fits data cube, returning data array as well as an array 
    of the velocity values.
    """
    o = fits.open(f)
    head = o[0].header
    data = o[0].data
    o.close()
    Nchan = data.shape[0]
    veldelt = head['CDELT3'] / 1000. #channel width (km/s)
    vel = np.arange(Nchan) * veldelt + head['CRVAL3'] / 1000.
    return data, vel

def cubeslice(data, vel, vcut=M33vcut, xcut=M33xcut, ycut=M33ycut, vdim=0, xdim=2,
        ydim=1):
    """
    Returns a slice of the datacube, cutting on velocity, x, and y.
    vind, xind, yind specify the dimension of data corresponding to vel, x, y.
    """
    #Convert a cut in velocity (vcut) to a cut on the velocity index (ii_vcut)
    if np.size(vcut) == 2:
        ii_vcut = np.ravel(np.where((vel > vcut[0]) & (vel < vcut[1])))
        ii_vcut = [ii_vcut[0], ii_vcut[-1]] #Pull out first/last index to include.

    if np.size(vcut) != 2:
        ii_vcut = [0, np.shape(data)[vdim]] #If vcut not defined, include all.
    if np.size(xcut) != 2:
        xcut = [0, np.shape(data)[xdim]]
    if np.size(ycut) != 2:
        ycut = [0, np.shape(data)[ydim]]

    cuts = [[],[],[]]
    cuts[vdim], cuts[xdim], cuts[ydim] = ii_vcut, xcut, ycut
    
    dataslice = data[cuts[0][0]:cuts[0][1], cuts[1][0]:cuts[1][1],
            cuts[2][0]:cuts[2][1]] 
    velslice = vel[ii_vcut[0]:ii_vcut[1]]

    return dataslice, velslice

def ave(data, vel, vcut=bgvcut, xcut=M33xcut, ycut=M33ycut):
    """
    Return the average value of data within the given cuts. 
    (defaults to background section around M33)
    """
    Tb, vel = cubeslice(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)
    return np.mean(Tb, axis=0)

def rms(data, vel, vcut=bgvcut, xcut=M33xcut, ycut=M33ycut):
    """
    Returns the std. deviation of the sky brightness temperature in the area of M33. 
    """
    #Slice up the data-cube to isolate the sky beneath M33.
    Tb, vel = cubeslice(data, vel, vcut=bgvcut, xcut=xcut,
            ycut=ycut)
    sigTb = np.sqrt(np.average((Tb - np.average(Tb, axis=0))**2, axis=0))

    return sigTb

def cutoff(data, sigdata, SNRcut=0.):
    data2 = data
    data2[data < SNRcut*sigdata] = 0.
    return data2

def column(data, vel, vcut=M33vcut, bgvcut=bgvcut, xcut=M33xcut, ycut=M33ycut,
        base_subtract=True):
    """
    Find the column density of neutral hydrogen in the region given. Integrate 
    Tb along velocity and multiply by constant.
    
    SNRcut: Use as a SNR cutoff below which velocity channels are discarded
    from integration.

    base_subtract: Find average background value and subtract from data before 
    doing SNR cutoff.
    """
   
    #object slice
    Tb, velslice = cubeslice(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)

    #Find std. dev. in background slice
    sigTb = rms(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut)

    if base_subtract:
        #Subtract off mean value in background from data.
        Tb = Tb - ave(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut) 

    dv = vel[1] - vel[0] #Channel width [km/s]
    Tb = cutoff(Tb, sigTb) #Set all low SNR points equal to 0.

    N_HI = 1.38e18 * np.sum(Tb, axis=0)*dv #Column density of HI [cm^-2]

    return N_HI

def HImass(data, vel, vcut=M33vcut, bgvcut=bgvcut, xcut=M33xcut, ycut=M33ycut, D=D_M33,
        beam=pixsize, base_subtract=True, SNRcut=False):
    """
    Compute the total mass of neutral hydrogen in the region specified.

    M_HI in beam (or pixel) = 0.343 X D^2 X beam^2 X integral(T_B(v)dv), where

    M_HI: [M_sun]
    D: Distance to emission [Mpc]
    beam: Size of beam/pixel [arcseconds] (In our case, 60 arcsec/pixel)
    T_B(v): Brightness temperature [K] at velocity v [km/s]

    The integral simplifies to integral(T_B(v)dv) = channel width * sum_vcut(T_B)
    where sum_vcut(T_B) is the sum of T_B over the specified velocity range.
    channel width: vel[1] - vel[0] [km/s]

    To compute the total mass of neutral hydrogen, simply sum M_HI/pixel over
    all pixels in the specified region.
    """
    Tb, velslice = cubeslice(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)
    dv = velslice[1] - velslice[0]


    if base_subtract:
        #Subtract off mean value in background from data.
        Tb = Tb - ave(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut) 
    if SNRcut:
        #Find std. dev. in background slice
        sigTb = rms(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut)
        Tb = cutoff(Tb, sigTb) #Set all low SNR points equal to 0.
 
    
    m_xy = 0.343 * D ** 2. * beam ** 2. * dv * np.sum(Tb, axis=0) #Mass at (x,y) 
    m_tot = np.sum(m_xy) #Solar Masses
    return m_tot

def central_velocity(data, vel, vcut=M33vcut, bgvcut=bgvcut, xcut=M33xcut, ycut=M33ycut,
        base_subtract=True, peak=False):
    """
    Find the central velocity of each pixel in the the region given. Use the
    cutoff method described in van Gorkom (1989).
    Here, the central velocity of the pixel is defined as the 1st moment of the velocity 
    distribution in the pixel. That is, the intensity-weighted mean velocity.
    """
    Tb, velslice = cubeslice(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)
    #Find std. dev. in background slice
    sigTb = rms(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut)

    if base_subtract:
        #Subtract off mean value in background from data.
        Tb = Tb - ave(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut) 
    Tb = cutoff(Tb, sigTb) #Set all low SNR points equal to 0.

    if peak:
        #Simply find the peak velocity.
        v_c = velslice[np.argmax(Tb, axis=0)]
   
    else:
        #Add 2nd and 3rd empty dimensions to velslice so it can be broadcast with 
        #Tb in order to carry out the weighted average below
        v = velslice[:, np.newaxis, np.newaxis]
        v_c = np.sum(Tb * v, axis=0) / np.sum(Tb, axis=0)

    return v_c 
    
def velocity_dispersion(data, vel, vcut=M33vcut, bgvcut=bgvcut, xcut=M33xcut, ycut=M33ycut,
        base_subtract=True, peak=False):
    """
    Find the velocity dispersion of each pixel, using the
    cutoff method described in van Gorkom (1989).
    Here, the velocity dispersion is defined as the 2nd central moment of the velocity 
    distribution in the pixel. 
    """
    Tb, velslice = cubeslice(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)
    #Find std. dev. in background slice
    sigTb = rms(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut)

    if base_subtract:
        #Subtract off mean value in background from data.
        Tb = Tb - ave(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut) 
    Tb = cutoff(Tb, sigTb) #Set all low SNR points equal to 0.

    if peak:
        #Simply find the peak velocity.
        v_c = velslice[np.argmax(Tb, axis=0)]
   
    else:
        #Add 2nd and 3rd empty dimensions to velslice so it can be broadcast with 
        #Tb in order to carry out the weighted average below
        v = velslice[:, np.newaxis, np.newaxis]
        v_c = np.sum(Tb * v, axis=0) / np.sum(Tb, axis=0)
    
    sig_v = np.sqrt(np.sum(Tb * (v - v_c)**2., axis=0)/np.sum(Tb, axis=0))

    return sig_v 

def pixspec(x, y, f=cube, vcut=M33vcut):
    """
    Return the HI spectrum of one pixel over the velocity range vcut.
    """
    data, vel = cuberead(f)
    Tb, vel = cubeslice(data, vel, vcut=vcut, xcut=[], ycut=[])
    Tb = Tb[:, y, x]
    return Tb, vel
    # plt.plot(vel, Tb, '+')
    # plt.show()

def avgspec(data, vel, vcut=M33vcut, bgvcut=bgvcut, xcut=M33xcut,
        ycut=M33ycut, vdim=0, xdim=2, ydim=1, median=False):
    """
    Return an average spectrum.
    """
    Tb, vel = cubeslice(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)
    if median:
        avgTb = np.median(Tb, axis=(xdim, ydim))
    else: 
        avgTb = np.sum(Tb, axis=(xdim, ydim)) / (Tb.shape[xdim] * Tb.shape[ydim])

    return vel, avgTb


def plotspec(f=cube, plotfile='M33spec.pdf', vcut=M33vcut, xcut=M33xcut,
        ycut=M33ycut):
    """
    Plot the average HI spectrum in the given region.
    """
    data, vel = cuberead(f)
    vel, avgTb = avgspec(data, vel, vcut=vcut, xcut=xcut, ycut=ycut)
    plt.clf()
    plt.plot(vel, avgTb, '+')
    plt.title('Average HI Spectrum of M33')
    plt.xlabel('Velocity (LSR) [km/s]')
    plt.ylabel(r'Average T$_B$')
    plt.tight_layout()
    plt.savefig(plotfile)

def plotrms(f=cube, plotfile='M33rms.pdf', vcut=bgvcut, xcut=M33xcut,
        ycut=M33ycut, vrange=[]):
    """
    Plots the sigma of the sky brightness temperature in the area of M33. 
    """
    data, vel = cuberead(f)

    #Slice up the data-cube to isolate the sky beneath M33.
    sigTb = rms(data, vel, vcut=bgvcut, xcut=xcut, ycut=ycut)
    
    plt.clf()
    plt.imshow(sigTb, interpolation='none')
    cbar = plt.colorbar()
    cbar.set_label(r'RMS($T_{B}$)')
    plt.show()
    plt.draw()
    plt.savefig(plotfile)
    # plt.save_fig(plotfile)

def plotN_HI(f=cube, plotfile='M33N_HI.pdf', vcut=M33vcut, xcut=M33xcut,
        ycut=M33ycut):
    """
    Plots the intensity (really brightness temperature) integrated over all 
    velocities of M33. This is the column density of Hydrogen modulo a constant
    (Nh [cm^-2] = 1.83 X 10^18 * dv * sum(Tb(v))

    1st attempt: Simply sum over all velocity channels within vcut.
    2nd attempt: Sum over all velocity channels within vcut, with Tb >
    3*sigma(Tb), where sigma(Tb) is the std. dev. of Tb in the background.
    """

    import matplotlib.colors as colors
    data, vel = cuberead(f)

    #Calculate column density of HI in cm^-2 
    N_HI = column(data, vel, vcut=vcut, bgvcut=bgvcut, xcut=xcut, ycut=ycut,
            base_subtract=True)

    plt.clf()
    # plt.contour(N_HI, 10)
    imgplot = plt.imshow(N_HI, interpolation='none')
    cbar = plt.colorbar()
    cbar.set_label(r'N$_{HI}$ [cm$^{-2}$]')
    plt.title('Column Density') 
    plt.show()
    plt.draw()
    plt.tight_layout()
    plt.savefig(plotfile)

def plotvc(f=cube, plotfile='M33v_c.pdf', vcut=M33vcut, xcut=M33xcut,
        ycut=M33ycut, peak=False, vrange=[-290, -70]):
    """
    Plots a map of the central velocity, or intensity-weighted average
    velocity, in the selected region.  
    """
    import matplotlib.colors as colors
    data, vel = cuberead(f)

    #Calculate the intensity-weighted mean velocity in km/s 
    v_c = central_velocity(data, vel, vcut=vcut, bgvcut=bgvcut, xcut=xcut, ycut=ycut,
            base_subtract=True, peak=peak)


    plt.clf()
    # plt.contour(N_HI, 10)
    imgplot = plt.imshow(v_c, interpolation='none', vmin=vrange[0],
            vmax=vrange[1])
    cbar = plt.colorbar()
    cbar.set_label(r'$v_c$ [km/s]')
    plt.title('Central Velocity of M33') 
    plt.show()
    plt.draw()
    plt.tight_layout()
    plt.savefig(plotfile)

def plotvsig(f=cube, plotfile='M33v_sig.pdf', vcut=M33vcut, xcut=M33xcut,
        ycut=M33ycut, peak=False):
    """
    Plots a map of the central velocity, or intensity-weighted average
    velocity, in the selected region.  
    """
    import matplotlib.colors as colors
    data, vel = cuberead(f)

    #Calculate the intensity-weighted mean velocity in km/s 
    v_sig = velocity_dispersion(data, vel, vcut=vcut, bgvcut=bgvcut, xcut=xcut, ycut=ycut,
            base_subtract=True, peak=peak)


    plt.clf()
    # plt.contour(N_HI, 10)
    imgplot = plt.imshow(v_sig, interpolation='none')
    cbar = plt.colorbar()
    cbar.set_label(r'$\sigma$ [km/s]')
    # plt.contour(v_sig, colors='black')
    plt.title('Velocity Dispersion of UFO') 
    plt.show()
    plt.draw()
    plt.tight_layout()
    plt.savefig(plotfile)


@custom_model_1d
def gauss_bkgrd(x, amp=1., mean=0., sig=1., c0=0., c1=0., c2=0.):
    """
    Define a fitting function using a sum of a gaussian (for HI line) and
    a linear background function.

    amp, mean, sig: amplitude, mean, and std dev. of gaussian.
    c0, c1: y-intercept, slope of linear background function. 
    """

    return amp * np.exp(-0.5 * ((x - mean) / sig)**2) + c0 + c1 * x + c2 * x**2
