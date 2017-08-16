from spectral_cube import SpectralCube
import astropy.units as u
import shell_model
import shells
import matplotlib.pyplot as plt
#Calculate various physical quantities from
#spectral cubes, spectra, and shell parameters.
nro_12co = "../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS.fits"
nro_13co = "../nro_maps/13CO_20170518_FOREST-BEARS_spheroidal_grid7.5_dV0.11kms_xyb_YS_regridto12CO.fits"  
def main():
    n=17
    shell_list = shells.get_shells()
    shell = shell_list[n]
    #Best-fit parameters from Shell Scores spreadsheet
    #https://docs.google.com/spreadsheets/d/1rq-UZuP2PbDR9wb6woCT65GQI74NdyIJ5aWQ3CXsDhI/
    r = 0.17 * u.pc
    dr = 0.05 * u.pc
    vexp = 4 * u.km/u.s
    v0 = 14.25 * u.km/u.s
    dist = 414*u.pc
    model_pars = {
        'dist':dist, # pc
        'pix_size':7.5*u.arcsec, # arcsec
        'vstep':0.099*u.km/u.s, # km/s
        'acen':shell.ra, # deg
        'dcen':shell.dec, # deg
        'thickness':0.0, # pc
        'fwhm':0.0, # km/s
        'beta':0.0, # spectral index
        'R':r, # pc
        'dr':dr, # pc
        'vexp':vexp, # km/s
        'depth_offset':0.0, # pc
        'vel_offset':0.0, # km/s
        'v0':v0, # km/s
        'ignore_cloud':1, #Ignore cloud.
        'method':'sample',
        'write_fits':False,
        'samples_per_voxel':27}
    shell_masked = extract_shell(cube_file=nro_12co, model_pars=model_pars)

    subcube_shell = SpectralCube.read(nro_12co).subcube(
        shell_masked.longitude_extrema[1],
        shell_masked.longitude_extrema[0],
        shell_masked.latitude_extrema[0],
        shell_masked.latitude_extrema[1])
    rms_shell = rms_map(cube=subcube_shell,
     velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s)

    plt.figure()
    plt.imshow(cube_shell_masked[40].data, interpolation='none')
    plt.savefig("test_mask.png")
    plt.figure()
    plt.imshow(rms_shell.value, interpolation='none')
    plt.colorbar(label='K')
    plt.title('12CO RMS in 40 emission-free channels around T Ori Shell')
    plt.savefig("test_rms.png")
    subcube_shell_13co = SpectralCube.read(nro_13co).subcube(
        shell_masked.longitude_extrema[1],
        shell_masked.longitude_extrema[0],
        shell_masked.latitude_extrema[0],
        shell_masked.latitude_extrema[1])

    subcube_shell_12co_correct = opacity_correct(
        subcube_shell, cube_thin=subcube_shell_13co, plot_ratio="average_ratio.png")

#Equation 1 in Arce+ 2011, from Rohlfs and Wilson 1996
def Tex(Tpeak, thick=True):
    """
    Find the excitation temperature given peak temperature
    of an optically thick line. 
    """
    if thick:
        Tex = 5.53 / np.log(1 + (5.53)/(Tpeak+0.82))
    else:
        raise("Line must be optically thick.")
    return Tex


def extract_shell(cube_file=nro_12co, model_pars=None, 
    mask_minimum=0.00001):
    """
    Return a masked cube from an observed spectral cube file,
    where the mask is True wherever a model shell cube
    with parameters given by `model_pars` dictionary is > mask_minimum.

    """
    model_cube = SpectralCube.read(shell_model.ppv_model(**model_pars))
    obs_cube = SpectralCube.read(cube_file).subcube(
                                model_cube.longitude_extrema[1],
                                model_cube.longitude_extrema[0],
                                model_cube.latitude_extrema[0],
                                model_cube.latitude_extrema[1],
                                model_cube.spectral_extrema[0],
                                model_cube.spectral_extrema[1])
    #Reset the cube wcs to the values corresponding to the subcube.
    obs_cube = SpectralCube(obs_cube.hdu.data, wcs=model_cube.wcs) * u.K
    shell_mask = model_cube > mask_minimum*u.dimensionless_unscaled
    obs_cube_masked = obs_cube.with_mask(shell_mask)
    #obs_array_masked = obs_cube_masked.filled_data[:,:,:]

    return obs_cube_masked

def rms_map(cube=None, velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s):
    """
    Returns 2D array of the standard deviation of a spectral cube,
    calculated only in the emission-free channels.
    """

    channel_range = [[cube.closest_spectral_channel(vpair[0]),
                  cube.closest_spectral_channel(vpair[1])]
                 for vpair in velocity_range]

    emissionless_channels = np.concatenate(
        [np.arange(c[0], c[1]+1) for c in channel_range])
    emissionless_cube = cube.unmasked_data[emissionless_channels,:,:]
    rms_map = np.nanstd(emissionless_cube, axis=0)
    return rms_map

def regrid(cube=None, new_axis=None, smooth=True):
    """
    Regrids a SpectralCube `cube` to a new spectral axis.
    To preserve Nyquist sampling, the new spectral axis must be coarser
    than the old one, and smooth must be set to True. new_axis should have 
    units of velocity.
    See http://spectral-cube.readthedocs.io/en/latest/smoothing.html
    """
    from astropy.convolution import Gaussian1DKernel
    if smooth:
        fwhm_factor = np.sqrt(8*np.log(2)) #(FWHM/std. dev.) of gaussian
        current_resolution = cube.spectral_axis[1] - cube.spectral_axis[0]
        target_resolution = new_axis[1] - new_axis[0]
        smooth_factor = target_resolution / current_resolution

        cube = cube.spectral_smooth(Gaussian1DKernel(smooth_factor/fwhm_factor))

    cube = cube.spectral_interpolate(new_axis, 
        suppress_smooth_warning=smooth)
    return cube

def mask_snr(cube=None, rms_map=None, snr=5., return_mask=False):
    """
    Returns the spectral cube with low significance voxels
    masked out, using a map of rms calculated in emission-free
    channels.
    """
    if return_mask:
        return (cube > snr * rms_map)
    else:
        return cube.with_mask(cube > snr * rms_map)

def ratio(cubes=[None, None], rms_maps=[None, None], return_uncertainty=True):
    """
    Returns spectral cube with ratio (uncertainties on the ratio) between
    two cubes. Uncertainity calculated with error propagation assuming
    the error on each cube is given by the rms in emission-free channels.
    """
    cube_ratio = cubes[0] / cubes[1]
    if return_uncertainty:
        cube_ratio_uncertainty = cube_ratio * np.sqrt(
            (rms_maps[0] / cubes[0]) ** 2. +
            (rms_maps[1] / cubes[1]) ** 2.)
        return cube_ratio, cube_ratio_uncertainty
    else:
        return cube_ratio

def average_spectrum(cube=None, weights=None, axis=(1,2), return_std=True):
    """
    Calculate the (weighted) average spectrum in a spectral cube
    Optionally calculate and return the (weighted) standard deviation
    spectrum. `weights` should be a numpy array with the same shape as `cube`
    `axis` denotes the axis/axes to average over. For a standard SpectralCube,
    the two spatial axes are (1,2).
    Returns
    -------
    average_spectrum : 1D array_like
    std_spectrum: 1D array_like, optional
    """
    average_spectrum = np.average(cube.filled_data[:,:,:],
     weights=weights, axis=axis)
    if return_std:
        resids_squared = (cube - average_spectrum[:,np.newaxis,np.newaxis])**2.
        std_spectrum = np.sqrt(
            np.average(resids_squared, weights=weights, axis=axis))
        return average_spectrum, std_spectrum
    else:
        return average_spectrum

def opacity_correct(cube_thick, cube_thin=None, abundance_ratio=62.,
    snr_cutoff=5., empty_velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s,
    regrid_cube=False, plot_ratio=None,
    fit=True, fit_order=2, **kwargs):
    """
    Correct an optically thick emission line cube using an (assumed) optically
    thin emission line cube. The abundance ratio betweeen the cube and cube_thin
    isotopologues is given as `abundance_ratio`.
    
    `regrid_cube`: Optionally regrid `cube_thin` to the velocity grid of `cube`, preserving
    Nyquist sampling.
    
    Uses method detailed in Zhang et al 2016 (c.f. Dunham+14, Arce+2001)
    """
    if regrid_cube:
        cube_thin = regrid(cube_thin, cube_thick.spectral_axis)

    rms_thick = rms_map(cube_thick, empty_velocity_range)
    rms_thin = rms_map(cube_thin, empty_velocity_range)

    cube_thick_masked = mask_snr(cube_thick, rms_thick,
     snr=snr_cutoff)
    cube_thin_masked = mask_snr(cube_thin, rms_thin,
     snr=snr_cutoff)

    ratio, sigma_ratio = ratio(
        cubes=[cube_thick_masked, cube_thin_masked],
        rms_maps=[rms_thick, rms_thin],
        return_uncertainty=True)


    weights = 1. / (sigma_ratio.filled_data[:,:,:]**2.)

    average_ratio, std_ratio = average_spectrum(
        cube=ratio, weights=weights)

    if fit:
        #Fit with quadratic
        notnan = ~np.isnan(average_ratio)
        vel = ratio.spectral_axis.to(u.km/u.s).value
        fit_coeff, fit_cov = np.polyfit(vel[notnan], average_ratio[notnan],
            fit_order, w=(1/std_ratio[notnan]), cov=True)

        fit_func = np.poly1d(fit_coeff)
        ratio_spectrum_fit = fit_func(vel)
        cube_correct = cube_thick * abundance_ratio / ratio_spectrum_fit #X_12,13 / (T_12/T_13)

    else:
        cube_correct = cube_thick * abundance_ratio / average_ratio

    if plot_ratio: 
        plt.plot(vel, average_ratio, 'o')
        if fit:
            xfit = np.linspace(vel[notnan][0], vel[notnan][-1], 1000)
            yfit = fit_func(xfit)
            plt.plot(xfit, yfit)
        plt.fill_between(vel,
            average_ratio+std_ratio, average_ratio-std_ratio,
                alpha=0.3)
        plt.ylabel(r"T$_{12}$ / T$_{13}$", size=16)
        plt.xlabel("Velocity [km / s]", size=16)
        #plt.ylim(1.3,3.9)
        plt.title("Weighted Average and Std. Dev. of 5sigma 12CO/13CO Ratio")
        plt.savefig(plot_ratio)

    return cube_correct

def column_density(opacity_correction=True):
    """
    For optically thick line, 
    """

def momentum(mass, velocity):
    return mass * velocity
def energy(mass, velocity):
    return 0.5 * mass * velocity ** 2.


if __name__ == '__main__':
    main()