import pyregion
import astropy.units as u


def main():
	"""
	Example call of calc_physics.
	"""
	nro_12co = "../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS_regrid0.11kms_reproj.fits"
	nro_13co = "../nro_maps/13CO_BEARS-FOREST_20170913_7.5grid_Spheroidal_Tmb_0.11kms_xy_YS.fits" 
	cube_12co = SpectralCube.read(nro_12co)
	cube_13co = SpectralCube.read(nro_13co)

	ra = 83.9292 * u.degree
	dec = -5.4631 * u.degree
	r = 0.35 * u.pc
	dr = 0.2 * u.pc                                                                            
	vexp = 3 * u.km/u.s
	v0 = 11 * u.km/u.s

	mass, momentum, energy = calc_physics(
    	ra=ra, dec=dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    	cube_12co=cube_12co, cube_13co=cube_13co, shell=True, plot=False, shell_snr_cutoff=5.)



def calc_physics(ra=None, dec=None, r=0.17*u.pc, dr=0.05*u.pc,
 vexp=4*u.km/u.s, v0=14*u.km/u.s, cube_12co=None, cube_13co=None,
 dist=414*u.pc, snr_cutoff=5., shell_snr_cutoff=3., shell=True,
 plot=False, linewidth_mode='fwhm',
 average_Tex=True):
    """
    shell_snr_cutoff is the sigma cutoff for extracting the shell
    voxels.
    """
    #print(cube_12co.header['CDELT3'])
    from spectral_cube.lower_dimensional_structures import Projection


    pix_size = (cube_12co.header['CDELT2']*u.deg).to(u.arcsec)
    vstep = (cube_12co.header['CDELT3']*(u.m/u.s)).to(u.km/u.s)

    if shell:
        model_pars = {
            'dist':dist, # pc
            'pix_size':pix_size, # arcsec
            'vstep':vstep, # km/s
            'acen':ra.to(u.deg), # deg
            'dcen':dec.to(u.deg), # deg
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

        #Extract 12co shell voxels using model.
        model_cube = SpectralCube.read(shell_model.ppv_model(**model_pars))
        #shell_masked, shell_mask = extract_shell(
        #    cube_file=cube_12co, model_pars=model_pars, return_mask=True)
        #Extract subcubes with same ra/dec range as shell voxel cube, but 
        #full velocity range.
        subcube_shell_12co = cube_12co.subcube(
            model_cube.longitude_extrema[1],
            model_cube.longitude_extrema[0],
            model_cube.latitude_extrema[0],
            model_cube.latitude_extrema[1])
        subcube_shell_13co = cube_13co.subcube(
            model_cube.longitude_extrema[1],
            model_cube.longitude_extrema[0],
            model_cube.latitude_extrema[0],
            model_cube.latitude_extrema[1])
        #print(subcube_shell_12co, subcube_shell_13co, model_cube)
        if plot:
            plt.figure()
            plt.subplot(131)
            plt.imshow(subcube_shell_12co.moment0().data)
            plt.subplot(132)
            plt.imshow(subcube_shell_13co.moment0().data)
            plt.subplot(133)
            plt.imshow(model_cube.moment0().data)
            plt.show()
    else:
        subcube_shell_12co = cube_12co
        subcube_shell_13co = cube_13co

    rms_12co = rms_map(cube=subcube_shell_12co)
    rms_13co = rms_map(cube=subcube_shell_13co)

    ###
    ### Use 13co if 3sigma detected, otherwise use corrected 12co if 3sigma detected.
    mask_use_13co = (subcube_shell_13co >= shell_snr_cutoff * rms_13co)
    mask_use_12co = (~mask_use_13co) & (subcube_shell_12co >= shell_snr_cutoff * rms_12co)   

    ### Excitation Temperature
    Tex = cube_Tex(subcube_shell_12co, average_first=True, plot=False,
        average=average_Tex)
    print("Tex is {}".format(Tex))
    ### Correct 12co for opacity.
    subcube_shell_12co_correct = opacity_correct(
        subcube_shell_12co, cube_thin=subcube_shell_13co,
        snr_cutoff=snr_cutoff, plot_ratio='ratio.png')
    #print(subcube_shell_12co_correct)

    subcube_shell_12co_correct = subcube_shell_12co_correct.with_mask(mask_use_12co)
    subcube_shell_13co = subcube_shell_13co.with_mask(mask_use_13co)
    
    # if plot:
    #     plt.figure()
    #     plt.subplot(121)
    #     plt.imshow(subcube_shell_12co_correct.moment0().data)
    #     plt.subplot(122)
    #     plt.imshow(subcube_shell_13co.moment0().data)
    #     plt.show()


    if shell:
        ### Extract shell voxels from opacity-corrected 12co
        shell_12co_correct = extract_shell(
            subcube_shell_12co_correct, keep_latlon=True, model_cube=model_cube)
        shell_13co = extract_shell(subcube_shell_13co,
            keep_latlon=True, model_cube=model_cube)
    else:
        shell_12co_correct = subcube_shell_12co_correct
        shell_13co = subcube_shell_13co

    # if plot:
    #     plt.figure()
    #     plt.subplot(121)
    #     plt.imshow(shell_12co_correct.moment0().data)
    #     plt.subplot(122)
    #     plt.imshow(shell_13co.moment0().data)
    #     plt.show()
    ###
    ### Use 13co if 3sigma detected, otherwise use corrected 12co if 3sigma detected.
    # mask_use_13co = (shell_13co >= shell_snr_cutoff * rms_13co)
    # mask_use_12co = (~mask_use_13co) & (subcube_shell_12co >= shell_snr_cutoff * rms_12co)


    ### Calculate column density of H2 from 13co where 13co is 3sig,
    ### otherwise use opacity-corrected 12co.

    
    shell_nH2_12co = column_density_H2(shell_12co_correct, Tex=Tex, molecule="12co")
    shell_nH2_13co = column_density_H2(shell_13co, Tex=Tex, molecule="13co")

    #print(shell_nH2_12co.shape, shell_nH2_13co.shape)
    print(np.nansum([shell_nH2_12co, shell_nH2_13co], axis=0), shell_nH2_12co.unit)
    shell_nH2 = Projection(np.nansum([shell_nH2_12co, shell_nH2_13co], axis=0),
     header=shell_nH2_12co.header, wcs=shell_nH2_12co.wcs, unit=shell_nH2_12co.unit,
     dtype=shell_nH2_12co.dtype, meta=shell_nH2_12co.meta, mask=shell_nH2_12co.mask)

    #print(shell_nH2.shape)

    # if plot:
    #     plt.figure()
    #     plt.subplot(131)
    #     plt.title("H2 from 12co")
    #     plt.imshow(shell_nH2_12co.data)
    #     plt.subplot(132)
    #     plt.title("H2 from 13co")
    #     plt.imshow(shell_nH2_13co.data)
    #     plt.subplot(133)
    #     plt.title("Total H2")
    #     plt.imshow(shell_nH2)


    ### Calculate Mass, Momentum, and Energy of Shell!
    if shell:
        shell_mass = mass(shell_nH2, distance=414*u.pc, molecule='H2',
         mass_unit=u.Msun)
        shell_momentum = momentum(shell_mass, vexp)
        shell_energy = energy(shell_mass, vexp)
        shell_luminosity = (shell_energy / (r / vexp)).to(u.erg/u.s)
    else:

        mass_map = mass(shell_nH2, distance=414*u.pc, molecule='H2',
            mass_unit=u.Msun, return_map=True)
        if linewidth_mode == "fwhm":
            vel_map = subcube_shell_13co.linewidth_fwhm()
        elif linewidth_mode == "sigma":
            vel_map = subcube_shell_13co.linewidth_sigma()
        elif linewidth_mode == "sigma3D":
            vel_map = 3.**0.5 * subcube_shell_13co.linewidth_sigma()

        # plt.figure()
        # plt.imshow(vel_map)
        # plt.show()
        #vel_average = np.nanmean(vel_map.value)
        shell_mass = u.Quantity(np.nansum(mass_map))
        shell_momentum = u.Quantity(np.nansum(mass_map * vel_map))
        #shell_momentum = vel
        shell_energy = 0.5 * u.Quantity(np.nansum(mass_map * vel_map * vel_map))
        #shell_energy = 0.5 * shell_mass
        #print(u.Quantity(0.5*shell_mass*u.Quantity(np.nanmean(vel_map))**2.).to(u.erg))
        if plot:
            from aplpy import FITSFigure
            from astropy.io import fits
            hdu = shell_nH2.hdu
            hdu.data = np.log10(shell_nH2.data)
            fig = FITSFigure(hdu.data)
            fig.show_colorscale(vmin=21.5, vmax=23.65, interpolation='none')
            fig.add_colorbar()
            # plt.imshow(np.log10(shell_nH2.value), interpolation='none',
            #     vmin=21.5, vmax=23.648)
            # plt.colorbar()
            # plt.title("log(n(H2) [cm^-2])")

            #plt.show()
            plt.savefig("../subregions/berneregion_lognH2.png")

            hdu = fits.open("../berne_Nh_2.fits")[0]
            hdu.data = np.log10(hdu.data)
            fig = FITSFigure(hdu.data)
            fig.show_colorscale(vmin=21.5, vmax=23.65, interpolation='none')
            fig.add_colorbar()

            #plt.show()
            plt.savefig("../subregions/berne_lognH2.png")

            # plt.figure()
            # plt.imshow(vel_map.to(u.km/u.s).value, interpolation='none')
            # #energy_map = (0.5*mass_map*vel_map*vel_map).to(u.erg).value
            # vels = vel_map.to(u.km/u.s).value.flatten()
            # plt.figure()
            # plt.hist(vels[vels>0], normed=True, log=True, bins=40)
            # plt.xlabel("Velocity FWHM (km/s)")
            # plt.ylabel("PDF")
            # plt.title('Velocity FWHM PDF')
            # plt.show()
        # plt.figure()
        # plt.imshow(vel_map.data)
        # plt.colorbar()
        # plt.show()

    return shell_mass, shell_momentum, shell_energy
    
   
#Equation 1 in Arce+ 2011, from Rohlfs and Wilson 1996
def Tex(Tpeak, thick=True):
    """
    Find the excitation temperature given peak temperature
    of an optically thick line. 
    """
    from astropy.units.core import UnitsError
    if thick:
        try:
            Tex = 5.53 / np.log(1 + (5.53)/(Tpeak+0.82))
        except UnitsError:
            Tex = 5.53*u.K / np.log(1 + (5.53*u.K)/(Tpeak+0.82*u.K))
    else:
        raise("Line must be optically thick.")
    return Tex

def cube_Tex(cube, thick=True, snr_cutoff=0,
 empty_velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s,
 average=True, average_first=False, plot=None):
    
    if snr_cutoff:
        rms = rms_map(cube, empty_velocity_range)
        cube = snr_mask(cube, rms, snr_cutoff)
    if average:
        if average_first:
            average_spec = average_spectrum(cube, return_std=False)
            Tpeak = np.max(average_spec)
        else:
            Tpeak = np.average(cube.max(axis=0))
    else:
        Tpeak = cube.max(axis=0)

    if plot:
        plt.figure()
        if average_first:
            vel = cube.spectral_axis.to(u.km/u.s).value
            plt.plot(vel, average_spec, label='Average Spectrum')
            plt.plot([vel[0], vel[-1]], [Tpeak.value, Tpeak.value], '--', label='Tpeak')
            plt.xlabel("velocity [km/s]")
            plt.ylabel("T [K]")
        else:
            plt.imshow(cube.max(axis=0).data, interpolation='none')
            plt.title("Tpeak Map: Average is {}".format(Tpeak))
            plt.colorbar(label="K")
        plt.savefig(plot)
    #print(Tpeak)
    return Tex(Tpeak, thick=thick)




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

    ratio, sigma_ratio = cube_ratio(
        cubes=[cube_thick_masked, cube_thin_masked],
        rms_maps=[rms_thick, rms_thin],
        return_uncertainty=True)

    #print(ratio.size, ratio.flattened().size)
    weights = 1. / (sigma_ratio.filled_data[:,:,:]**2.)
    #print(weights)
    average_ratio, std_ratio = average_spectrum(
        cube=ratio, weights=weights)
    #print(average_ratio, std_ratio, weights)

    if fit:
        #Fit with quadratic
        notnan = ~np.isnan(average_ratio)
        vel = ratio.spectral_axis.to(u.km/u.s).value
        #print(vel[notnan], average_ratio[notnan])
        fit_coeff, fit_cov = np.polyfit(vel[notnan], average_ratio[notnan],
            fit_order, w=(1/std_ratio[notnan]), cov=True)

        fit_func = np.poly1d(fit_coeff)
        ratio_spectrum_fit = fit_func(vel)
        cube_correct = (
        cube_thick * abundance_ratio / ratio_spectrum_fit[:,np.newaxis,np.newaxis]) #X_12,13 / (T_12/T_13)

    else:
        cube_correct = cube_thick * abundance_ratio / average_ratio[:,np.newaxis,np.newaxis]

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


def column_density_H2(cube, Tex,
    molecule="12co", transition="1-0",
    opacity_correction=True, Qrot_order=100,
    beam_filling_factor=1., moment0=True):
    """
    Calculate the column density of a molecule from an 
    optically thin emission-line cube.

    For optically thin lines only! Must correct 
    optically thick line for opacity FIRST.
    
    Can return the column density per channel (dn/dv) or 
    the total column density with a moment0 of dn/dv.
    """
    import astropy.constants as const
    import astropy.units as u
    if molecule == "12co":
        if transition == "1-0":
            #From Zhang et al. 2016
            nu_ul = 115.271 * u.GHz
            A_ul = 7.203E-8 * (1/u.s) 
            g_u = 3 # 2*J_u + 1
            E_u_k = 5.53 * u.K
            B0_k = 2.765 * u.K
            X_factor = 1e-4 # 12CO/H2
    if molecule == "13co":
        if transition == "1-0":
            #From Zhang et al. 2016
            nu_ul = 110.201 * u.GHz
            A_ul = 6.294E-8 * (1/u.s) 
            g_u = 3 # 2*J_u + 1
            E_u_k = 5.29 * u.K
            B0_k = 2.645 * u.K
            X_factor = 1e-4 / 62. # 13CO/H2 = 1.61e-6
            #X_factor = 1.43e-6 #From Berné 2014

    factor = (
        (8*np.pi*const.k_B*nu_ul**2.)\
        / (const.h*const.c**3.*A_ul*g_u)
             )
    ### APPLY THE ABUNDANCE RATIO BETWEEN CO/H2. CHECK THIS
    factor = factor * Qrot_partial(Tex, B0_k, N=Qrot_order)\
     * np.exp(E_u_k/Tex) / beam_filling_factor / X_factor


    if moment0:
        return (factor * cube.moment0()).to(1/(u.cm**2.))
    else:
        return (factor * cube).decompose()

def Qrot_partial(Tex, B0_k=2.765*u.K, N=20):
    """
    Calculate Partial sum of partition function
    at a given excitation temperature. 
    B_0/k depends on the transition:
        12CO (1-0): B_0/k = 2.765 K
    """
    #Tex = u.Quantity(Tex, u.K)
    Qrot = np.zeros(Tex.shape)
    for J in range(N+1):
        Qrot += (2*J + 1) * np.exp(
            -B0_k*J*(J+1)/Tex)
    print(Qrot)
    return Qrot



def mass(column_density, distance=414*u.pc, molecule='H2',
    return_map=False, mass_unit=u.Msun):
    """
    
    """
    if molecule == 'H2':
        mass_per_molecule = 2.34e-24*u.gram

    pixel_angle = abs(column_density.header['CDELT2']) * u.deg
    pixel_area = (pixel_angle.to(u.radian).value * distance)**2.
    #n_pixels = nH2[~np.isnan(nH2)].size
    mass_map = (column_density * mass_per_molecule * pixel_area).to(mass_unit)
    if return_map:
        return mass_map
    else:
        return u.Quantity(mass_map.nansum())



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

def cube_ratio(cubes=[None, None], rms_maps=[None, None], return_uncertainty=True):
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
