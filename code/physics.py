from spectral_cube import SpectralCube
import astropy.units as u
import shell_model
import shells
import matplotlib.pyplot as plt
import numpy as np
import glob
#Calculate various physical quantities from
#spectral cubes, spectra, and shell parameters.
nro_12co = "../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS_regrid0.11kms_reproj.fits"
nro_13co = "../nro_maps/13CO_BEARS-FOREST_20170913_7.5grid_Spheroidal_Tmb_0.11kms_xy_YS.fits" 
nro_13co_divided = "../nro_maps/13CO_20170518_FOREST-BEARS_spheroidal_grid7.5_dV0.11kms_xyb_YS_regridto12CO_divide1.4.fits" 
best_shells = [3,6,9,11,17,18,21,24,25,30,37,38]
def main():
    dist = 414*u.pc
    plot_physicsrange(column=1, plotname="massrange_all.png")
    plot_physicsrange(column=2, plotname="momentumrange_all.png")
    plot_physicsrange(column=3, plotname="energyrange_all.png")
    # all_low_tables = glob.glob("shell*properties*low_13co_1.4.txt")
    # all_hi_tables = glob.glob("shell*properties*hi_13co_1.4.txt")
    # robust_low_tables = glob.glob("shell[369]_properties_low_13co_1.4.txt")\
    #  + glob.glob("shell1[178]_properties_low_13co_1.4.txt")\
    #  + glob.glob("shell2[145]_properties_low_13co_1.4.txt")\
    #  + glob.glob("shell3[078]_properties_low_13co_1.4.txt")
    # robust_hi_tables = glob.glob("shell[369]_properties_hi_13co_1.4.txt")\
    #  + glob.glob("shell1[178]_properties_hi_13co_1.4.txt")\
    #  + glob.glob("shell2[145]_properties_hi_13co_1.4.txt")\
    #  + glob.glob("shell3[078]_properties_hi_13co_1.4.txt")

    # all_low_tables = glob.glob("shell*properties*low.txt")
    # all_hi_tables = glob.glob("shell*properties*hi.txt")
    # robust_low_tables = glob.glob("shell[369]_properties_low.txt")\
    #  + glob.glob("shell1[178]_properties_low.txt")\
    #  + glob.glob("shell2[145]_properties_low.txt")\
    #  + glob.glob("shell3[078]_properties_low.txt")
    # robust_hi_tables = glob.glob("shell[369]_properties_hi.txt")\
    #  + glob.glob("shell1[178]_properties_hi.txt")\
    #  + glob.glob("shell2[145]_properties_hi.txt")\
    #  + glob.glob("shell3[078]_properties_hi.txt")
    # mass_all = hist_physics(table_list=all_low_tables+all_hi_tables,
    #     column=1, table_list_shaded=robust_low_tables+robust_hi_tables,
    #     plotname='hist_mass_all.png')
    # momentum_all = hist_physics(table_list=all_low_tables+all_hi_tables,
    #     column=2, table_list_shaded=robust_low_tables+robust_hi_tables,
    #     plotname='hist_momentum_all.png')
    # energy_all = hist_physics(table_list=all_low_tables+all_hi_tables,
    #     column=3, table_list_shaded=robust_low_tables+robust_hi_tables,
    #     plotname='hist_energy_all.png')

    # mass_all_low = hist_physics(table_list=all_low_tables, column=1,
    #     plotname='hist_mass_low_all_new13co.png', table_list_shaded=robust_low_tables
    #     #xlim=[-50,300]
    #     )
    # mass_all_hi = hist_physics(table_list=all_hi_tables, column=1,
    #     plotname='hist_mass_hi_all_new13co.png', table_list_shaded=robust_hi_tables
    #     #xlim=[-50,300]
    #     )
    # momentum_all_low = hist_physics(table_list=all_low_tables, column=2,
    #     plotname='hist_momentum_low_all_new13co.png', table_list_shaded=robust_low_tables)
    # momentum_all_hi = hist_physics(table_list=all_hi_tables, column=2,
    #     plotname='hist_momentum_hi_all_new13co.png', table_list_shaded=robust_hi_tables)
    # energy_all_low = hist_physics(table_list=all_low_tables, column=3,
    #     plotname='hist_energy_low_all_new13co.png', table_list_shaded=robust_low_tables)
    # energy_all_hi = hist_physics(table_list=all_hi_tables, column=3,
    #     plotname='hist_energy_hi_all_new13co.png', table_list_shaded=robust_hi_tables)

    # print("Total for all shells:\n Mass - {} to {} Msun\n\
    #  Momentum - {} to {} Msun km/s\n\
    #  Energy - {} to {} erg\n".format(
    #     mass_all_low,mass_all_hi,
    #     momentum_all_low,momentum_all_hi,
    #     energy_all_low,energy_all_hi))


    # mass_robust_low = hist_physics(table_list=robust_low_tables, column=1,
    #     plotname='hist_mass_low_robust_new13co.png')
    # mass_robust_hi = hist_physics(table_list=robust_hi_tables, column=1,
    #     plotname='hist_mass_hi_robust_new13co.png')
    # momentum_robust_low = hist_physics(table_list=robust_low_tables, column=2,
    #     plotname='hist_momentum_low_robust_new13co.png')
    # momentum_robust_hi = hist_physics(table_list=robust_hi_tables, column=2,
    #     plotname='hist_momentum_hi_robust_new13co.png')
    # energy_robust_low = hist_physics(table_list=robust_low_tables, column=3,
    #     plotname='hist_energy_low_robust_new13co.png')
    # energy_robust_hi = hist_physics(table_list=robust_hi_tables, column=3,
    #     plotname='hist_energy_hi_robust_new13co.png')

    # print("Total for most confident 12 shells:\n Mass - {} to {} Msun\n\
    #  Momentum - {} to {} Msun km/s\n\
    #  Energy - {} to {} erg\n".format(
    #     mass_robust_low,mass_robust_hi,
    #     momentum_robust_low,momentum_robust_hi,
    #     energy_robust_low,energy_robust_hi))
    # cube_12co = SpectralCube.read(nro_12co)
    # cube_13co = SpectralCube.read(nro_13co)
    # #cube_13co = SpectralCube.read(nro_13co_divided)

    # shell_list = shells.get_shells()

    # for n in [19,33,41]:
    #     #shell = shell_list[n]
    #     shell = shell_list[n-1]


    #     params = np.loadtxt("shell_parameters_full.txt")
    #     params = params[params[:,0] == 1.*n, 1:][0]
    #     r_best, r_sig = params[0], params[1]
    #     dr_best, dr_sig = params[2], params[3]
    #     vexp_best, vexp_sig = params[4], params[5]
    #     v0_best, v0_sig = params[6], params[7]

    #     dv = 0.1
    #     v0_sample = np.arange(v0_best-v0_sig, v0_best+v0_sig+dv, dv)
    #     N = v0_sample.size
    #     print(N)

    #     ### Loop over several v0 values for the minimum r, dr, vexp.
    #     properties_low = np.empty((N, 4))
    #     r = (r_best - r_sig) * u.pc
    #     dr = (dr_best - dr_sig) * u.pc
    #     vexp = (vexp_best - vexp_sig) * u.km/u.s
        
    #     for i,v0 in enumerate(v0_sample):

    #         v0 = v0 * u.km/u.s   
    #         # v0 = 14.25*u.km/u.s 
    #         print(shell.ra, shell.dec, r, dr, vexp, v0) 
    #         try:
    #             mass, momentum, energy = calc_physics(
    #             ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    #             cube_12co=cube_12co, cube_13co=cube_13co)
    #             properties_low[i] = [v0.value, mass.value, momentum.value, energy.value]

    #             print("Shell Physical Properties:")
    #             print("------------------------------")
    #             print("Mass = {}".format(mass))
    #             print("Expansion Velocity = {}".format(vexp))
    #             print("Momentum = {}".format(momentum))
    #             print("Energy = {}".format(energy))
    #         except ValueError:
    #             print("Shell {} failed due to mismatched data shape.".format(n))
        
    #     np.savetxt("shell{}_properties_low.txt".format(n), properties_low)

    #     ### Loop over several v0 values for the mid r, dr, vexp.
    #     properties_mid = np.empty((N, 4))
    #     r = (r_best) * u.pc
    #     dr = (dr_best) * u.pc
    #     vexp = (vexp_best) * u.km/u.s
   
    #     for i,v0 in enumerate(v0_sample):

    #         v0 = v0 * u.km/u.s   
    #         # v0 = 14.25*u.km/u.s 
    #         print(shell.ra, shell.dec, r, dr, vexp, v0) 
    #         try:
    #             mass, momentum, energy = calc_physics(
    #             ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    #             cube_12co=cube_12co, cube_13co=cube_13co)
    #             properties_mid[i] = [v0.value, mass.value, momentum.value, energy.value]

    #             print("Shell Physical Properties:")
    #             print("------------------------------")
    #             print("Mass = {}".format(mass))
    #             print("Expansion Velocity = {}".format(vexp))
    #             print("Momentum = {}".format(momentum))
    #             print("Energy = {}".format(energy))
    #         except ValueError:
    #             print("Shell {} failed due to mismatched data shape.".format(n))
        
    #     np.savetxt("shell{}_properties_mid.txt".format(n), properties_mid)



    #     ### Loop over v0 values with the maximum r, dr, vexp.
    #     properties_hi = np.empty((N, 4))
    #     r = (r_best + r_sig) * u.pc
    #     dr = (dr_best + dr_sig) * u.pc
    #     vexp = (vexp_best + vexp_sig) * u.km/u.s
    #     for i,v0 in enumerate(v0_sample):

    #         v0 = v0 * u.km/u.s   
    #         print(shell.ra, shell.dec, r, dr, vexp, v0) 
    #         try:
    #             mass, momentum, energy = calc_physics(
    #             ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    #             cube_12co=cube_12co, cube_13co=cube_13co)
    #             properties_hi[i] = [v0.value, mass.value, momentum.value, energy.value]

    #             print("Shell Physical Properties:")
    #             print("------------------------------")
    #             print("Mass = {}".format(mass))
    #             print("Expansion Velocity = {}".format(vexp))
    #             print("Momentum = {}".format(momentum))
    #             print("Energy = {}".format(energy))
    #         except ValueError:
    #             print("Shell {} failed due to mismatched data shape.".format(n))
        


    #     np.savetxt("shell{}_properties_hi.txt".format(n), properties_hi)
    # # plt.figure()
    # plt.imshow(subcube_shell_12co_correct.moment0().data, interpolation='none')
    # plt.colorbar()
    # plt.title("Opacity-Corrected Integrated 12CO in K*m/s")
    # plt.savefig("corrected_12co.png")

def hist_physics(table_list=None, table_list_shaded=None, mode='median', column=1, plotname="hist_mass_low_all.png",
    return_total=True, bins='auto', xlim=None):
    """
    table_list_shaded gives shell property tables that I want to shade in the histogram.
    """
    #print(table_list)
    x = []
    x_shaded = []
    for t in table_list:
        if table_list_shaded and t in table_list_shaded:
            x_shaded_sample = np.loadtxt(t)[:,column]
            if mode == 'median':
                x_shaded = np.append(x_shaded, np.median(x_shaded_sample))
        x_sample = np.loadtxt(t)[:,column]
        if mode == 'median':
            x = np.append(x, np.median(x_sample))

    plt.figure()
    n, b, patches = plt.hist(
        [x, x_shaded], histtype='step',
        bins=bins,# stacked=True,
        label=["All shells", "Best 12 shells"],
        facecolor='black', edgecolor='black')
    hatches = ['', '']
    fills = [False,True]
    for patch_set, hatch, fill in zip(patches, hatches, fills):
        plt.setp(patch_set, hatch=hatch)
        plt.setp(patch_set, fill=fill)

    if column == 1:
        plt.xlabel(r"Mass [$M_\odot$]")
    if column == 2:
        plt.xlabel(r"Momentum [$M_\odot$ km/s]")
    if column == 3:
        plt.xlabel(r"Kinetic Energy [erg]")
    if xlim:
        plt.xlim(xlim)

    plt.ylabel("count")
    plt.legend()

    plt.savefig(plotname)

    if return_total:
        return np.sum(x)

def plot_physicsrange(low_name="_properties_low", mid_name="_properties_mid", hi_name="_properties_hi", name_tail=".txt",
    all_n=np.arange(1,44), best_n=best_shells, mode='median',
    column=1, plotname='massrange_all.png'):
    import matplotlib.lines as mlines
    plt.figure()

    for n in all_n:
        print(n, "shell{}{}{}".format(n, low_name, name_tail))
        low = np.loadtxt("shell{}{}{}".format(n, low_name, name_tail))[:,column]
        mid = np.loadtxt("shell{}{}{}".format(n, mid_name, name_tail))[:,column]
        hi = np.loadtxt("shell{}{}{}".format(n, hi_name, name_tail))[:,column]
        if mode == 'median':
            low = np.median(low)
            mid = np.median(mid)
            hi = np.median(hi)
        print(low,mid,hi)

        if n == 31:
            plt.annotate("No CO Shell", (0.3, 31),
             xycoords=('figure fraction', 'data'), verticalalignment='center')

        elif n in best_n:
            plt.plot([np.min([low,mid,hi]), np.max([low,mid,hi])], [n,n],
                color='k', ls='-')
            plt.plot(mid, n, marker='o', color='k')
            

        else:
            plt.plot([np.min([low,mid,hi]), np.max([low,mid,hi])], [n,n],
                color='k', ls=':')
            plt.plot(mid, n, marker='o', markerfacecolor='white', color='k')

            

    if column == 1:
        plt.xlabel(r"Mass [$M_\odot$]")
    if column == 2:
        plt.xlabel(r"Momentum [$M_\odot$ km/s]")
    if column == 3:
        plt.xlabel(r"Kinetic Energy [erg]")

    plt.ylabel("Shell Number")
    best_line = mlines.Line2D([], [],
     color='k', marker='o', linestyle='solid', label='Best 12 shells')
    other_line = mlines.Line2D([], [],
     color='k',marker='o', markerfacecolor='white', linestyle='dashed', label='Other shells')
    plt.legend(handles=[best_line, other_line], loc='best')
    plt.savefig(plotname, dpi=300)



def calc_physics(ra=None, dec=None, r=0.17*u.pc, dr=0.05*u.pc,
 vexp=4*u.km/u.s, v0=14*u.km/u.s, cube_12co=None, cube_13co=None,
 dist=414*u.pc, snr_cutoff=5.):
    #print(cube_12co.header['CDELT3'])
    pix_size = (cube_12co.header['CDELT2']*u.deg).to(u.arcsec)
    vstep = (cube_12co.header['CDELT3']*(u.m/u.s)).to(u.km/u.s)
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
    print(subcube_shell_12co, subcube_shell_13co, model_cube)
    ### Excitation Temperature
    Tex = cube_Tex(subcube_shell_12co, average_first=True)
    
    ### Correct 12co for opacity.
    subcube_shell_12co_correct = opacity_correct(
        subcube_shell_12co, cube_thin=subcube_shell_13co,
        snr_cutoff=snr_cutoff)
    print(subcube_shell_12co_correct)
    ### Extract shell voxels from opacity-corrected 12co
    shell_12co_correct = extract_shell(
        subcube_shell_12co_correct, keep_latlon=True, model_cube=model_cube)

    ### Calculate column density of H2 from opacity-corrected 12co
    shell_nH2 = column_density_H2(shell_12co_correct, Tex=Tex)

    ### Calculate Mass, Momentum, and Energy of Shell!
    shell_mass = mass(shell_nH2, distance=414*u.pc, molecule='H2',
     mass_unit=u.Msun)
    shell_momentum = momentum(shell_mass, vexp)
    shell_energy = energy(shell_mass, vexp)

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


def extract_shell(cube_file=nro_12co, model_pars=None,
 mask_minimum=0.00001, return_mask=False, model_cube=None,
 keep_latlon=False):
    """
    Return a masked cube from an observed spectral cube file,
    where the mask is True wherever a model shell cube
    with parameters given by `model_pars` dictionary is > mask_minimum.

    """
    if not model_cube:
        model_cube = SpectralCube.read(shell_model.ppv_model(**model_pars))
    if type(cube_file) == str:
        if keep_latlon:
            obs_cube = SpectralCube.read(cube_file).spectral_slab(
                model_cube.spectral_extrema[0],
                model_cube.spectral_extrema[1])
        else:
            obs_cube = SpectralCube.read(cube_file).subcube(
                                model_cube.longitude_extrema[1],
                                model_cube.longitude_extrema[0],
                                model_cube.latitude_extrema[0],
                                model_cube.latitude_extrema[1],
                                model_cube.spectral_extrema[0],
                                model_cube.spectral_extrema[1])
    else:
        if keep_latlon:
            obs_cube = cube_file.spectral_slab(
                model_cube.spectral_extrema[0],
                model_cube.spectral_extrema[1])
        else:
            obs_cube = cube_file.subcube(
                                model_cube.longitude_extrema[1],
                                model_cube.longitude_extrema[0],
                                model_cube.latitude_extrema[0],
                                model_cube.latitude_extrema[1],
                                model_cube.spectral_extrema[0],
                                model_cube.spectral_extrema[1])
    print("Before changing wcs, obs_cube shape: ", obs_cube.shape)
    #Reset the cube wcs to the values corresponding to the subcube.
    obs_cube = SpectralCube(obs_cube.hdu.data, wcs=model_cube.wcs) * u.K
    shell_mask = model_cube > mask_minimum*u.dimensionless_unscaled
    print(model_cube.shape, obs_cube.shape)
    obs_cube_masked = obs_cube.with_mask(shell_mask)
    #obs_array_masked = obs_cube_masked.filled_data[:,:,:]
    if return_mask:
        return obs_cube_masked, shell_mask
    else:
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

def average_spectrum(cube=None, weights=1., axis=(1,2), return_std=True,
    ignore_nan=True):
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
    weighted_data = cube.filled_data[:,:,:] * weights
    average_spectrum = np.nanmean(weighted_data, axis=axis)
    #average_spectrum = np.average(cube.filled_data[:,:,:],
    # weights=weights, axis=axis)

    if return_std:
        resids_squared = (cube.filled_data[:,:,:] - average_spectrum[:,np.newaxis,np.newaxis])**2. * weights
        std_spectrum = np.sqrt(
            np.nanmean(resids_squared, axis=axis))
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
            X_factor = 1e-4 # CO/H2

    factor = (
        (8*np.pi*const.k_B*nu_ul**2.)\
        / (const.h*const.c**3.*A_ul*g_u)
             )

    factor *= Qrot_partial(Tex, B0_k, N=Qrot_order)
    factor *= np.exp(E_u_k/Tex) / beam_filling_factor

    ### APPLY THE ABUNDANCE RATIO BETWEEN CO/H2. CHECK THIS
    factor /= X_factor

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
    Tex = u.Quantity(Tex, u.K)
    Qrot = 0.
    for J in range(N+1):
        Qrot += (2*J + 1) * np.exp(
            -B0_k*J*(J+1)/Tex)
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




def momentum(mass, velocity, unit=u.Msun*(u.km/u.s)):
    return (mass * velocity).to(unit)
def energy(mass, velocity, unit=u.erg):
    return (0.5 * mass * velocity ** 2.).to(unit)


if __name__ == '__main__':
    main()