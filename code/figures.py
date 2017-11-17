from spectral_cube import SpectralCube
import pyregion
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.wcs import WCS
#from wcsaxes import WCSAxes
import aplpy
from aplpy import FITSFigure
import numpy as np
from astroquery.simbad import Simbad
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib
import astropy.units as u

import shells
import physics
import shell_model
from astropy.io import ascii
from astropy.io import fits

orion_dist = 414*u.pc #pc
nro_12co = "../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS_regrid0.11kms_reproj.fits"
nro_13co = "../nro_maps/13CO_BEARS-FOREST_20170913_7.5grid_Spheroidal_Tmb_0.11kms_xy_YS.fits" 
best_shells = [3,6,9,11,17,18,21,24,25,30,36,37]
north_shells = [18,19,20,21,22,23,24,29,40]
central_shells = [16,17,26,30,36,38,39]
south_shells = [3,4,5,6,7,15,28,33,34,35]
l1641_shells = [1,2,8,9,10,11,12,13,14,25,27,31,32,37,41,42]

mips_l1641_file = '../catalogs/MIPS_L1641a_24um.fits'
mips_onc_file = '../catalogs/MIPS_ONC_24um.fits'

irac1_l1641_file = '../catalogs/IRAC_L1641_ch1_merged_clean.fits'
irac1_onc_file = '../catalogs/IRAC_ONC_ch1_merged_clean.fits'

irac2_l1641_file = '../catalogs/IRAC_L1641_ch2_merged_clean.fits'
irac2_onc_file = '../catalogs/IRAC_ONC_ch2_merged_clean.fits'

irac4_l1641_file = '../catalogs/IRAC_L1641_ch4_merged_clean_northup.fits'
irac4_onc_file = '../catalogs/IRAC_ONC_ch4_merged_clean_northup.fits'

planck_herschel_file = '../catalogs/planck_herschel.fits'

region_file = '../shell_candidates/AllShells.reg'
vrange_file = '../shell_candidates/AllShells_vrange.txt'
parameter_file = "shell_parameters_full.txt"

obaf_file = 'stars_obaf.txt'
yso_file = "../catalogs/spitzer_orion.fit"
hops_file = "../catalogs/hops.fits"

hst_rgbpars = dict([
['stretch_r','arcsinh'],
['stretch_g','arcsinh'],
['stretch_b','arcsinh'],
['vmin_r',0.03],
['vmax_r',4.55],
['vmin_g',0.03],
['vmax_g',5.35],
['vmin_b',0.03],
['vmax_b',3.08]])

def main():
    
    # plot_overview(plotname="../paper/figs/12co_nroonly_peak_full_shells.png", show_shells=True,
 #     dist=orion_dist, vmin=None, vmax=None, scalebar_color='black', scale_factor = 1.,
 #     title=r"", shells_highlight=best_shells, circle_style='dotted', circle_linewidth=1.5,
 #     scalebar_pc=2., tick_color="gray" #,recenter=False, ra=83.99191, dec=-5.6611303, radius=0.117325
 #     )

    shell_list = shells.get_shells(region_file=region_file, velocity_file=vrange_file)
    shell_parameters = ascii.read(parameter_file) 
    obaf = ascii.read(obaf_file)
    obaf_ra, obaf_dec, obaf_id, obaf_sptype = np.array(obaf['RA']), np.array(obaf['DEC']),\
     np.array([" ".join([word for word in line.split() if '*' not in word]).strip("b'") for line in obaf['MAIN_ID']]),\
     np.array([sp.strip("b'") for sp in obaf['SP_TYPE']])

    yso = fits.open(yso_file)[1].data
    yso_ra, yso_dec, yso_label = yso['RAJ2000'], yso['DEJ2000'], yso['Cl']

    hops = fits.open(hops_file)[1].data
    hops_ra, hops_dec, hops_label = hops['RAJ2000'], hops['DEJ2000'], hops['Class']


    ir_l1641_hdu = fits.open(irac4_l1641_file)[0]
    ir_onc_hdu = fits.open(irac4_onc_file)[0]

    #spec_cube = SpectralCube.read(nro_13co)
    spec_cube = SpectralCube.read(nro_12co)
    molecule_name = "12CO"
    low_sigma, hi_sigma, step_sigma = 10, 22, 2

    ra_grid = spec_cube.spatial_coordinate_map[1].to(u.deg).value
    dec_grid = spec_cube.spatial_coordinate_map[0].to(u.deg).value
    vel_grid = spec_cube.spectral_axis
    pad_factor = 1.8

    #plot_overview(show_shells=True)
    #plot_overview(plotname="12co_nro_peak.png", show_shells=False)
    #plot_overview(cube="/Volumes/Untitled/13co_pix_2.cm.fits", plotname="13co_combined_peak.png", show_shells=False)
    #return
    #channel_vmax = [12.9, 14]
    for nshell in [38]:
        shell = shell_list[nshell-1]
        #ra, dec, radius = shell.ra.value, shell.dec.value, shell.radius.value
        ra, dec, radius = shell.ra.value, shell.dec.value, ((206265. / 3600.) * shell_parameters['r'][nshell-1]*u.pc / orion_dist).value
        thickness = ((206265. / 3600.) * shell_parameters['dr'][nshell-1]*u.pc / orion_dist).value

        l1641_xy = WCS(ir_l1641_hdu).wcs_world2pix(ra, dec, 0)
        #print(l1641_xy)
        if  (l1641_xy[0] >= 0) & (l1641_xy[0] <= ir_l1641_hdu.shape[1]) & \
           (l1641_xy[1] >= 0) & (l1641_xy[1] <= ir_l1641_hdu.shape[0]):
            ir_hdu = ir_l1641_hdu
        else:
            ir_hdu = ir_onc_hdu

        #Extract sub_cube around shell.
        subcube_mask = (abs(ra_grid - ra) < radius * pad_factor) &\
               (abs(dec_grid - dec) < radius * pad_factor)
        sub_cube = spec_cube.with_mask(subcube_mask).minimal_subcube()
        sigma = np.nanmedian(physics.rms_map(sub_cube).value)
        #shell.vmin, shell.vmax = 12.*u.km/u.s, 13.4*u.km/u.s
        sub_cube = sub_cube.spectral_slab(shell.vmin, shell.vmax)

        #Integrate between vmin and vmax.
        peak_hdu = sub_cube.max(axis=0).hdu
        mom0_hdu = sub_cube.moment0().hdu

        mask_inshell = (abs(ra_grid - ra) < radius) &\
               (abs(dec_grid - dec) < radius)
        subcube_inshell = spec_cube.with_mask(mask_inshell).minimal_subcube().spectral_slab(shell.vmin, shell.vmax)
        mom0_hdu_inshell = subcube_inshell.moment0().hdu

        #Calculate contour levels.
        empty_channel = spec_cube.closest_spectral_channel(500*u.Unit('m/s'))
        # sigma = np.nanstd(spec_cube[empty_channel][subcube_mask]).value
        #print("sigma: {}".format(sigma))
        delta_vel = (sub_cube.spectral_extrema[1] - sub_cube.spectral_extrema[0]).value
        #print("delta_vel: {}".format(delta_vel))
        mom0_sigma = sigma * delta_vel
        print(sigma, mom0_sigma)     
        #contour_levels = np.linspace(5.*mom0_sigma, np.nanmax(mom0_hdu_inshell.data), 12)

        # Shell 17 12co mom0 contours:
        # 
        #
        #

        #Shell 17 13co mom0 contours:
        # low_sigma, hi_sigma, step_sigma = 10, 22, 2

        #Shell 37 12co mom0 contours:
        # low_sigma, hi_sigma, step_sigma = 48, 60, 3

        #Shell 37 13co mom0 contours:
        # low_sigma, hi_sigma, step_sigma = 20.01, 40.01, 4

        #Shell 42 13co mom0 contours:
        #low_sigma, hi_sigma, step_sigma = 3, 8, 1

        #Shell 42 12co mom0 contours:
        low_sigma, hi_sigma, step_sigma = 35, 75, 8     

        contour_levels = np.arange(low_sigma*mom0_sigma, (hi_sigma+1)*mom0_sigma, step_sigma*mom0_sigma)
        #13co contours:
        #contour_levels = np.linspace(20*sigma, 100*sigma, 8)
        print("Contours drawn at {} or {} sigma.".format(contour_levels, contour_levels / mom0_sigma))

        #Get source coordinates.
        print(obaf_ra, obaf_dec, radius)
        obaf_in_shell = ((obaf_ra - ra) ** 2. + (obaf_dec - dec) ** 2.) <= (radius*pad_factor) ** 2.
        print(obaf_in_shell)
        yso_in_shell = ((yso_ra - ra) ** 2. + (yso_dec - dec) ** 2.) <= (radius) ** 2.

        yso_ra, yso_dec, yso_label = yso_ra[yso_in_shell], yso_dec[yso_in_shell], yso_label[yso_in_shell]
        obaf_ra, obaf_dec, obaf_id = obaf_ra[obaf_in_shell], obaf_dec[obaf_in_shell], obaf_id[obaf_in_shell]
        print(yso_ra, yso_dec, yso_label)
        source_ra, source_dec, source_labels = [obaf_ra, yso_ra], [obaf_dec, yso_dec], [obaf_id, None]
        print(source_ra, source_dec, source_labels)
        # source_ra, source_dec, source_labels = obaf_ra, obaf_dec, obaf_id
        #source_labels[~obaf_in_shell] = "" 
        print(ra, dec, radius)
        plot_stamp(map=ir_hdu,
            ra=ra, dec=dec, radius=radius,
            thickness=thickness, rgb=None,
            circle_color='red', scalebar_pc=0.1, nan_color='black',
            pad_factor=pad_factor, contour_map=mom0_hdu, contour_levels=contour_levels, contour_color='white',
            plotname='../paper/figs/mom0{}shell{}_{}{}to{}_stamp.png'.format('8micron', nshell, molecule_name, str(shell.vmin.value).replace('.','_'), str(shell.vmax.value).replace('.','_')),
            return_fig=False, vmin=0, vmax=70000,
            stretch='linear', plot_simbad_sources=False, dist=orion_dist,
            auto_scale=False, auto_scale_mode='median', auto_scale_pad_factor=0.8, auto_scale_nsigma=5.,
            cbar_label="Counts", cmap='inferno', show_colorbar=False,
            # source_ra=[obaf_ra, yso_ra], source_dec=[obaf_dec, yso_dec],
            source_colors=['white', 'cyan'], source_markers=['*', '+'], source_sizes=[300,50],
            source_edge_colors=['black','cyan'],
            # source_colors='white', source_markers='*', source_sizes=300,
            # source_labels=[obaf_id, None], dpi=200
            source_ra=source_ra, source_dec=source_dec, 
            source_labels=source_labels, label_offset=[-0.005, -0.004], dpi=200
            )

        #cube_file = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
        
        
        # plot_channels(cube=nro_12co, ra=ra, dec=dec, radius=radius,
        #     source_lists=None, stretch='linear', pad_factor=pad_factor,
        #     vel_min=11*u.km/u.s, vel_max=13.5*u.km/u.s,
        #     # vel_min=shell.vmin.value, vel_max=shell.vmax.value,
        #     plotname='../paper/figs/12co_channels_shell'+str(nshell)+'.png', chan_step=1, plot_simbad_sources=False,
        #     vmin=None, vmax=None, max_chans=12, cmap='gist_yarg', text_color='red',
        #     xspacing_in_s=10.,
        #     #cbar_label="Counts",
        #     source_ra=[obaf_ra, yso_ra], source_dec=[obaf_dec, yso_dec],
        #     source_colors=['white', 'cyan'], source_edge_colors=['black','cyan'],
        #      source_markers=['*', '+'],
        #     source_sizes=[160,15], dpi=300)
        shell_list = shells.get_shells()
        shell = shell_list[nshell-1]

        model_pars = {
            'dist':414*u.pc, # pc
            'pix_size':7.5*u.arcsec, # arcsec
            'vstep':0.110*u.km/u.s, # km/s
            'acen':shell.ra, # deg
            'dcen':shell.dec, # deg
            'R':0.22*u.pc, # pc
            'dr':0.2*u.pc, # pc
            'vexp':2.5*u.km/u.s, # km/s
            'v0':13.6*u.km/u.s, # km/s
            'samples_per_voxel':27}
        
        model_pars['R'] = shell_parameters['r'][nshell-1] * u.pc
        # model_pars['R'] =  0.05 * u.pc
        model_pars['dr'] = shell_parameters['dr'][nshell-1]* u.pc
        # model_pars['dr'] = 0.1 * u.pc
        model_pars['vexp'] = shell_parameters['vexp'][nshell-1] * u.km / u.s
        model_pars['v0'] = shell_parameters['vsys'][nshell-1] * u.km / u.s

        # fig = plot_pv(plot_file='../paper/figs/12co_pvaverage_shell'+str(nshell)+'.png',
        #     cube_file=nro_12co,
        #     model_pars=model_pars, normalize=False,
        #     contour_levels=20, pv_length_in_radii=5,
        #     pv_width_in_pixels=3., pv_angle=0*u.deg,
        #     pv_angle_step=90*u.deg, average_pv_obs=True,
        #     remove_first_contour=True)
        # pv = plot_pv(cube=cube_file, ra_center=shell.ra, dec_center=shell.dec,
        #      vel=[shell.vmin - 1*u.km/u.s, shell.vmax + 1*u.km/u.s], length=shell.radius*4.,
        #      width=7.5*u.arcsec, angle=angle,
        #      pad_factor=1., plotname='12co_pv_shell'+str(nshell)+'_angle'+str(angle.value)+'.png',
        #      stretch='linear', auto_scale=True, dpi=900.)


def plot_channels(cube=None, shellShape=None, ra=None, dec=None, radius=None, circle_color='white',
    pad_factor=1., contour_map=None, source_ra=None, source_dec=None, source_lists=None, source_sizes=None,
    source_markers=None, cmap='viridis', text_color='red', xspacing_in_s=None,
    source_ra_colnames='RAJ2000', tick_color = 'gray',
    source_dec_colnames='DEJ2000', source_colors='blue', source_edge_colors='blue',
     plotname='shell_stamp.png', return_subplot=True,
    stretch='linear', vmin=20, vmax=50, plot_simbad_sources=True, simbad_type='star', simbad_color='green',
    vel_min=0, vel_max=10, vel_unit=u.km/u.s, n_cols=3, max_chans=16,
    chan_step=1, auto_scale=True, dpi=300):
    """
    Parameters
    ----------
    cube : SpectralCube or str, optional
        The cube to be used, either a fits file str of a SpectralCube object.
    shellShape : pyregion.Shape, optional
        Shape object denoting the shell's center and radius.
    ra : float, optional
        If a Shape object is not specified, use ra/dec/radius in degrees.
    dec : float, optional
        Declination of shell in degrees
    radius : float, optional
        Radius optionalf shell in degrees.
    circle_color : str, optional
        Description
    pad_factor : float, optional
        Size of stamp to plot, where size = radius * pad_factor
    contour_map : same as map, optional
        Optionally plot contours on top of the primary map.
    source_ra : list or ndarray, optional
        Source ra in degrees
    source_dec : list or ndarray, optional
        Source dec in degrees
    source_lists : table or list of tables
        Can be fits table, or HDU, which specifies the source catalogs to
        use for plotting sources.
    source_ra_colnames : str, optional
        Description
    source_dec_colnames : str, optional
        Description
    source_colors : str, optional
        Description
    plotname : str, optional
        Plot to write out.
    return_subplot : bool, optional
        Return a subplot that can be used in multipanel plots.
    stretch : str, optional
        Description
    vmin : int, optional
        Description
    vmax : int, optional
        Description
    plot_simbad_sources : bool, optional
        Description
    simbad_type : str, optional
        Description
    simbad_color : str, optional
        Description
    vel_min : float or Quantity, optional
        Description
    vel_max : float or Quantity, optional
        Description
    vel_unit : Unit or str, optional
        Description
    n_cols : int, optional
        Description
    chan_step : int, optional
        Description
    auto_scale : bool, optional
        Description
    
    Raises
    ------
    Exception
        Description
    
    """
    if plot_simbad_sources:
        print("Finding simbad sources of type {} within {} of ra {} and dec {}.".format(
            simbad_type, radius, ra, dec))
        simbad_table = simbad_sources(ra, dec, radius, unit='deg')
        simbad_coords = coord.SkyCoord(simbad_table['RA'], simbad_table['DEC'],
            unit=(u.hourangle, u.deg))

    try:
       # if cube is not already a SpectralCube, but is a fits file.
        spec_cube = SpectralCube.read(cube)
    except ValueError:
        # `spec_cube` is a SpectralCube object
        cube = cube.hdu
    #Find the min and max spectral channel.
    try:
        chan_min = spec_cube.closest_spectral_channel(vel_min)
        chan_max = spec_cube.closest_spectral_channel(vel_max)
    except AttributeError:
        try:
            chan_min = spec_cube.closest_spectral_channel(vel_min * vel_unit)
            chan_max = spec_cube.closest_spectral_channel(vel_max * vel_unit)
        except AttributeError:
            chan_min = spec_cube.closest_spectral_channel(vel_min * u.Unit(vel_unit))
            chan_max = spec_cube.closest_spectral_channel(vel_max * u.Unit(vel_unit))

    #n_chan = chan_max - chan_min + 1
    n_chan = len(np.arange(chan_min, chan_max + 1, chan_step))
    
    while n_chan > max_chans:
        chan_step += 1
        n_chan = len(np.arange(chan_min, chan_max + 1, chan_step))

    # The number of rows necessary to fit all n channels.
    n_cols = 4
    n_rows =  np.ceil(n_chan / n_cols)


    #For the auto color scaling to the min and max intensities inside the shell.
    ra_grid = spec_cube.spatial_coordinate_map[1].to(u.deg).value
    dec_grid = spec_cube.spatial_coordinate_map[0].to(u.deg).value
    shell_mask = (ra_grid - ra) ** 2. + (dec_grid - dec) ** 2. < radius ** 2.

    #Create a figure and add subplots.
    fig = plt.figure(figsize=(8, 8*n_rows/n_cols))
    for ith_subplot, chan in enumerate(range(chan_min, chan_max + 1, chan_step)):
        print("Plotting channel {} out of {}.".format(chan, n_chan))
        subplot = FITSFigure(cube, figure=fig,
            slices=[chan],
            subplot=(n_rows, n_cols, ith_subplot+1))
        #Establish whether this subplot is on the left edge or bottom.
        #test True if this subplot on the left edge.
        left_edge = (ith_subplot % n_cols) == 0 
        #test True if a subplot will not appear below this one.
        bottom_edge = (ith_subplot + n_cols + 1) > n_chan

        subplot.recenter(ra, dec, radius*pad_factor) #Pan/zoom to shell 
        subplot.tick_labels.set_yformat("dd:mm")
        if xspacing_in_s:
            subplot.ticks.set_xspacing(xspacing_in_s*15/3600.)
            


        subplot.tick_labels.set_xformat("hh:mm:ss")
        subplot.tick_labels.set_style('plain')
        if ith_subplot != 0:
            #Hide all Y-labels except for first, left edge.
            subplot.hide_ytick_labels()
            subplot.hide_yaxis_label()

        # if not left_edge:
        #     subplot.hide_ytick_labels()
        #     subplot.hide_yaxis_label()
        if ith_subplot != n_cols * (n_rows - 1):
            #Hide all X-labels except for bottom left corner.
            subplot.hide_xtick_labels()
            subplot.hide_xaxis_label()

        # if not bottom_edge:
        #     subplot.hide_xtick_labels()
        #     subplot.hide_xaxis_label()


        #AUTO SCALING
        if auto_scale:
            shell_pixels = spec_cube[chan][shell_mask].value
            vmin, vmax = np.nanmin(shell_pixels), np.nanmax(shell_pixels)
        subplot.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap=cmap)     

        #COLORBAR
        #subplot.add_colorbar()
        #cb = subplot.colorbar
        #cb.set_axis_label_text(r'T$_{MB}$ [K]')

        #SHELL OUTLINE
        subplot.show_circles(ra, dec, radius, linestyle='dashed', edgecolor=circle_color,
            facecolor='none')

            #POINT SOURCES
        if source_ra and source_dec:
            #Use source_ra/dec if these parameters are set.
            try:
                subplot.show_markers(source_ra, source_dec,
                 c=source_colors, marker=source_markers, s=source_sizes,
                 edgecolor=source_edge_colors)
            except TypeError:
                #If more than one source list to be plotted with different markers
                #source_ra, source_dec, source_colors must be nested lists or ndarrays 
                #with same shape.
                #print(source_ra, source_dec, source_colors, source_markers, source_sizes)
                for i in range(len(source_colors)):
                    subplot.show_markers(source_ra[i], source_dec[i], c=source_colors[i],
                     marker=source_markers[i], s=source_sizes[i],
                     edgecolor=source_edge_colors[i])

        elif source_lists and type(source_lists) is not list:
            try:
                #source_lists is a single fits file string
                source_hdu = fits.open(source_lists)[-1] #Select last HDU in HDUList
            except OSError:
                #source_lists is a single hdu object
                source_hdu = source_lists

            source_ra = source_hdu.data[source_ra_colnames]
            source_dec = source_hdu.data[source_dec_colnames]
            subplot.show_markers(source_ra, source_dec, edgecolor=source_edge_colors[i])

        elif source_lists and type(source_lists) is list:
            raise Exception("source_lists of type(list) not implemented.")
            if type(source_lists[0]) is str:
            # source_lists is a list of fits file strings
                raise Exception("Souce_lists is a list of fits file strings not implemented.")
            else:
                # source_lists is a list of hdu objects
                raise Exception("Souce_lists is a list of hdu objects not implemented.")

        if plot_simbad_sources:
            simbad_table = simbad_sources(ra, dec, radius, unit='deg')
            simbad_coords = coord.SkyCoord(simbad_table['RA'], simbad_table['DEC'],
                unit=(u.hourangle, u.deg))

            subplot.show_markers(simbad_coords.ra, simbad_coords.dec, edgecolor=simbad_color)

            ### LABEL THE VELOCITY OF EACH CHANNEL MAP
        
        subplot.add_label(0.6, 0.9,
            str(np.round(spec_cube.spectral_axis[chan].to('km/s'), 1).value)+' km/s',
            color=text_color, relative=True)
        subplot.ticks.set_color(tick_color)
        subplot.ticks._ax1.tick_params(direction='in', which='both', axis='both')
        subplot.ticks._ax2.tick_params(direction='in', which='both', axis='both')

        #subplot.close()
    #Make the figure prettier.
    #fig.tight_layout(h_pad=0, w_pad=0)
    frac = 0.08
    fig.subplots_adjust(bottom=frac, left=frac * n_rows / n_cols, 
        top=1 - frac, right=1 - (frac * n_rows / n_cols),
        wspace=0., hspace=0.)
    #fig.subplots_adjust(wspace=None, hspace=None)
    fig.canvas.draw()
    # fig.ticks.set_color(tick_color)
    fig.savefig(plotname, dpi=dpi, bbox_inches='tight')
#    fig.close()
   

def plot_pv(plot_file=None, cube_file="../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS_regrid0.11kms_reproj.fits",
    regionfile="../shell_candidates/AllShells.reg", show_model=True,
    pv_width_in_pixels=3., pv_length_in_radii=3., pv_angle=0*u.deg, pv_angle_step=90*u.deg,
    average_pv_obs=True, average_pv_model=False, model_pars=None, cmap='gist_yarg',
    normalize=True, contour_levels=20, draw_radius=True, dist=orion_dist, pixel=7.5*u.arcsec,
    fontsize=16, remove_first_contour=False, vmin=None, vmax=None, pmin=0.25, pmax=99.75):
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
    from matplotlib.ticker import FuncFormatter
    radius_in_pixels = (model_pars['R'] / dist) * (206265 * u.arcsec) / pixel
    pad_pixels = (pv_length_in_radii - 2) * radius_in_pixels / 2

    try:
        cube_model = SpectralCube.read(shell_model.ppv_model(pad_pixels=pad_pixels, **model_pars))
    except:
        cube_model = SpectralCube.read(shell_model.ppv_model(pad_pixels=pad_pixels))

    head_model = cube_model.header
    cube_obs = SpectralCube.read(cube_file).subcube(
                                cube_model.longitude_extrema[1],
                                cube_model.longitude_extrema[0],
                                cube_model.latitude_extrema[0],
                                cube_model.latitude_extrema[1],
                                cube_model.spectral_extrema[0],
                                cube_model.spectral_extrema[1])
    cube_obs = cube_obs.with_spectral_unit('km/s')
    cube_model = cube_model.with_spectral_unit('km/s')
    print(cube_obs.header, cube_model.header)
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

    if show_model:
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

    # fig = plt.figure(figsize=(20,10))

    # title = """angle = {}, length = {} , width = {},
    #     r = {} pc, dr = {} pc, v0 = {} km/s, vexp = {} km/s""".format(
    #     pv_angle.round(2), pv_length.to(u.arcmin).round(2), pv_width.to(u.arcsec).round(2),
    #     head_model['r'], head_model['dr'], head_model['v0'], head_model['vexp'])

    fig = FITSFigure(pv_obs)
    fig.show_colorscale(aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, pmin=pmin, pmax=pmax)

    #pv_model.data /= np.nanmax(pv_model.data)
    if show_model:
        fig.show_contour(pv_model, levels=contour_levels, layer='contour')
    if remove_first_contour:
        c = fig.get_layer('contour')
        c.collections[0].remove()
        plt.draw()


    #fig.add_colorbar()

    # figpv.set_title(title)

    # figmom0 = FITSFigure(cube_obs.moment0().hdu, figure=fig, subplot=((1,2,2)))
    # figmom0.show_grayscale()
    # #figmom0.show_contour(cube_model.moment0().hdu, levels=int(contour_levels/2))
    # figmom0.set_title(title)

    # if draw_radius:
    #     r_degrees = (head_model['R'] / head_model['DIST']) * 360. / (2 * np.pi) 
    #     dr_degrees = (head_model['DR'] / head_model['DIST']) * 360. / (2 * np.pi)

    #     figmom0.show_circles(ra_center.value, dec_center.value, r_degrees)
    #     figmom0.show_circles(ra_center.value, dec_center.value, r_degrees - dr_degrees / 2.,
    #      linestyle='--', edgecolor='red')
    #     figmom0.show_circles(ra_center.value, dec_center.value, r_degrees + dr_degrees / 2.,
    #      linestyle='--', edgecolor='red')
    plt.ylabel("Velocity [km/s]", size=fontsize)
    plt.xlabel("Offset [deg]", size=fontsize)
    dv_kms = cube_obs.spectral_axis[1].to(u.km/u.s).value - cube_obs.spectral_axis[0].to(u.km/u.s).value
    v0_kms = cube_obs.spectral_axis[0].to(u.km/u.s).value
    print(dv_kms, v0_kms)

    FormatKMS = FuncFormatter(lambda x,y: '{0:g}'.format(round(v0_kms + (x-1)*dv_kms,5)))
    fig.image.axes.yaxis.set_major_formatter(FormatKMS)

    if plot_file:
        fig.save(plot_file)

    return fig



def plot_overview(cube=nro_12co, subregion_file="../subregions/my_subregions.reg",
 region_file='../shell_candidates/AllShells.reg', mode='peak', plotname='12co_peak_shells.png',
 interactive=False, show_shells=True, shells_highlight=None, dist=orion_dist, vmin=None, vmax=None,
 scalebar_color="white", scalebar_pc = 1., scale_factor=1., pmin=0.25,
 pmax=99.75, cbar_label=r"Peak T$_{\rm MB}$ [K]", tick_color="white",
 circle_color='white', circle_linewidth=1, circle_style="solid", return_fig=False, show=True,
 title=r"$^{12}$CO Peak T$_{MB}$", recenter=False, ra=None, dec=None, radius=None, show_subregions=True,
 subregion_color='red', subregion_linewidth=1, subregion_style='dashed'):
    """
    Show full image with all shells.
    
    Parameters
    ----------
    cube : str, optional
        Description
    region_file : str, optional
        Description
    mode : str, optional
        Description
    plotname : str, optional
        Description
    interactive : bool, optional
        Description
    show_shells : bool, optional
        Description
    
    """
    try:
        cube = SpectralCube.read(cube)
    except ValueError:
        pass

    if mode == "peak":
        image = (cube.max(axis=0) * scale_factor).hdu
        

    if mode == "mom0":
        image = (cube.moment0() * scale_factor).hdu



    fig = FITSFigure(image)
    if show:
        fig.show_colorscale(cmap='viridis', vmin=vmin, vmax=vmax, pmin=pmin,
                pmax=pmax, interpolation='none')
    fig.tick_labels.set_yformat("dd:mm")
    fig.tick_labels.set_xformat("hh:mm")

    fig.ticks.set_color(tick_color)
    #fig.hide_yaxis_label()
    #fig.hide_ytick_labels()
    plt.title(title)
    plt.xlabel("RA (J2000)")
    plt.ylabel("DEC (J2000)")

    if show_shells:
        shell_list = get_shells(region_file=region_file)
        for i, shell in enumerate(shell_list):
            if shells_highlight:
                if i+1 in shells_highlight:
                    fig.show_circles(shell.ra.value, shell.dec.value, shell.radius.value, linestyle='solid', edgecolor=circle_color,
                        facecolor='none', linewidth=3)
                else:
                    fig.show_circles(shell.ra.value, shell.dec.value, shell.radius.value, linestyle=circle_style, edgecolor=circle_color,
                        facecolor='none', linewidth=circle_linewidth)
            else:
                fig.show_circles(shell.ra.value, shell.dec.value, shell.radius.value, linestyle=circle_style, edgecolor=circle_color,
                    facecolor='none', linewidth=circle_linewidth)

    if show_subregions:
        region_list = pyregion.open(subregion_file)
        label_list = [r.comment for r in region_list]
        
        for region in region_list:

            reg_ra, reg_dec, reg_width, reg_height, reg_number = region.coord_list
            label = region.comment
            if label == "North":
                vertices = np.array([[reg_ra + reg_width/2., reg_ra - reg_width/2.],
                                     [reg_dec - reg_height/2., reg_dec - reg_height/2.]])
            if label == "Central":
                vertices = np.array([[reg_ra + reg_width/2., reg_ra - reg_width/2.],
                                     [reg_dec - reg_height/2., reg_dec - reg_height/2.]])
            if label == "South":
                l1641n_ra = region_list[label_list.index('L1641N')].coord_list[0]
                l1641n_width = region_list[label_list.index('L1641N')].coord_list[2]
                vertices = np.array([[reg_ra + reg_width/2., l1641n_ra - l1641n_width/2.],
                                     [reg_dec - reg_height/2., reg_dec - reg_height/2.]])
            fig.add_label(84.2500, reg_dec+reg_height/2. - 0.06, label, color=subregion_color)
            fig.show_lines([vertices],
             color=subregion_color, linewidth=subregion_linewidth,
             linestyle=subregion_style)



    #RECENTER
    if recenter:
        fig.recenter(ra, dec, radius)


    #SCALEBAR
    fig.add_scalebar(206265 * scalebar_pc / (dist.to(u.pc).value * 3600), color=scalebar_color,
        )
    fig.scalebar.set_label("{} pc".format(scalebar_pc))

    fig.add_colorbar()
    cb = fig.colorbar
    cb.set_axis_label_text(cbar_label)

    if return_fig:
        return fig
    else:
        fig.save(plotname, dpi=600)

def plot_stamp(map=None, fig=None, shell=None, ra=None, dec=None, radius=None, thickness=None,
    circle_color='red', rgb=None, scalebar_pc=0.2,
    pad_factor=1.5, contour_map=None, contour_levels=5, contour_color='white',
    source_ra=None, source_dec=None, source_lists=None, source_ra_colnames='RAJ2000',
    source_dec_colnames='DEJ2000', source_colors='cyan', source_edge_colors=None,
    source_markers=None, source_sizes=None,
    source_labels=None, label_size=15, label_offset=[0, 0],
    plotname='shell_stamp.png', return_fig=False, nan_color='white',
    stretch='linear', vmin=0, vmax=3000, plot_simbad_sources=True, simbad_type='star', simbad_color='cyan',
    dist=orion_dist , cbar_label=r'T$_{MB}v$ [K m/s]', show_colorbar=True,
    auto_scale=True, auto_scale_mode='min/max', auto_scale_nsigma=1.,
    auto_scale_pad_factor=None, dpi=300, cmap='viridis'):
    """
    Parameters
    ----------
    map : str, see below
        If str, denotes a fits image.
        Can also pass in any HDU  objects that can
        be used with aplpy.FITSFigure:
    
        astropy.io.fits.PrimaryHDU
        astropy.io.fits.ImageHDU pyfits.PrimaryHDU
        pyfits.ImageHDU 

    fig : None, optional
        Description
    shell : None, optional
        Shell object.
    ra : float, optional
        If a Shape object is not specified, use ra/dec/radius in degrees.
    dec : float, optional
        Declination of shell in degrees
    radius : float, optional
        Radius optionalf shell in degrees.
    circle_color : str, optional
        Description
    pad_factor : float, optional
        Size of stamp to plot, where size = radius * pad_factor
    contour_map : same as map, optional
        Optionally plot contours on top of the primary map.
    contour_levels : float, optional
        Description
    source_ra : list or ndarray, optional
        Source ra in degrees
    source_dec : list or ndarray, optional
        Source dec in degrees
    source_lists : table or list of tables
        Can be fits table, or HDU, which specifies the source catalogs to
        use for plotting sources.
    source_ra_colnames : str, optional
        Description
    source_dec_colnames : str, optional
        Description
    source_colors : str, optional
        Description
    plotname : str, optional
        Plot to write out.
    return_fig : bool, optional
        Description
    stretch : str, optional
        Description
    vmin : int, optional
        Description
    vmax : int, optional
        Description
    plot_simbad_sources : bool, optional
        Description
    simbad_type : str, optional
        Description
    simbad_color : str, optional
        Description
    dist : TYPE, optional
        Description
    cbar_label : str, optional
        Description
    auto_scale : bool, optional
        If True, run the auto_scale_mode to determine the pixel value scale.
    auto_scale_mode : str, optional
        If 'min/max', use vmin and vmax as absolute ceilin and floor display values,
        and use the local min/max inside the plot area if more restrictive than vmin/vmax.
    auto_scale_pad_factor : None, optional
        If float, use this to pad the area around the shell
        for the auto scaling instead of the full display area.
    
    Raises
    ------
    Exception
        Description
    
    Deleted Parameters
    ------------------
    shellShape : pyregion.Shape, optional
        Shape object denoting the shell's center and radius.
    return_subplot : bool, optional
        Return a subplot that can be used in multipanel plots.
    """
    #map_list = list(map)

    try:
        ra, dec, radius = shell.ra.value, shell.dec.value, shell.radius.value
    except:
        pass

        try:
            #If map is a str fits filename.
            hdu = fits.open(map)[0]
        except OSError:
            hdu = map

        wcs = WCS(hdu)

    # for map in map_list:
    #     try:
    #         #If map is a str fits filename.
    #         hdu = fits.open(map)[0]
    #     except OSError:
    #         hdu = map

    #     wcs = WCS(hdu)
    #     pix_xy = wcs.all_world2pix(ra, dec, 0)
    #     #Check if the center of the shell is contained within this map, if not, move onto
    #     #the next map in the list.
    #     if (pix_xy[0] >= 0) & (pix_xy[0] <= hdu.shape[1]) \
    #     & (pix_xy[1] >= 0) & (pix_xy[1] <= hdu.shape[0]):
    #         break


    try:
        #IF a FITSFigure is not passed into function as fig, but a map is specified.
        fig = aplpy.FITSFigure(map)
        fig.recenter(ra, dec, radius*pad_factor) #Pan/zoom to shell
    except:
        raise

    fig.recenter(ra, dec, radius*pad_factor) #Pan/zoom to shell 
    
    fig.tick_labels.set_yformat("dd:mm")
    fig.tick_labels.set_xformat("hh:mm:ss")

    #Auto Scaling

    if auto_scale:
        ra_grid, dec_grid = shells.worldgrid(hdu, returnorder='radec', returnunit='deg')
        #print(ra_grid, dec_grid)
        if auto_scale_pad_factor:
            mask = (abs(ra_grid.value - ra) < radius * auto_scale_pad_factor) &\
               (abs(dec_grid.value - dec) < radius * auto_scale_pad_factor)
        else:
            mask = (abs(ra_grid.value - ra) < radius * pad_factor) &\
               (abs(dec_grid.value - dec) < radius * pad_factor)
        mask_pixels = hdu.data[mask]
        #print(mask_pixels)

        if auto_scale_mode == 'min/max':
            if (vmin is not None) & (vmax is not None):
                vmin, vmax = np.nanmax([np.nanmin(mask_pixels), vmin]), np.nanmin([np.nanmax(mask_pixels), vmax])
            else:
                vmin, vmax = np.nanmin(mask_pixels), np.nanmax(mask_pixels)

        if auto_scale_mode == 'median':
            vmin, vmax = np.nanmin(mask_pixels), np.nanmedian(mask_pixels) + auto_scale_nsigma*np.nanstd(mask_pixels)

    #print(vmin, vmax)
    if rgb:
        fig.show_rgb(rgb)
    else:   
        fig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap=cmap, interpolation='none')


    #CONTOURS
    if contour_map:
        contour = fig.show_contour(contour_map, levels=contour_levels, colors=contour_color)
        #print("Contour levels: {}".format(contour_levels))
    #COLORBAR
    if show_colorbar:
        fig.add_colorbar()
        cb = fig.colorbar
        cb.set_axis_label_text(cbar_label)

    #SHELL OUTLINE
    print(type(ra), type(dec), type(radius))
    fig.show_circles(ra, dec, radius, linestyle='solid', edgecolor=circle_color,
        facecolor='none', linewidth=3)

    if thickness:
        fig.show_circles(ra, dec, radius-thickness/2., linestyle='dashed', edgecolor=circle_color,
        facecolor='none', linewidth=3)
        fig.show_circles(ra, dec, radius+thickness/2., linestyle='dashed', edgecolor=circle_color,
        facecolor='none', linewidth=3)




    #SCALEBAR
    fig.add_scalebar(206265 * scalebar_pc / (dist.to(u.pc).value * 3600), color='white',
        corner='bottom')
    fig.scalebar.set_label("{} pc".format(scalebar_pc))

    #POINT SOURCES
    if np.any(source_ra[0]) and np.any(source_ra[1]):
        print("Using source_ra/dec these parameters are set.")
        try:
            fig.show_markers(source_ra, source_dec,
             c=source_colors, marker=source_markers,
             edgecolor=source_edge_colors, s=source_sizes, label=source_labels)

            for j in range(len(source_ra)):
                fig.add_label(source_ra[j] + label_offset[0], source_dec[j] + label_offset[1], source_labels[j],
                 color=source_colors, size=label_size)

        except TypeError:
            #If more than one source list to be plotted with different markers
            #source_ra, source_dec, source_colors must be nested lists or ndarrays 
            #with same shape.
            #print(source_ra, source_dec, source_colors, source_markers, source_sizes)
            print("""
                #If more than one source list to be plotted with different markers
            #source_ra, source_dec, source_colors must be nested lists or ndarrays 
            #with same shape.
            #print(source_ra, source_dec, source_colors, source_markers, source_sizes)
                """)
            for i in range(len(source_colors)):
                # in_shell = (source_ra[i]**2. + source_dec[i]**2.) <= radius**2.
                # source_ra[i] = source_ra[i][in_shell]
                # source_dec[i] = source_dec[i][in_shell]
                print("Showing markers.")
                fig.show_markers(source_ra[i], source_dec[i], c=source_colors[i],
                 marker=source_markers[i], s=source_sizes[i],
                 edgecolor=source_edge_colors[i])   

                for j in range(len(source_ra[i])):
                    try:
                        print("Adding {} out of {} labels.".format(j+1, len(source_ra[i])))
                        fig.add_label(source_ra[i][j] + label_offset[0],
                         source_dec[i][j] + label_offset[1], source_labels[i][j],
                         color=source_colors[i], size=label_size)
                    except TypeError:
                        pass

    elif source_lists and type(source_lists) is not list:
        try:
            #source_lists is a single fits file string
            source_hdu = fits.open(source_lists)[-1] #Select last HDU in HDUList
        except OSError:
            #source_lists is a single hdu object
            source_hdu = source_lists

        source_ra = source_hdu.data[source_ra_colnames]
        source_dec = source_hdu.data[source_dec_colnames]
        fig.show_markers(source_ra, source_dec, edgecolor=source_colors)

    elif source_lists and type(source_lists) is list:
        raise Exception("source_lists of type(list) not implemented.")
        if type(source_lists[0]) is str:
        # source_lists is a list of fits file strings
            raise Exception("Souce_lists is a list of fits file strings not implemented.")
        else:
            # source_lists is a list of hdu objects
            raise Exception("Souce_lists is a list of hdu objects not implemented.")

    if plot_simbad_sources:
        print("plot_simbad_sources = True, so plotting Simbad sources.")
        simbad_table = simbad_sources(ra, dec, radius, unit='deg')
        simbad_coords = coord.SkyCoord(simbad_table['RA'], simbad_table['DEC'],
            unit=(u.hourangle, u.deg))

        fig.show_markers(simbad_coords.ra, simbad_coords.dec, edgecolor=simbad_color)

    fig.recenter(ra, dec, radius*pad_factor) #Pan/zoom to shell 
    if return_fig:
        return fig
    else:
        fig.save(plotname, dpi=dpi)
    fig.close()

if __name__ == '__main__':
    main()