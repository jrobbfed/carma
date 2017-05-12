#! 
"""Summary
"""
from spectral_cube import SpectralCube
import pyregion
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.wcs import WCS
from wcsaxes import WCSAxes
from astropy.io import (fits, ascii)
from astropy.table import Table, Column
#from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import aplpy
from aplpy import FITSFigure
import numpy as np
from astroquery.simbad import Simbad
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib

orion_dist = 414*u.pc #pc
shells_score3 = [3, 6, 9, 11, 17, 18, 21, 24, 25, 30, 37]

#fwhm_to_


def paper_figures(dir='~/carma/paper/figs/'):
    pass

def main():
    """Summary
    
    Returns
    -------
    TYPE
        Description
    """
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

    matplotlib.rc('font', **font)

    ###Plot overviews

#     plot_overview(cube='../combined_maps/12co_pix_2.cm.fits', plotname="12co_combined_peak_full.png",
#      show_shells=False,title=r"Combined $^{12}$CO Peak T$_{MB}$",
#      dist=orion_dist, vmin=None, vmax=None, scalebar_color='white',
#      scalebar_pc=1.,recenter=False, ra=83.99191, dec=-5.6611303, radius=0.117325)

#     plot_overview(plotname="12co_nroonly_peak_full.png", show_shells=False,
#      dist=orion_dist, vmin=None, vmax=None, scalebar_color='black', scale_factor =
#      1.,title=r"NRO $^{12}$CO Peak T$_{MB}$",
#      scalebar_pc=1,recenter=False, ra=83.99191, dec=-5.6611303, radius=0.117325)

    # plot_overview(cube='../combined_maps/12co_pix_2.cm.fits', plotname="12co_combined_mom0_cometary.png",
    #  show_shells=False, title=r"Combined Integrated $^{12}$CO",
    #  dist=orion_dist, scalebar_color='white', pmax=93., mode='mom0',
    #  scale_factor=1./1000,
    #  scalebar_pc=0.2,recenter=True, ra=83.99191, dec=-5.6611303, radius=0.117325)

    # plot_overview(plotname="12co_nroonly_mom0_cometary.png", show_shells=False,
    #  dist=orion_dist, scalebar_color='white', pmax=93., mode='mom0',
    #  scale_factor=1./1000, title=r"NRO Integrated $^{12}$CO",
    #  scalebar_pc=0.2,recenter=True, ra=83.99191, dec=-5.6611303, radius=0.117325)

    plot_overview(cube='../combined_maps/12co_pix_2.cm.fits', plotname="12co_combined_peak_full_shells.png",
     show_shells=True, shells_highlight=shells_score3, title=r"Combined $^{12}$CO Peak T$_{MB}$",
     dist=orion_dist, vmin=None, vmax=None, scalebar_color='white', circle_style='dotted',
     scalebar_pc=1.,recenter=False, ra=83.99191, dec=-5.6611303, radius=0.117325)

    return

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
    shell_list = get_shells(region_file=region_file, velocity_file=vrange_file)

    obaf_file = 'stars_obaf.txt'
    yso_file = "../catalogs/spitzer_orion.fit"

    obaf = ascii.read(obaf_file)
    obaf_ra, obaf_dec, obaf_label = np.array(obaf['RA']), np.array(obaf['DEC']), np.array([sp.strip("b'") for sp in obaf['SP_TYPE']])
    yso = fits.open(yso_file)[1].data
    yso_ra, yso_dec, yso_label = yso['RAJ2000'], yso['DEJ2000'], yso['Cl']

    # for nshell in range(19,43):
    #     shell = shell_list[nshell-1]
    #     ra, dec, radius = shell.ra.value, shell.dec.value, shell.radius.value

    #     #Check whether shell is in each mips image coverage.
    #     l1641_xy = WCS(mips_l1641_hdu).all_world2pix(ra, dec, 0)
    
    #     if  (l1641_xy[0] >= 0) & (l1641_xy[0] <= mips_l1641_hdu.shape[1]) & \
    #         (l1641_xy[1] >= 0) & (l1641_xy[1] <= mips_l1641_hdu.shape[0]):
    #         hdu = mips_l1641_hdu
    #     else:
    #         hdu = mips_onc_hdu

    #     plot_stamp(map=hdu, ra=ra, dec=dec, radius=radius, circle_color='red',
    #         pad_factor=1.5, contour_map=None, contour_levels=5., source_ra=None, source_dec=None, source_lists=None, 
    #         source_colors='cyan', plotname='{}shell{}_stamp.png'.format('MIPS',nshell), return_fig=False,
    #         stretch='linear', plot_simbad_sources=False, dist=orion_dist, cbar_label=r'counts',
    #         auto_scale=True, auto_scale_mode='min/max', auto_scale_pad_factor=1., vmin=0, vmax=3000)


    #cube_file = '../nro_maps/13CO_20161011_FOREST-BEARS_xyb_spheroidal_dV0.11kms_YS.fits'
    cube_file = '../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits'
    ir_l1641_hdu = fits.open(irac4_l1641_file)[0]
    ir_onc_hdu = fits.open(irac4_onc_file)[0]

    spec_cube = SpectralCube.read(cube_file)
    ra_grid = spec_cube.spatial_coordinate_map[1].to(u.deg).value
    dec_grid = spec_cube.spatial_coordinate_map[0].to(u.deg).value
    vel_grid = spec_cube.spectral_axis
    pad_factor = 1.5

    #plot_overview(show_shells=True)
    #plot_overview(plotname="12co_nro_peak.png", show_shells=False)
    #plot_overview(cube="/Volumes/Untitled/13co_pix_2.cm.fits", plotname="13co_combined_peak.png", show_shells=False)
    #return
    channel_vmax = [12.9, 14]
    for nshell in [30,]:
        shell = shell_list[nshell-1]
        ra, dec, radius = shell.ra.value, shell.dec.value, shell.radius.value

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
        sub_cube = spec_cube.with_mask(subcube_mask).minimal_subcube().spectral_slab(shell.vmin, shell.vmax)

        #Integrate between vmin and vmax.
        mom0_hdu = sub_cube.moment0().hdu

        # mask_inshell = (abs(ra_grid - ra) < radius) &\
        #        (abs(dec_grid - dec) < radius)
        # subcube_inshell = spec_cube.with_mask(mask_inshell).minimal_subcube().spectral_slab(shell.vmin, shell.vmax)
        # mom0_hdu_inshell = subcube_inshell.moment0().hdu

        #Calculate contour levels.
        empty_channel = spec_cube.closest_spectral_channel(500*u.Unit('m/s'))
        sigma = np.nanstd(spec_cube[empty_channel][subcube_mask]).value
        #print("sigma: {}".format(sigma))
        delta_vel = (sub_cube.spectral_extrema[1] - sub_cube.spectral_extrema[0]).value
        #print("delta_vel: {}".format(delta_vel))
        mom0_sigma = sigma * delta_vel
        #print(mom0_sigma)     
        #contour_levels = np.linspace(5.*mom0_sigma, np.nanmax(mom0_hdu_inshell.data), 12)
        contour_levels = np.linspace(28.*mom0_sigma, 45.*mom0_sigma, 6  )

        #Get source coordinates.


        # plot_stamp(map=ir_hdu, ra=ra, dec=dec, radius=radius, circle_color='red',
        #     pad_factor=pad_factor, contour_map=mom0_hdu, contour_levels=contour_levels, contour_color='white',
        #     plotname='{}shell{}_{}{}to{}_stamp.png'.format('8µm', nshell, "12CO", shell.vmin.value, shell.vmax.value),
        #     return_fig=False,
        #     stretch='linear', plot_simbad_sources=False, dist=orion_dist,
        #     auto_scale=True, auto_scale_mode='median', auto_scale_pad_factor=0.8, auto_scale_nsigma=4.,
        #     cbar_label="Counts", cmap='inferno',
        #     source_ra=[obaf_ra, yso_ra], source_dec=[obaf_dec, yso_dec],
        #     source_colors=['white', 'red'], source_markers=['*', 'None'], source_sizes=[300,50],
        #     source_labels=[obaf_label, yso_label], dpi=300
        #     )

        #cube_file = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
        
        
        plot_channels(cube=cube_file, ra=ra, dec=dec, radius=radius,
            source_lists=None, stretch='linear', pad_factor=1.5, vel_min=shell.vmin.value, vel_max=14.,
            plotname='12co_channels_shell'+str(nshell)+'.png', chan_step=2, plot_simbad_sources=False,
            vmin=None, vmax=None, max_chans=12,
            #cbar_label="Counts",
            source_ra=[obaf_ra, yso_ra], source_dec=[obaf_dec, yso_dec],
            source_colors=['white', 'red'], source_markers=['*', 'None'], source_sizes=[200,15], dpi=300)

        # angle = 90*u.deg
        # pv = plot_pv(cube=cube_file, ra_center=shell.ra, dec_center=shell.dec,
        #      vel=[shell.vmin - 1*u.km/u.s, shell.vmax + 1*u.km/u.s], length=shell.radius*4.,
        #      width=7.5*u.arcsec, angle=angle,
        #      pad_factor=1., plotname='12co_pv_shell'+str(nshell)+'_angle'+str(angle.value)+'.png',
        #      stretch='linear', auto_scale=True, dpi=900.)
    


    #simbad_brightstars(output='stars_obaf.txt', output_format='ascii', replace_ra='deg')


    #movie(test=False, labels=False) 

    #cube_file = '../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits'
    # region_file = '../nro_maps/SouthShells.reg'
    # N = 2 # Number of shell candidates to plot
    # shell_list = get_shells(region_file=region_file)

    # cube_file = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
    # #cube_file = "../nro_maps/13CO_20161011_FOREST-BEARS_xyb_spheroidal_dV0.11kms_YS.fits"

    # for n in range(1,2):
    #     shell = shell_list[n]
    #     for deg in np.linspace(0, 180, 13):
    #         angle = deg*u.deg
    #         pv = plot_pv(cube=cube_file, ra_center=shell.ra, dec_center=shell.dec,
    #             vel=[4*u.km/u.s, 8*u.km/u.s], length=shell.radius*2.*4.,
    #             width=7.5*u.arcsec, angle=105*u.deg,
    #             pad_factor=1., plotname='12co_pv_shell'+str(n+1)+'_angle'+str(angle.value)+'morev.png', return_subplot=True,
    #             stretch='linear', auto_scale=True)

    # cube_file = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
    # for n in range(N):
    #     shell = shell_list[n]
    #     plot_channels(cube=cube_file, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
    #         source_lists=None, stretch='linear', pad_factor=1.5, vel_min=shell.vmin.value, vel_max=shell.vmax.value,
    #         plotname='12co_channels_shell'+str(n+1)+'.png', chan_step=2, plot_simbad_sources=True, simbad_color='blue')

    # cube_file = "../nro_maps/13CO_20161011_FOREST-BEARS_xyb_spheroidal_dV0.11kms_YS.fits"
    # for n in range(N):
    #     shell = shell_list[n]
    #     plot_channels(cube=cube_file, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
    #         source_lists=None, stretch='linear', pad_factor=1.5, vel_min=shell.vmin.value, vel_max=shell.vmax.value,
    #         plotname='13co_channels_shell'+str(n+1)+'.png', chan_step=2, plot_simbad_sources=True, simbad_color='blue')
#--------------------------------------------------------------------       
    # cube_file = "../catalogs/IRAC_L1641_ch4_merged_clean.fits"
    # #source_file = '../catalogs/spitzer_orion.fit'
    # hdu = fits.open(cube_file)[0]
    # #std = np.nanstd(np.array(hdu.data)[100:200, 100:200])
    # wcs = WCS(hdu)
    # x_grid = np.indices(np.shape(hdu.data))[1]
    # y_grid = np.indices(np.shape(hdu.data))[0]

    # for n in range(N):
    #     shell = shell_list[n]
    #     #AUTO SCALE
            
    #     center = wcs.all_world2pix(shell.ra, shell.dec, 0) 
    #     radius = abs(wcs.all_world2pix(shell.ra + shell.radius, shell.dec, 0)[0] - center[0])
    #     shell_mask = (x_grid - center[0]) ** 2. + (y_grid - center[1]) ** 2. < radius ** 2.
    #     shell_data = hdu.data[shell_mask]
    #     vmin, vmax = np.nanmin(shell_data), np.nanmedian(shell_data) + 2*np.nanstd(shell_data)
    #     plot_stamp(map=hdu, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
    #         source_lists=None, stretch='linear', pad_factor=1.5, vmin=vmin, vmax=vmax,
    #         plotname='8micron_shell'+str(n+1)+'.png', cbar_label='')

    # cube_file = "../catalogs/MIPS_L1641a_24um.fits"
    # #source_file = '../catalogs/spitzer_orion.fit'
    # hdu = fits.open(cube_file)[0]
    # #std = np.nanstd(np.array(hdu.data)[100:200, 100:200])
    # wcs = WCS(hdu)
    # x_grid = np.indices(np.shape(hdu.data))[1]
    # y_grid = np.indices(np.shape(hdu.data))[0]

    # for n in range(N):
    #     shell = shell_list[n]
    #     #AUTO SCALE
            
    #     center = wcs.all_world2pix(shell.ra, shell.dec, 0) 
    #     radius = abs(wcs.all_world2pix(shell.ra + shell.radius, shell.dec, 0)[0] - center[0])
    #     shell_mask = (x_grid - center[0]) ** 2. + (y_grid - center[1]) ** 2. < radius ** 2.
    #     shell_data = hdu.data[shell_mask]
    #     vmin, vmax = np.nanmin(shell_data), np.nanmedian(shell_data) + 2*np.nanstd(shell_data)
    #     plot_stamp(map=hdu, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
    #         source_lists=None, stretch='linear', pad_factor=1.5, vmin=vmin, vmax=vmax,
    #         plotname='24micron_shell'+str(n+1)+'.png', cbar_label='')

#------------------------------------------------------------------       

#     cube_file = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
#     cube = SpectralCube.read(cube_file)
#     hdu = cube.moment0().hdu
#     #std = np.nanstd(np.array(hdu.data)[100:200, 100:200])
#     for n in range(N):
#         shell = shell_list[n]
#         #AUTO SCALE
#         ra_grid = cube.spatial_coordinate_map[1]
#         dec_grid = cube.spatial_coordinate_map[0]
#         shell_mask = (ra_grid - shell.ra) ** 2. + (dec_grid - shell.dec) ** 2. < shell.radius ** 2.
#         shell_data = hdu.data[shell_mask]
#         vmin, vmax = np.nanmin(shell_data), np.nanmax(shell_data)

#         plot_stamp(map=hdu, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
#             source_lists=None, stretch='linear', pad_factor=1.5, vmin=vmin, vmax=vmax,
#             plotname='12co_mom0_shell'+str(n+1)+'.png')
# #--------------------------------------------------------------------       
#     cube_file = "../nro_maps/13CO_20161011_FOREST-BEARS_xyb_spheroidal_dV0.11kms_YS.fits"
#     cube = SpectralCube.read(cube_file)
#     hdu = cube.moment0().hdu
#     #std = np.nanstd(np.array(hdu.data)[100:200, 100:200])
#     for n in range(N):
#         shell = shell_list[n]
#         ra_grid = cube.spatial_coordinate_map[1]
#         dec_grid = cube.spatial_coordinate_map[0]
#         shell_mask = (ra_grid - shell.ra) ** 2. + (dec_grid - shell.dec) ** 2. < shell.radius ** 2.
#         shell_data = hdu.data[shell_mask]
#         vmin, vmax = np.nanmin(shell_data), np.nanmax(shell_data)
#         plot_stamp(map=hdu, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
#             source_lists=None, stretch='linear', pad_factor=1.5, vmin=vmin, vmax=vmax,
#             plotname='13co_mom0_shell'+str(n+1)+'.png')


    # interactive=False
    # show_shells=True

    # plot_overview(cube=cube, mode=mode, region_file=region_file, plotname=plotname,
    #     interactive=interactive, show_shells=show_shells)
class Shell(object):
    """Summary
    
    Attributes
    ----------
    center : TYPE
        Description
    dec : TYPE
        Description
    ra : TYPE
        Description
    radius : TYPE
        Description
    vmax : TYPE
        Description
    vmin : TYPE
        Description
    """
    def __init__(self, ra=0., dec=0., radius=0., vmin=0., vmax=0.,
     ra_unit='degree', dec_unit='degree', radius_unit='degree', vunit='km/s'):
        """Summary
        
        Parameters
        ----------
        ra : float, optional
            Description
        dec : float, optional
            Description
        radius : float, optional
            Description
        vmin : float, optional
            Description
        vmax : float, optional
            Description
        ra_unit : str, optional
            Description
        dec_unit : str, optional
            Description
        radius_unit : str, optional
            Description
        vunit : str, optional
            Description
        """
        self.ra = coord.Angle(ra, ra_unit)
        self.dec = coord.Angle(dec, dec_unit)
        self.center = coord.SkyCoord(ra, dec, unit=(ra_unit, dec_unit))
        self.radius = coord.Angle(radius, radius_unit)
    
        self.vmin = vmin * u.Unit(vunit)
        self.vmax = vmax * u.Unit(vunit)

def worldgrid(hdu, wcs=None, origin=0, returnorder='radec', returnunit=None):
    """
    Returns an array of world coordinates corresponding to pixels in an image.
    Adapted for the 2D case from spectral_cube.SpectralCube.world()
    
    Parameters
    ----------
    hdu : astropy.io.fits.PrimaryHDU
        astropy.io.fits.ImageHDU pyfits.PrimaryHDU
        pyfits.ImageHDU or numpy.ndarray
    
        If an HDU, then the wcs will be found, or data ndarray
        and corresponding wcs object can be passed explicitly.
    
    wcs : None or astropy.wcs.WCS, optional
        If hdu is a ndarray, then wcs should be the corresponding WCS object.WCS 
    origin : int, optional
        0 for numpy-type indexing, 1 for fits-type indexing
    returnorder : str, optional
        'radec' or 'decra' for the order in which the coordinates are returned
    returnunit : str, optional
        The name of the astropy.unit to use for returning the world coordinates.
        If not specified, then use the units in wcs.
    
    """
    try:
        wcs = WCS(hdu)
    except AttributeError:
        #It's ok, if hdu is a data ndarray.
        pass

    inds = np.ogrid[[slice(0,s) for s in hdu.shape]]
    inds = np.broadcast_arrays(*inds)
    inds = inds[::-1] # numpy -> wcs axis order (RA/DEC)
    shp = inds[0].shape
    inds = np.column_stack([i.ravel() for i in inds])
    print(inds)
    try:
        world = wcs.wcs_pix2world(inds, 0).T
        print(world)
    except AttributeError:
        raise("Needs a valid wcs.")

    world = [w.reshape(shp) for w in world] #1D -> ND

    if returnunit:
        world = [w * u.Unit(returnunit) 
                for i, w in enumerate(world)]
    else:
        world = [w * u.Unit(wcs.wcs.cunit[i]) 
                for i, w in enumerate(world)]

    if returnorder == 'radec':
        return world
    elif returnorder == 'decra':
        return world[::-1]


def movie(file='../combined_maps/12co_01_90_noboundary.fits', test=False, labels=True):
    f = fits.open(file)
    header = f[0].header
    f.close()
    vel_min = header['CRVAL3'] - header['CDELT3'] * (header['CRPIX3'] - 1)
    if test:
        Ni = 59
        Nf = 61
    else:
        Ni = 0
        Nf = header['NAXIS3']
    for i in range(Ni, Nf):
        fig = FITSFigure(file, dimensions=[0,1], slices=[i, 0])
        vel_kms = round((vel_min + header['CDELT3']*i) / 1000., 2)
        fig.show_colorscale(cmap='viridis')
        if labels:
            fig.set_tick_labels_xformat("hh:mm")
            fig.add_label(84.25,-5, str(vel_kms)+" km/s",size=15)
        else:
            #fig.set_theme('pretty')
            fig.set_nan_color('black')
            #fig.hide_grid()
            fig.hide_axis_labels()
            fig.hide_tick_labels()
            fig.ticks.hide()
        fig.save('12co_movie/12co_carmanro_chan0'+str(i+1)+'_'+str(vel_kms)+'_kms.png', dpi=600) 
        fig.close()

def shell_mask(cube, shell):
    """
    Returns a boolean mask for the pixels inside shell.
    
    Parameters
    ----------
    cube : TYPE
        Description
    shell : TYPE
        Description
    """
    try:
        cube = SpectralCube.read(cube)
    except ValueError:
        pass
    wcs = cube.wcs
    ra_grid = cube.spatial_coordinate_map[1]
    dec_grid = cube.spatial_coordinate_map[0]

    center = wcs.all_world2pix(shell.ra, shell.dec, 0)
    radius = abs(wcs.all_world2pix(shell.ra + shell.radius, shell.dec, 0)[0] - center[0])
    mask = (x_grid - center[0]) ** 2. + (y_grid - center[1]) ** 2. < radius ** 2. 
    return mask

def simbad_brightstars(image_file="../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits",
    brighter_than='G0', extra_criteria="(ra < 84.4 | dec < -6.66)", otypes="Star",
    replace_ra='hourangle', replace_dec='deg', add_sptype_letter_column=True,
    output=None, output_format='fits'):
    """
    Returns a table of SIMBAD stars brighter than the given spectral type found within the footprint of
    the input cube or image.
    
    Parameters
    ----------
    image_file : str, optional
        Input fits cube or 2-d image with world coordinates in header.
    brighter_than : str, optional
        Stars selected will all have spectral types earlier than this.
    extra_criteria : str, optional
        In the form of "criteriaA & criteriaB & criteriaC". 
    otypes : str, optional
        The SIMBAD object type identifier to use.
        "Star" finds all objects of the Star category
             and all subcategories of Star like YSOs.
    replace_ra : str, optional
        If not None, replace the right ascension in table with values in this unit.
        Must be an astropy.unit name.
    replace_dec : str, optional
        If not None, replace the declination in table with values in this unit.
        Must be an astropy.unit name.
    add_sptype_letter_column : bool, optional
        If True, add a column to the table consisting of a single letter indicating the spectral type.
    output : str, optional
        If None, return table, else `output` is name of the file to write the star table to.
    
    """
    try:
        wcs = WCS(image_file).celestial #Drop non-celestial axes (like velocity and stokes). 
    except:
        raise("image_file must be a fits image or cube with wcs in header.")

    footprint = wcs.calc_footprint()

 
    ### ra_min/max, dec_min/max need to be in degrees.
    ### In the fits headers I have they are, but this may not always be true.
    ###
    ra_min, ra_max = footprint[:,0].min(), footprint[:,0].max()
    dec_min, dec_max = footprint[:,1].min(), footprint[:,1].max()

    s = Simbad()
    s.add_votable_fields('sptype')

    if extra_criteria:
        stars = s.query_criteria("ra > {} & ra < {} & dec > {} & dec < {} & sptypes < {} & {}".format(
            ra_min, ra_max, dec_min, dec_max, brighter_than, extra_criteria), otypes="Star")
    else:
        stars = s.query_criteria("ra > {} & ra < {} & dec > {} & dec < {} & sptypes < {}".format(
            ra_min, ra_max, dec_min, dec_max, brighter_than), otypes="Star")

    stars_coord = coord.SkyCoord(stars['RA'], stars['DEC'], unit=(u.hourangle, u.deg))

    if replace_ra:
        stars.replace_column('RA', Column(stars_coord.ra, name='RA', unit=replace_ra))
    if replace_dec:
        stars.replace_column('DEC', Column(stars_coord.dec, name='DEC', unit=replace_dec))

    if add_sptype_letter_column:
        stars.add_column(Column([sptype[0] for sptype in stars['SP_TYPE'].astype('str')], name='SP_LETTER', unit='str'))

    if output:
        stars.write(output, format=output_format)##
    else:
        return stars



def simbad_sources(ra, dec, radius, unit='degree', object_type="Star", fields='sptype'):
    """Returns the simbad sources within `radius` around [`ra`, `dec`].
    
    Parameters
    ----------
    ra : float
        Right ascencion in units of `unit`
    dec : float
        Declination in units of `unit`
    radius : float
        Radius in units of `unit`
    unit : str, optional
        Astropy angular unit to use for ra, dec, and radius.
    object_type : str, optional
        SIMBAD identifier type.
    
    Returns
    -------
    astropy Table
        Table of the sources within the requested region of the requested object type.
    
    Deleted Parameters
    ------------------
    type : str, optional
        SIMBAD identifier type.
    """
    from astroquery.simbad import Simbad
    c = coord.SkyCoord(ra, dec, unit=unit)
    r = radius * u.Unit(unit)
    s = Simbad()
    s.add_votable_fields(fields)
    result_table = s.query_criteria("region(circle, " + str(c.ra.to(u.deg).value) + " "
        + '+'*(c.dec.value >= 0) + str(c.dec.to(u.deg).value) + ', ' + str(r.to(u.deg).value) + 'd)', 
        otypes=object_type)
    return result_table
    #result_table = Simbad.query_region(c, radius=r) 
    
def data_cutout(data=None, wcs=None, ra=None, dec=None, radius=None, unit='deg', shape='circle'):
    """
    Takes a WCS object and the dimensions of a shape in degrees and returns the pixel based mask
    for the WCS-based array. 
    
    Parameters
    ----------
    data : None, optional
        Description
    wcs : None, optional
        Description
    ra : None, optional
        Description
    dec : None, optional
        Description
    radius : None, optional
        Description
    unit : str, optional
        Description
    shape : str, optional
        Description
    """
    from astropy.coordinates import Angle, SkyCoord
    if shape == 'circle':
        from regions.shapes.circle import CirclePixelRegion
        from regions import CircleSkyRegion
        sky_reg = CircleSkyRegion(SkyCoord(ra, dec, unit=unit), Angle(radius, unit))

    pix_reg = sky_reg.to_pix(wcs)

    mask.cutout(da)
def get_shells(velocity_file='../shell_candidates/AllShells_vrange.txt',
<<<<<<< HEAD
    region_file='../shell_candidates/AllShells.reg',
=======
    region_file='../shell_candidates/SouthShells.reg',
>>>>>>> e604d1952303030a5d070108bfd382af1ffe382c
    ra_col="ra", dec_col="dec", radius_col="radius", vmin_col='vmin', vmax_col='vmax',
    ra_unit='deg', dec_unit='deg', radius_unit='deg', v_unit='km/s'):
    """
    Read a ds9 region file and return a list of Shell objects corresponding to the
    regions in the file. 
    
    Parameters
    ----------
    velocity_file : str, optional
        Description
    region_file : str, optional
        Description
    ra_col : str, optional
        Description
    dec_col : str, optional
        Description
    radius_col : str, optional
        Description
    vmin_col : str, optional
        Description
    vmax_col : str, optional
        Description
    ra_unit : str, optional
        Description
    dec_unit : str, optional
        Description
    radius_unit : str, optional
        Description
    v_unit : str, optional
        Description
    """
    shell_list = []
    try:
        region_list = pyregion.open(region_file)
        vel_table = ascii.read(velocity_file)
        vmin_list, vmax_list = vel_table[vmin_col], vel_table[vmax_col]
    except ValueError:
        raise

    for i, region in enumerate(region_list):
        ra, dec, radius = region.coord_list[0], region.coord_list[1], region.coord_list[2]
        shell_list += [Shell(ra, dec, radius, vmin_list[i], vmax_list[i])]


    return shell_list

def plot_overview(cube='../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits',
 region_file='../nro_maps/AllShells.reg', mode='peak', plotname='12co_peak_shells.png',
 interactive=False, show_shells=True, shells_highlight=None, dist=orion_dist, vmin=None, vmax=None,
 scalebar_color="white", scalebar_pc = 1., scale_factor=1., pmin=0.25,
 pmax=99.75, cbar_label=r"T$_{MB}$ v [K km/s]",
 circle_color='white', circle_linewidth=1, circle_style="solid", return_fig=False, show=True,
 title=r"$^{12}$CO Peak T$_{MB}$", recenter=False, ra=None, dec=None, radius=None):
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

    #RECENTER
    if recenter:
    	fig.recenter(ra, dec, radius)


    #SCALEBAR
    fig.add_scalebar(206265 * scalebar_pc / (dist.to(u.pc).value * 3600), color=scalebar_color)
    fig.scalebar.set_label("{} pc".format(scalebar_pc))

    fig.add_colorbar()
    cb = fig.colorbar
    cb.set_axis_label_text(cbar_label)

    if return_fig:
        return fig
    else:
        fig.save(plotname, dpi=600)
def plot_stamp(map=None, fig=None, shell=None, ra=None, dec=None, radius=None, circle_color='red',
    pad_factor=1.5, contour_map=None, contour_levels=5, contour_color='white',
    source_ra=None, source_dec=None, source_lists=None, source_ra_colnames='RAJ2000',
    source_dec_colnames='DEJ2000', source_colors='cyan', source_markers=None, source_sizes=None,
    source_labels=None,
    plotname='shell_stamp.png', return_fig=False,
    stretch='linear', vmin=0, vmax=3000, plot_simbad_sources=True, simbad_type='star', simbad_color='cyan',
    dist=orion_dist, cbar_label=r'T$_{MB}v$ [K m/s]',
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
        ra_grid, dec_grid = worldgrid(hdu, returnorder='radec', returnunit='deg')
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
    fig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap=cmap, interpolation='none')

    #CONTOURS
    if contour_map:
        contour = fig.show_contour(contour_map, levels=contour_levels, colors=contour_color)
        #print("Contour levels: {}".format(contour_levels))
    #COLORBAR
    fig.add_colorbar()
    cb = fig.colorbar
    cb.set_axis_label_text(cbar_label)

    #SHELL OUTLINE
    fig.show_circles(ra, dec, radius, linestyle='dashed', edgecolor=circle_color,
        facecolor='none', linewidth=5)

    #SCALEBAR
    fig.add_scalebar(206265 * 0.2 / (dist.to(u.pc).value * 3600), color='white')
    fig.scalebar.set_label("0.2 pc")

    #POINT SOURCES
    if source_ra and source_dec:
        #Use source_ra/dec if these parameters are set.
        try:
            fig.show_markers(source_ra, source_dec,
             c=source_colors, marker=source_markers, s=source_sizes)
        except TypeError:
            #If more than one source list to be plotted with different markers
            #source_ra, source_dec, source_colors must be nested lists or ndarrays 
            #with same shape.
            #print(source_ra, source_dec, source_colors, source_markers, source_sizes)
            for i in range(len(source_colors)):
                fig.show_markers(source_ra[i], source_dec[i], c=source_colors[i],
                 marker=source_markers[i], s=source_sizes[i])


    elif source_lists and type(source_lists) is not list:
        try:
            #source_lists is a single fits file string
            source_hdu = fits.open(source_lists)[-1] #Select last HDU in HDUList
        except OSError:
            #source_lists is a single hdu object
            source_hdu = source_lists

        source_ra = source_hdu.data[source_ra_colnames]
        source_dec = source_hdu.data[source_dec_colnames]
        fig.show_markers(source_ra, source_dec, edgecolor=source_colors, label=source_labels)

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

        fig.show_markers(simbad_coords.ra, simbad_coords.dec, edgecolor=simbad_color)

    if return_fig:
        return fig
    else:
        fig.save(plotname, dpi=dpi)
    fig.close()

def plot_channels(cube=None, shellShape=None, ra=None, dec=None, radius=None, circle_color='white',
    pad_factor=1., contour_map=None, source_ra=None, source_dec=None, source_lists=None, source_sizes=None,
    source_markers=None,
    source_ra_colnames='RAJ2000',
    source_dec_colnames='DEJ2000', source_colors='blue', plotname='shell_stamp.png', return_subplot=True,
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
        #subplot.ticks.set_xspacing(0.0001)
        subplot.tick_labels.set_xformat("hh:mm")
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
        subplot.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='viridis')     

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
                 c=source_colors, marker=source_markers, s=source_sizes)
            except TypeError:
                #If more than one source list to be plotted with different markers
                #source_ra, source_dec, source_colors must be nested lists or ndarrays 
                #with same shape.
                #print(source_ra, source_dec, source_colors, source_markers, source_sizes)
                for i in range(len(source_colors)):
                    subplot.show_markers(source_ra[i], source_dec[i], c=source_colors[i],
                     marker=source_markers[i], s=source_sizes[i])

        elif source_lists and type(source_lists) is not list:
            try:
                #source_lists is a single fits file string
                source_hdu = fits.open(source_lists)[-1] #Select last HDU in HDUList
            except OSError:
                #source_lists is a single hdu object
                source_hdu = source_lists

            source_ra = source_hdu.data[source_ra_colnames]
            source_dec = source_hdu.data[source_dec_colnames]
            subplot.show_markers(source_ra, source_dec, edgecolor=source_colors)

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
            color='white', relative=True)
        #subplot.close()
    #Make the figure prettier.
    #fig.tight_layout(h_pad=0, w_pad=0)

    frac = 0.08
    fig.subplots_adjust(bottom=frac, left=frac * n_rows / n_cols, 
        top=1 - frac, right=1 - (frac * n_rows / n_cols),
        wspace=0., hspace=0.)
    #fig.subplots_adjust(wspace=None, hspace=None)
    fig.canvas.draw()
    fig.savefig(plotname, dpi=dpi, bbox_inches='tight')
    #fig.close()

def channel_slicer(cube=None, ra=None, dec=None, radius=None,
                   title=None, pad_factor=2., chan_init=0,
                   stretch='linear', circle_color='blue',
                   circle_style='-'):
    spec_cube = SpectralCube.read(cube)
    vel_list = spec_cube.spectral_axis

    #For the auto color scaling to the min and max intensities in the zoomed region.
    ra_grid = spec_cube.spatial_coordinate_map[1].to(u.deg).value
    dec_grid = spec_cube.spatial_coordinate_map[0].to(u.deg).value
    #shell_mask = (ra_grid - ra) ** 2. + (dec_grid - dec) ** 2. < (radius) ** 2.
    subcube_mask = (abs(ra_grid - ra) < radius * pad_factor) &\
               (abs(dec_grid - dec) < radius * pad_factor)
    sub_cube = spec_cube.with_mask(subcube_mask).minimal_subcube()

    #auto color scaling
    subcube_pixels = sub_cube[chan_init].value
    vmin, vmax = np.nanmin(subcube_pixels), np.nanmax(subcube_pixels)
    print(vmin, vmax)


    #center plot on the shell
    fig = plt.figure()
    subplot = FITSFigure(sub_cube.hdu, figure=fig, slices=[chan_init], auto_refresh=True)
    subplot.set_title("{} @ {}".format(title, vel_list[chan_init]))
    #subplot.recenter(ra, dec, radius*pad_factor)

    #Make aplpy grayscale plot with nice color scaling/stretch
    subplot.show_grayscale(stretch=stretch, vmin=vmin, vmax=vmax)     

    subplot.tick_labels.set_yformat("dd:mm")
    subplot.tick_labels.set_xformat("hh:mm")
    subplot.tick_labels.set_style('plain')

#Set up slider to change spectral channels.
    axcolor = 'lightgoldenrodyellow'
    ax = fig.add_axes([0.25, 0.95, 0.65, 0.03], axisbg=axcolor)
    slider = Slider(ax, 'Channel', 0, vel_list.shape[0] - 1,
                    valinit=chan_init, valfmt='%i')

    #Show the shell as circle.
    subplot.show_circles(ra, dec, radius, edgecolor=circle_color)

    def update(val):
        ind = int(slider.val)
        
        #subplot = FITSFigure(hdu, figure=fig, slices=[chan], auto_refresh=True)
        #subplot.recenter(ra, dec, radius*pad_factor)
        #subplot.
        subplot.set_title("{} @ {}".format(title, vel_list[ind])) 
        subcube_pixels = sub_cube[ind].value
        vmin, vmax = np.nanmin(subcube_pixels), np.nanmax(subcube_pixels)
        subplot.image.set_clim(vmin, vmax)
        subplot.image.set_data(sub_cube[ind].hdu.data)
        
        #fig.canvas.draw()
        
        #subplot.show_grayscale(stretch=stretch, vmin=vmin, vmax=vmax)
        #subplot.set_data()
        #fig.canvas.draw()

#def write_velmin():
#    current_chan = int(slider.val)
#def write_velmax():
#    current_chan = inte(slider.val)
    
        
    slider.on_changed(update)

#button_velmin.on_click(write_velmin)
#button_velmax.on_click(write_velmax)

    plt.show()
def pv_slice(cube=None, ra_center=None, dec_center=None,
    width=1*u.arcsec, angle=0*u.deg, length=1*u.arcmin, return_path=False):
    """Summary
    
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
    angle : TYPE, optional
        Description
    length : TYPE, optional
        Description
    return_path : bool, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    from pvextractor import PathFromCenter, extract_pv_slice
    center = coord.SkyCoord(ra_center, dec_center)
    # print(center)
    path = PathFromCenter(center=center, length=length, angle=angle, width=width)
    pv_slice = extract_pv_slice(cube, path)
    if return_path:
        return (pv_slice, path)
    else:
        return pv_slice


def plot_pv(cube=None, ra_center=None, dec_center=None, vel=[None, None],
    #ra_unit='deg', dec_unit='deg',  vel_unit=u.km/u.s,
    width=1*u.arcsec, angle=0*u.deg, length=1*u.arcmin,
    pad_factor=1., plotname='shell_pv.png', return_subplot=True,
    stretch='linear', auto_scale=True, dpi=300.):
    """<FRESHLY_INSERTED>"""
    try:
       # if cube is not already a SpectralCube, but is a fits file.
        spec_cube = SpectralCube.read(cube)
    except ValueError:
        # `spec_cube` is a SpectralCube object
        cube = cube.hdu
    #Find the min and max spectral channel.
    
    chan_min = spec_cube.closest_spectral_channel(vel[0])
    chan_max = spec_cube.closest_spectral_channel(vel[1])

    pv = pv_slice(cube, ra_center, dec_center, width, angle, length)


    fig = FITSFigure(pv)
    #fig.tick_labels.set_yformat()
    
    if vel[0] is not None:
        #Only show the velocities of the shell.
        p_lim, v_lim = fig.pixel2world(list(fig._ax1.get_xlim()), list(fig._ax1.get_ylim()))
        p_unit = fig._header['CUNIT1']
        v_unit = fig._header['CUNIT2']

        #These are floats with the corresponding units of the figure axes.
        v_center = ((vel[1] + vel[0]) / 2.).to(u.Unit(v_unit)).value
        v_range = (vel[1] - vel[0]).to(u.Unit(v_unit)).value
        #Position axis is kept the same.
        p_center = (p_lim[1] + p_lim[0]) / 2.
        p_range = p_lim[1] - p_lim[0]
        fig.recenter(p_center, v_center, width=p_range, height=v_range)

    #fig.tick_labels.set_xformat('%4.2f')
    fig.show_colorscale(aspect='auto', cmap='viridis')
    fig.add_colorbar()
    cb = fig.colorbar
    cb.set_axis_label_text(r'T$_{MB}$ [K]')

    fig.save(plotname, dpi=dpi)
    

def subcubes_from_ds9(cube, region_file='../nro_maps/AllShells.reg', pad_factor=1., shape='exact'):
    """
    Extracts subcubes using the ds9 region file.
    
    Parameters
    ----------
    cube : SpectralCube, str
        The cube to be chopped. Must be type spectral_cube.SpectralCube or str filename.
    region_file : str
        Path to a ds9 region file.
    pad_factor : float, optional
        Expand the subcube around the region by this factor.
    shape : {'square', 'exact'}
        The shape of the subcube returned. 'square' returns the
        smallest square subcube that contains the region.
        'exact' returns only the pixels contained within the region.
    
    Returns
    -------
    subcubes: list of SpectralCube of SpectralCube
    """
    from spectral_cube import SpectralCube
    import pyregion

    try:
        #If cube is a str filename, read a SpectralCube.
        cube = SpectralCube.read(cube)
    except ValueError:
        pass

    if shape == 'square':
        import astropy.units as u
        subcube_list = []
        region_list = pyregion.open(region_file)
        for region in region_list:
            half_width = region.coord_list[2] * pad_factor * u.deg
            ra_center = region.coord_list[0] * u.deg
            dec_center = region.coord_list[1] * u.deg
            ra_range = [ra_center - half_width, ra_center + half_width]
            dec_range = [dec_center - half_width, dec_center + half_width]
            #print(ra_range, dec_range)
            subcube_list.append(cube.subcube(ra_range[1], ra_range[0], dec_range[0], dec_range[1]))
    if shape == 'exact':
        region_list = pyregion.open(region_file)
        subcube_list = []
        for region in region_list:
            
            if pad_factor != 1.:
                new_string = '{};{}({},{},{}")'.format(region.coord_format, region.name,
                                region.coord_list[0], region.coord_list[1],
                                region.coord_list[2]*3600.*pad_factor)
                region = pyregion.parse(new_string)[0]
                
            subcube_list.append(cube.subcube_from_ds9region(pyregion.ShapeList([region])))
    if len(subcube_list) == 1:
        return subcube_list[0]
    else:
        return subcube_list

if __name__ == "__main__":
    main()
