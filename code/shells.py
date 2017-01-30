#! 
"""Summary
"""
from spectral_cube import SpectralCube
import pyregion
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.wcs import WCS
from wcsaxes import WCSAxes
from astropy.io import fits
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import aplpy
from aplpy import FITSFigure
import numpy as np

orion_dist = 414*u.pc #pc

def main():
    """Summary
    
    Returns
    -------
    TYPE
        Description
    """
    #cube_file = '../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits'
    region_file = '../nro_maps/SouthShells.reg'
    N = 2 # Number of shell candidates to plot
    shell_list = get_shells(region_file=region_file)

    cube_file = "../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits"
    #cube_file = "../nro_maps/13CO_20161011_FOREST-BEARS_xyb_spheroidal_dV0.11kms_YS.fits"

    for n in range(1,2):
        shell = shell_list[n]
        for deg in np.linspace(0, 180, 13):
            angle = deg*u.deg
            pv = plot_pv(cube=cube_file, ra_center=shell.ra, dec_center=shell.dec,
                vel=[4*u.km/u.s, 8*u.km/u.s], length=shell.radius*2.*4.,
                width=7.5*u.arcsec, angle=105*u.deg,
                pad_factor=1., plotname='12co_pv_shell'+str(n+1)+'_angle'+str(angle.value)+'morev.png', return_subplot=True,
                stretch='linear', auto_scale=True)

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

def movie(file='cube_57_90_noboundary.fits'):
    f = fits.open(file)
    header = f[0].header
    f.close()
    vel_min = header['CRVAL3'] - header['CDELT3'] * (header['CRPIX3'] - 1)
    N = header['NAXIS3']
    for i in range(N):
        fig = FITSFigure(file, dimensions=[0,1], slices=[i, 0])
        vel_kms = round((vel_min + header['CDELT3']*i) / 1000., 2)
        fig.show_grayscale()
        fig.set_tick_labels_xformat("hh:mm")
        fig.add_label(84.25,-5, str(vel_kms)+" km/s",size=15)
        fig.save('12co_carmanro_'+str(vel_kms)+'_kms.png', dpi=300) 


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




def get_shells(velocity_file='SouthShells_vrange.txt', region_file='../nro_maps/SouthShells.reg',
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

#    def subplot(source_catalog='../')

def plot_overview(cube='../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits',
 region_file='../nro_maps/SouthShells.reg', mode='peak', plotname='12co_peak_shells.png',
 interactive=False, show_shells=True, circle_color='blue', circle_linewidth=5, circle_style="dashed", return_fig=False, show=True):
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
        image = cube.max(axis=0)

    if mode == "mom0":
        image = cube.moment0().hdu

    fig = FITSFigure(image)
    if show:
        fig.show_grayscale()

    #plt.title(r"$^{12}$CO Peak T$_{MB}$")
    #plt.xlabel("RA (J2000)")
    #plt.ylabel("DEC (J2000)")

    if show_shells:
        shell_list = get_shells(region_file=region_file)
        for shell in shell_list:
            fig.show_circles(shell.ra.value, shell.dec.value, shell.radius.value, linestyle=circle_style, edgecolor=circle_color,
        facecolor='none', linewidth=circle_linewidth)


    if return_fig:
        return fig
    else:
        fig.save()
def plot_stamp(map=None, fig=None, shell=None, ra=None, dec=None, radius=None, circle_color='red',
    pad_factor=1., contour_map=None, source_ra=None, source_dec=None, source_lists=None, source_ra_colnames='RAJ2000',
    source_dec_colnames='DEJ2000', source_colors='cyan', plotname='shell_stamp.png', return_fig=False,
    stretch='linear', vmin=20, vmax=50, plot_simbad_sources=True, simbad_type='star', simbad_color='cyan',
    dist=orion_dist, cbar_label=r'T$_{MB}v$ [K m/s]'):
    """
    Parameters
    ----------
    map : str, see below
        If str, denotes a fits image.
        Can also pass in any HDU or WCS objects that can
        be used with aplpy.FITSFigure:
    
        astropy.io.fits.PrimaryHDU
        astropy.io.fits.ImageHDU pyfits.PrimaryHDU
        pyfits.ImageHDU astropy.wcs.WCS
        np.ndarray RGB image with AVM meta-data
    
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
    
    Raises
    ------
    Exception
        Description
    
    Deleted Parameters
    ------------------
    shellShape : pyregion.Shape, optional
        Shape object denoting the shell's center and radius.
    """
    if map and not fig:
        #IF a FITSFigure is not passed into function as fig, but a map is specified.
        fig = aplpy.FITSFigure(map)

    try:
        ra, dec, radius = shell.ra.value, shell.dec.value, shell.radius.value
    except:
        pass

    fig.recenter(ra, dec, radius*pad_factor) #Pan/zoom to shell 
    fig.tick_labels.set_yformat("dd:mm")
    fig.tick_labels.set_xformat("hh:mm:ss")

    fig.show_grayscale(stretch=stretch, vmin=vmin, vmax=vmax)

    #COLORBAR
    fig.add_colorbar()
    cb = fig.colorbar
    cb.set_axis_label_text(cbar_label)

    #SHELL OUTLINE
    fig.show_circles(ra, dec, radius, linestyle='dashed', edgecolor=circle_color,
        facecolor='none', linewidth=5)

    #SCALEBAR
    fig.add_scalebar(206265 * 0.2 / (dist.to(u.pc).value * 3600))
    fig.scalebar.set_label("0.2 pc")

    #POINT SOURCES
    if source_ra and source_dec:
        #Use source_ra/dec if these parameters are set.
        try:
            fig.show_markers(source_ra, source_dec, edgecolor=source_colors)
        except ValueError:
            #If more than one source list to be plotted with different markers
            #source_ra, source_dec, source_colors must be nested lists or ndarrays 
            #with same shape.
            for ra_list, dec_list, color in source_ra, source_dec, source_colors:
                fig.show_markers(ra_list, dec_list, edgecolor=color)

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
        simbad_table = simbad_sources(ra, dec, radius, unit='deg')
        simbad_coords = coord.SkyCoord(simbad_table['RA'], simbad_table['DEC'],
            unit=(u.hourangle, u.deg))

        fig.show_markers(simbad_coords.ra, simbad_coords.dec, edgecolor=simbad_color)

    if return_fig:
        return fig
    else:
        fig.save(plotname)
    fig.close()

def plot_channels(cube=None, shellShape=None, ra=None, dec=None, radius=None, circle_color='red',
    pad_factor=1., contour_map=None, source_ra=None, source_dec=None, source_lists=None, source_ra_colnames='RAJ2000',
    source_dec_colnames='DEJ2000', source_colors='blue', plotname='shell_stamp.png', return_subplot=True,
    stretch='linear', vmin=20, vmax=50, plot_simbad_sources=True, simbad_type='star', simbad_color='green',
    vel_min=0, vel_max=10, vel_unit=u.km/u.s, n_cols=3, chan_step=1, auto_scale=True):
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
    
    while n_chan > 16:
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
        subplot.tick_labels.set_xformat("hh:mm:ss")
        subplot.tick_labels.set_style('plain')
        if not left_edge:
            subplot.hide_ytick_labels()
            subplot.hide_yaxis_label()
        if not bottom_edge:
            subplot.hide_xtick_labels()
            subplot.hide_xaxis_label()


        #AUTO SCALING
        if auto_scale:
            shell_pixels = spec_cube[chan][shell_mask].value
            vmin, vmax = np.nanmin(shell_pixels), np.nanmax(shell_pixels)
        subplot.show_grayscale(stretch=stretch, vmin=vmin, vmax=vmax)     


        #COLORBAR
        #subplot.add_colorbar()
        #cb = subplot.colorbar
        #cb.set_axis_label_text(r'T$_{MB}$ [K]')

        #SHELL OUTLINE
        subplot.show_circles(ra, dec, radius, linestyle='dashed', edgecolor='green',
            facecolor='none')

        #POINT SOURCES
        if source_ra and source_dec:
            #Use source_ra/dec if these parameters are set.
            try:
                subplot.show_markers(source_ra, source_dec, edgecolor=source_colors)
            except ValueError:
                #If more than one source list to be plotted with different markers
                #source_ra, source_dec, source_colors must be nested lists or ndarrays 
                #with same shape.
                for ra_list, dec_list, color in source_ra, source_dec, source_colors:
                    subplot.show_markers(ra_list, dec_list, edgecolor=color)

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
            subplot.show_markers(simbad_coords.ra, simbad_coords.dec, edgecolor=simbad_color)

            ### LABEL THE VELOCITY OF EACH CHANNEL MAP
        subplot.add_label(0.7, 0.9,
            str(np.round(spec_cube.spectral_axis[chan].to('km/s'), 1)),
            color='red', relative=True)
        #subplot.close()
    #Make the figure prettier.
    #fig.tight_layout(h_pad=0, w_pad=0)

    frac = 0.08
    fig.subplots_adjust(bottom=frac, left=frac * n_rows / n_cols, 
        top=1 - frac, right=1 - (frac * n_rows / n_cols),
        wspace=0., hspace=0.)
    #fig.subplots_adjust(wspace=None, hspace=None)
    fig.canvas.draw()
    fig.savefig(plotname, dpi=300, bbox_inches='tight')
    #fig.close()

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
    stretch='linear', auto_scale=True):
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
    fig.show_grayscale(aspect='auto')
    fig.save(plotname)
    

def subcubes_from_ds9(cube, region_file='../nro_maps/SouthShells.reg', pad_factor=1., shape='exact'):
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
