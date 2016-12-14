"""Summary
"""
from spectral_cube import SpectralCube
import pyregion
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.wcs import WCS
from wcsaxes import WCSAxes
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
import aplpy
import numpy as np

def main():
    """Summary
    
    Returns
    -------
    TYPE
        Description
    """
    #cube_file = '../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits'
    cube_file = "../catalogs/MIPS_L1641a_24um.fits"
    source_file = '../catalogs/spitzer_orion.fit'
    #cube = SpectralCube.read(cube_file)

    #hdu = cube.hdu
    hdu = fits.open(cube_file)[-1]
    #map_hdu = cube.max(axis=0).hdu
    #mom0 = cube.moment0()
    std = np.nanstd(np.array(hdu.data)[100:200, 100:200])

    #mom0hdu = mom0.hdu

    #mom0_fig = aplpy.FITSFigure(mom012co_hdu)
    #mode = 'peak'
    region_file = '../nro_maps/SouthShells.reg'
    #plotname='12co_peak.png'
    n = 10
    shell_list = get_shells(region_file=region_file)
    shell = shell_list[n]

    #std = mom012
    #plot_stamp(fig=mom012co_fig, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
    #    source_lists=source_file, stretch='linear', pad_factor=1.5, vmin=2*std, vmax=20*std)


    for n in range(12):
        shell = shell_list[n]
        plot_stamp(map=hdu, ra=shell.ra.value, dec=shell.dec.value, radius=shell.radius.value,
            source_lists=source_file, stretch='linear', pad_factor=1.5, vmin=0, vmax=800,
            plotname='mips_shell'+str(n+1)+'.png')


    # interactive=False
    # show_shells=True

    # plot_overview(cube=cube, mode=mode, region_file=region_file, plotname=plotname,
    #     interactive=interactive, show_shells=show_shells)

def get_shells(region_file='../nro_maps/SouthShells.reg'):
    """
    Read a ds9 region file and return a list of Shell objects corresponding to the
    regions in the file. 
    
    Parameters
    ----------
    region_file : str, optional
        Description
    """
    shell_list = []
    region_list = pyregion.open(region_file)

    for region in region_list:
        ra, dec, radius = region.coord_list[0], region.coord_list[1], region.coord_list[2]
        shell_list += [Shell(ra, dec, radius)]

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
     ra_unit='degree', dec_unit='degree', radius_unit='degree', vunit='kms'):
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
        if vunit == 'kms':
            self.vmin = vmin * u.km / u.s
            self.vmax = vmax * u.km / u.s

#    def subplot(source_catalog='../')

def plot_overview(cube='../nro_maps/12CO_20161002_FOREST-BEARS_spheroidal_xyb_grid7.5_0.099kms.fits',
 region_file='../nro_maps/SouthShells.reg', mode='peak', plotname='12co_peak_shells.png',
 interactive=False, show_shells=False):
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

    fig = plt.figure()
    wcs = WCS(image.header)
    ax = WCSAxes(fig, [0.1,0.1,0.8,0.8], wcs=wcs) 
    fig.add_axes(ax)      
    imgplot = plt.imshow(image.data, cmap=cm.gray, origin='lower', interpolation='none',
        vmin=0., vmax=100)
    cb = plt.colorbar()
    cb.set_label(r'K [T$_{MB}$]')
    plt.title(r"$^{12}$CO Peak")

    if show_shells:
        r = pyregion.open(region_file).as_imagecoord(image.header)
        patch_list, artist_list = r.get_mpl_patches_texts()

        for p in patch_list:
            ax.add_patch(p)
        for t in artist_list:
            ax.add_artist(t)

        pass

    if interactive:
        plt.show()
    else:
        plt.savefig(plotname)

def plot_stamp(map=None, fig=None, shellShape=None, ra=None, dec=None, radius=None, circle_color='red',
    pad_factor=1., contour_map=None, source_ra=None, source_dec=None, source_lists=None, source_ra_colnames='RAJ2000',
    source_dec_colnames='DEJ2000', source_colors='blue', plotname='shell_stamp.png', return_subplot=True,
    stretch='linear', vmin=20, vmax=50):
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
    
    shellShape : pyregion.Shape, optional
        Shape object denoting the shell's center and radius.
    ra : float, optional
        If a Shape object is not specified, use ra/dec/radius in degrees.
    dec : float, optional
        Declination of shell in degrees
    radius : float, optional
        Radius optionalf shell in degrees.
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
    
    Raises
    ------
    Exception
        Description
    """
    if map and not fig:
        #IF a FITSFigure is not passed into function as fig, but a map is specified.
        fig = aplpy.FITSFigure(map)

    fig.recenter(ra, dec, radius*pad_factor) #Pan/zoom to shell 
    fig.tick_labels.set_yformat("dd:mm")
    fig.tick_labels.set_xformat("hh:mm:ss")
    fig.show_grayscale(stretch=stretch, vmin=vmin, vmax=vmax)

    #COLORBAR
    fig.add_colorbar()
    cb = fig.colorbar
    cb.set_axis_label_text(r'T$_{MB}$ [K]')

    #SHELL OUTLINE
    fig.show_circles(ra, dec, radius, linestyle='dashed', edgecolor=circle_color,
        facecolor='none')

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

    elif type(source_lists) is not list:
        try:
            #source_lists is a single fits file string
            source_hdu = fits.open(source_lists)[-1] #Select last HDU in HDUList
        except OSError:
            #source_lists is a single hdu object
            source_hdu = source_lists
        source_ra = source_hdu.data[source_ra_colnames]
        source_dec = source_hdu.data[source_dec_colnames]
        fig.show_markers(source_ra, source_dec, edgecolor=source_colors)

    elif type(source_lists) is list:
        raise Exception("source_lists of type(list) not implemented.")
        if type(source_lists[0]) is str:
        # source_lists is a list of fits file strings
            raise Exception("Souce_lists is a list of fits file strings not implemented.")
        else:
            # source_lists is a list of hdu objects
            raise Exception("Souce_lists is a list of hdu objects not implemented.")


    fig.save(plotname)


def plot_channels():
    """Summary
    
    Returns
    -------
    TYPE
        Description
    """
    pass

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
