"""Summary
"""
from shells import pv_slice

def analytic_profile(r, R, dR, radius_type="middle"):
    """Summary
    
    Parameters
    ----------
    r : TYPE
        Description
    R : TYPE
        Description
    dR : TYPE
        Description
    radius_type : str, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    
    Raises
    ------
    Exception
        Description
    """
    #From Beaumont & Williams (2010), assumes bare optically thin shell without cloud.
    #r: impact parameter from center of shell
    #R: radius of shell (measured at R_type of rim)
    #dR: thickness of shell
    #R_type: Where in the rim R is measured at. 
    if R_type == "middle":
        pass
    elif R_type == "inner":
        R = R + dR / 2.
    elif R_type == "outer":
        R = R - dR / 2.
    else:
        raise Exception('radius_type must be one of "inner", "middle", or "outer".')

    if r < R - dR/2.:
        profile = 2*R * ((1 + dR/R)**2. - (r/R)**2.)**0.5 -\
                2*R * (1 - (r/R)**2.)**0.5

    elif r >= R and r < R + dR/2.:
        profile = 2*R * ((1 + dR/R)**2. - (r/R)**2.)**0.5

    elif r >= R + dR/2.:
        profile = 0

    return profile

def radial_profile(array, center=None, mode='average',
                   nbins=10, returnSEM=True, removeNaN=True):
    """Return a radial profile of `array`. The radial profile can
    be binned in `nbins` radial bins and then return a statistic given
    by `mode` such as the average in each radial bin. Or, if `mode` = 'alongx',
    returns the unbinned profile through the `center` of the array along the
    x-direction.
    
    Parameters
    ----------
    array : TYPE
        Description
    center : None, optional
        Description
    mode : str, optional
        Description
    nbins : int, optional
        Description
    returnSEM : bool, optional
        Description
    removeNaN : bool, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    
    Raises
    ------
    Exception
        Description
    """
    from scipy.stats import binned_statistic, sem
    import numpy as np
    y, x = np.indices(array.shape)
    if center is None:
        center = [(array.shape[1] - 1) / 2., (array.shape[0] - 1) / 2.]
    r = np.sqrt((center[0] - x) ** 2. + (center[1] - y) ** 2.)
    if removeNaN:
        r, array = r[np.isfinite(array)], array[np.isfinite(array)]
# Profile is the summed pixel values in each bin of radius.
# profile, rbins = np.histogram(r, weights=array, bins=nbins)
# print(r.max())
    if mode == 'average':
        # Divide by the number of pixels in each bin
        # profile = profile / np.histogram(r, bins=nbins)[0]
        profile, rbins, binnumber = binned_statistic(
            r, array, 'mean', bins=nbins)
    else:
        raise Exception("mode {} not implemented".format(mode))
    rbin_centers = (rbins[1:] + rbins[:-1]) / 2.
    if returnSEM:
        "Return the standard error on the mean of each bin."
        SEM = binned_statistic(
            r, array, lambda y: sem(y), bins=nbins)[0]
        return profile, SEM, rbin_centers
    else:
        return profile, rbin_centers
