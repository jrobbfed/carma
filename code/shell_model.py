import numpy as np

def analytic_profile(r, R, dR, radius_type="middle"):
    #From Beaumont & Williams (2010), assumes bare optically thin shell without cloud.
    #r: impact parameter from center of shell
    #R: radius of shell (measured at R_type of rim)
    #dR: thickness of shell
    #R_type: Where in the rim R is measured at. 
    
    if radius_type == "middle":
        R = R - dR / 2.
    elif radius_type == "inner":
        pass
    elif radius_type == "outer":
        R = R - dR
    else:
        raise Exception('radius_type must be one of "inner", "middle", or "outer".')

    r = np.asarray(r)
    profile = np.zeros_like(r)
    inside = (r < R)
    #print(inside)
    profile[inside] = \
                2*R * ((1 + dR/R)**2. - (r[inside]/R)**2.)**0.5 -\
                2*R * (1 - (r[inside]/R)**2.)**0.5

    on = (r >= R) & (r < R + dR)
    #print(on)
    profile[on] = \
                2*R * ((1 + dR/R)**2. - (r[on]/R)**2.)**0.5

    outside = (r >= R + dR)
    #print(outside)
    profile[outside] = 0

    return profile

def radial_profile(array, center=None, mode='average',
                   nbins=10, returnSEM=True, removeNaN=True):
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
        rbin_centers = (rbins[1:] + rbins[:-1]) / 2.
    elif mode == 'along x':
        print(np.shape(array), np.shape(x), np.shape(y))
        print(y, np.floor(center[1])) 
        profile = array[(y == np.floor(center[1])) & (x - np.floor(center[0]) >= 0)]
        rbin_centers = r[(y == np.floor(center[1])) & (x - np.floor(center[0]) >= 0)]
        returnSEM = False #wouldn't make sense in this context.
    else:
        raise Exception("mode {} not implemented".format(mode))

    if returnSEM:
        "Return the standard error on the mean of each bin."
        SEM = binned_statistic(
            r, array, lambda y: sem(y), bins=nbins)[0]
        return profile, SEM, rbin_centers
    else:
        return profile, rbin_centers
