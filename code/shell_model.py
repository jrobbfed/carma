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
