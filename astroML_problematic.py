import numpy as np

## AstroML files that didn't work with PyInstaller

def binned_statistic(x, values, statistic='mean',
                     bins=10, range=None):
    """
    Compute a binned statistic for a set of data.

    This is a generalization of a histogram function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values within each bin.

    Parameters
    ----------
    x : array_like
        A sequence of values to be binned.
    values : array_like
        The values on which the statistic will be computed.  This must be
        the same shape as x.
    statistic : string or callable, optional
        The statistic to compute (default is 'mean').
        The following statistics are available:

          * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
          * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
          * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
          * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
          * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.

    bins : int or sequence of scalars, optional
        If `bins` is an int, it defines the number of equal-width
        bins in the given range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    range : (float, float), optional
        The lower and upper range of the bins.  If not provided, range
        is simply ``(x.min(), x.max())``.  Values outside the range are
        ignored.

    Returns
    -------
    statistic : array
        The values of the selected statistic in each bin.
    bin_edges : array of dtype float
        Return the bin edges ``(length(statistic)+1)``.

    Notes
    -----
    All but the last (righthand-most) bin is half-open.  In other words, if
    `bins` is::

      [1, 2, 3, 4]

    then the first bin is ``[1, 2)`` (including 1, but excluding 2) and the
    second ``[2, 3)``.  The last bin, however, is ``[3, 4]``, which *includes*
    4.

    Examples
    --------
    >>> binned_statistic([1,2,1], bins=[0,1,2,3], 'count')
    (array([0, 2, 1]), array([0, 1, 2, 3]))\

    See Also
    --------
    np.histogram, binned_statistic_2d, binned_statistic_dd
    """
    try:
        N = len(bins)
    except TypeError:
        N = 1

    if N != 1:
        bins = [np.asarray(bins, float)]

    medians, edges = binned_statistic_dd([x], values, statistic,
                                         bins, range)

    return medians, edges[0]


def binned_statistic_2d(x, y, values, statistic='mean',
                        bins=10, range=None):
    """
    Compute a bidimensional binned statistic for a set of data.

    This is a generalization of a histogram2d function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values within each bin.

    Parameters
    ----------
    x : array_like
        A sequence of values to be binned along the first dimension.
    y : array_like
        A sequence of values to be binned along the second dimension.
    values : array_like
        The values on which the statistic will be computed.  This must be
        the same shape as x.
    statistic : string or callable, optional
        The statistic to compute (default is 'mean').
        The following statistics are available:

          * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
          * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
          * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
          * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
          * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.

    bins : int or [int, int] or array-like or [array, array], optional
        The bin specification:

          * the number of bins for the two dimensions (nx=ny=bins),
          * the number of bins in each dimension (nx, ny = bins),
          * the bin edges for the two dimensions (x_edges=y_edges=bins),
          * the bin edges in each dimension (x_edges, y_edges = bins).

    range : array_like, shape(2,2), optional
        The leftmost and rightmost edges of the bins along each dimension
        (if not specified explicitly in the `bins` parameters):
        [[xmin, xmax], [ymin, ymax]]. All values outside of this range will be
        considered outliers and not tallied in the histogram.

    Returns
    -------
    statistic : ndarray, shape(nx, ny)
        The values of the selected statistic in each two-dimensional bin
    xedges : ndarray, shape(nx + 1,)
      The bin edges along the first dimension.
    yedges : ndarray, shape(ny + 1,)
      The bin edges along the second dimension.

    See Also
    --------
    np.histogram2d, binned_statistic, binned_statistic_dd
    """

    # This code is based on np.histogram2d
    try:
        N = len(bins)
    except TypeError:
        N = 1

    if N != 1 and N != 2:
        xedges = yedges = np.asarray(bins, float)
        bins = [xedges, yedges]

    medians, edges = binned_statistic_dd([x, y], values, statistic,
                                         bins, range)

    return medians, edges[0], edges[1]


def binned_statistic_dd(sample, values, statistic='mean',
                        bins=10, range=None):
    """
    Compute a multidimensional binned statistic for a set of data.

    This is a generalization of a histogramdd function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values within each bin.

    Parameters
    ----------
    sample : array_like
        Data to histogram passed as a sequence of D arrays of length N, or
        as an (N,D) array.
    values : array_like
        The values on which the statistic will be computed.  This must be
        the same shape as x.
    statistic : string or callable, optional
        The statistic to compute (default is 'mean').
        The following statistics are available:

          * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
          * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
          * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
          * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
          * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.

    bins : sequence or int, optional
        The bin specification:

          * A sequence of arrays describing the bin edges along each dimension.
          * The number of bins for each dimension (nx, ny, ... =bins)
          * The number of bins for all dimensions (nx=ny=...=bins).

    range : sequence, optional
        A sequence of lower and upper bin edges to be used if the edges are
        not given explicitely in `bins`. Defaults to the minimum and maximum
        values along each dimension.

    Returns
    -------
    statistic : ndarray, shape(nx1, nx2, nx3,...)
        The values of the selected statistic in each two-dimensional bin
    edges : list of ndarrays
        A list of D arrays describing the (nxi + 1) bin edges for each
        dimension

    See Also
    --------
    np.histogramdd, binned_statistic, binned_statistic_2d
    """
    if type(statistic) == str:
        if statistic not in ['mean', 'median', 'count', 'sum']:
            raise ValueError('unrecognized statistic "%s"' % statistic)
    elif callable(statistic):
        pass
    else:
        raise ValueError("statistic not understood")

    # This code is based on np.histogramdd
    try:
        # Sample is an ND-array.
        N, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = np.atleast_2d(sample).T
        N, D = sample.shape

    nbin = np.empty(D, int)
    edges = D * [None]
    dedges = D * [None]

    try:
        M = len(bins)
        if M != D:
            raise AttributeError('The dimension of bins must be equal '
                                 'to the dimension of the sample x.')
    except TypeError:
        bins = D * [bins]

    # Select range for each dimension
    # Used only if number of bins is given.
    if range is None:
        smin = np.atleast_1d(np.array(sample.min(0), float))
        smax = np.atleast_1d(np.array(sample.max(0), float))
    else:
        smin = np.zeros(D)
        smax = np.zeros(D)
        for i in np.arange(D):
            smin[i], smax[i] = range[i]

    # Make sure the bins have a finite width.
    for i in np.arange(len(smin)):
        if smin[i] == smax[i]:
            smin[i] = smin[i] - .5
            smax[i] = smax[i] + .5

    # Create edge arrays
    for i in np.arange(D):
        if np.isscalar(bins[i]):
            nbin[i] = bins[i] + 2  # +2 for outlier bins
            edges[i] = np.linspace(smin[i], smax[i], nbin[i] - 1)
        else:
            edges[i] = np.asarray(bins[i], float)
            nbin[i] = len(edges[i]) + 1  # +1 for outlier bins
        dedges[i] = np.diff(edges[i])

    nbin = np.asarray(nbin)

    # Compute the bin number each sample falls into.
    Ncount = {}
    for i in np.arange(D):
        Ncount[i] = np.digitize(sample[:, i], edges[i])

    # Using digitize, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    outliers = np.zeros(N, int)
    for i in np.arange(D):
        # Rounding precision
        decimal = int(-np.log10(dedges[i].min())) + 6
        # Find which points are on the rightmost edge.
        on_edge = np.where(np.around(sample[:, i], decimal)
                           == np.around(edges[i][-1], decimal))[0]
        # Shift these points one bin to the left.
        Ncount[i][on_edge] -= 1

    # Compute the sample indices in the flattened statistic matrix.
    ni = nbin.argsort()
    shape = []
    xy = np.zeros(N, int)
    for i in np.arange(0, D - 1):
        xy += Ncount[ni[i]] * nbin[ni[i + 1:]].prod()
    xy += Ncount[ni[-1]]

    result = np.empty(nbin.prod(), float)

    if statistic == 'mean':
        result.fill(np.nan)
        flatcount = np.bincount(xy, None)
        flatsum = np.bincount(xy, values)
        a = np.arange(len(flatcount))
        result[a] = flatsum
        result[a] /= flatcount
    elif statistic == 'count':
        result.fill(0)
        flatcount = np.bincount(xy, None)
        a = np.arange(len(flatcount))
        result[a] = flatcount
    elif statistic == 'sum':
        result.fill(0)
        flatsum = np.bincount(xy, values)
        a = np.arange(len(flatsum))
        result[a] = flatsum
    elif statistic == 'median':
        result.fill(np.nan)
        for i in np.unique(xy):
            result[i] = np.median(values[xy == i])
    elif callable(statistic):
        try:
            null = statistic([])
        except:
            null = np.nan
        result.fill(null)
        for i in np.unique(xy):
            result[i] = statistic(values[xy == i])

    # Shape into a proper matrix
    result = result.reshape(np.sort(nbin))
    for i in np.arange(nbin.size):
        j = ni.argsort()[i]
        result = result.swapaxes(i, j)
        ni[i], ni[j] = ni[j], ni[i]

    # Remove outliers (indices 0 and -1 for each dimension).
    core = D * [slice(1, -1)]
    result = result[core]

    if (result.shape != nbin - 2).any():
        raise RuntimeError('Internal Shape Error')
    return result, edges

sigmaG_factor = 0.74130110925280102
def sigmaG(a, axis=None, overwrite_input=False, keepdims=False):
    """Compute the rank-based estimate of the standard deviation
    Parameters
    ----------
    a : array_like
        Array containing numbers whose mean is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the means are computed. The default is to compute
        the mean of the flattened array.
    overwrite_input : bool, optional
       If True, then allow use of memory of input array `a` for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted.
       Default is False. Note that, if `overwrite_input` is True and the
       input is not already an array, an error will be raised.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.
    Returns
    -------
    median : ndarray, see dtype parameter above
        array containing the median values
    sigmaG : ndarray, see dtype parameter above.
        array containing the robust estimator of the standard deviation
    See Also
    --------
    median_sigmaG : robust rank-based estimate of mean and standard deviation
    Notes
    -----
    This routine uses a single call to ``np.percentile`` to find the
    quartiles along the given axis, and uses these to compute the
    sigmaG, a robust estimate of the standard deviation sigma:
    sigmaG = 0.7413 * (q75 - q25)
    where 0.7413 ~ 1 / (2 sqrt(2) erf^-1(0.5))
    """
    q25, q75 = np.percentile(a, [25, 75],
                             axis=axis,
                             overwrite_input=overwrite_input)
    sigmaG = sigmaG_factor * (q75 - q25)

    if keepdims:
        if axis is None:
            newshape = a.ndim * (1,)
        else:
            newshape = np.asarray(a.shape)
            newshape[axis] = 1

        sigmaG = sigmaG.reshape(newshape)

    return sigmaG

def convert_to_stdev(logL):
    """
    Given a grid of log-likelihood values, convert them to cumulative
    standard deviation.  This is useful for drawing contours from a
    grid of likelihoods.
    """
    sigma = np.exp(logL)

    shape = sigma.shape
    sigma = sigma.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(sigma)[::-1]
    i_unsort = np.argsort(i_sort)

    sigma_cumsum = sigma[i_sort].cumsum()
    sigma_cumsum /= sigma_cumsum[-1]

    return sigma_cumsum[i_unsort].reshape(shape)
