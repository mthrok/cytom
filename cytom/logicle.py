from __future__ import division

import logging
from sys import float_info

import numpy as np
from scipy.optimize import brenth

_LG = logging.getLogger(__name__)


###############################################################################
def markertype(channel_name):
    """Convert channel name into marker type."""
    channel_name = channel_name.upper()
    if channel_name[:2] in ('FS', 'SS'):
        return 'SCATTER'
    if channel_name.startswith('TIME'):
        return 'TIME'
    if channel_name.endswith('-A'):
        return 'FLUO_AREA'
    return 'FLUO_NON_AREA'


def biexponential_(x, a, b, c, d, f, w):
    return a * np.exp(b * (x - w)) - c * np.exp(-d * (x - w)) + f


def biexponential_transform(X, a, b, c, d, f, w, max_ite=5000):
    # TODO: Find a faster way
    Y = np.zeros(X.shape)
    for i in range(len(X)):
        x = X[i]
        # We can probably get rid of this part.
        step = 0.5
        for _ in range(max_ite):
            step *= 1.5
            bf1 = biexponential_(-step, a, b, c, d, f, w)
            bf2 = biexponential_(step, a, b, c, d, f, w)
            if (bf1 - x) * (bf2 - x) <= 0:
                break
        else:
            _LG.info('Failed to find range.')
        Y[i] = brenth(lambda x_: biexponential_(x_, a, b, c, d, f, w)-x, -step, step)
    return Y


def logicle_(X, w, r, d, scale):
    def _func(p, w):  # TODO: Give appropriate name
        return -w + 2*p*np.log(p)/(p+1)

    d = d * np.log(10)
    scale = scale / d
    p = brenth(_func, 0.99, 1000, args=(w,))  # If w > 6, then this will raise
    a = r * np.exp(-(d-w))
    b = 1
    c = a * p * p
    d = 1 / p
    f = a * (p*p - 1)
    # Apply biexponential transform here.
    y = biexponential_transform(X, a, b, c, d, f, w)
    y = y * scale
    y[y < 0] = 0
    return y


def logicle_w(w, cutoff=-111, r=262144, d=4.5, scale=1):
    if w > d:
        w = d
    x = np.asarray([cutoff])
    return logicle_(x, w, r, d, scale)[0] - float_info.epsilon**0.6


def logicle(data, r=262144, d=4.5, scale=4096, cutoff=-111, w=-1):
    """Apply logicle transformation to data

    Args:
      data(NumPy NDArray): 1D data
      r(int)  : The upper bound of data. E.g., 10,000 for
        common 4 decade data and 262,144 for an 18 bit data.
      d(float): The total decade for displaying.
      scale(float): The upper range after transformation.
      cutoff(float): The minimum value to retain. Data points smaller than
        this are truncated.
      w(float): The decade of linear region, including negative region.
        If negative value is given, automated is used.
    """
    if w > d:
        raise ValueError(
            'Negative range decades(w) must be '
            'smaller than total number of decades(d)')
    if w < 0:
        w = brenth(logicle_w, 0, d, args=(cutoff, ))
    _LG.debug('params: r={}, d={}, scale={}, cutoff={}, w={}'
              .format(r, d, scale, cutoff, w))
    return logicle_(data, w, r, d, scale)
