import logging

import numpy as np
from numpy.linalg import inv

_LG = logging.getLogger(__name__)


def parse_compensation_matrix(fcs):
    """Parse compensation matrix from FCS file
    Args:
      fcs (FCMeasurement): Input FCS data
    Returns:
      NumPy NdArray: Compensation matrix
      list: Channel names
    """
    row_spill = fcs.meta['SPILL'].split(',')
    dim = int(row_spill[0])
    channel_names = row_spill[1:1+dim]
    mat = np.array(list(map(float, row_spill[1+dim:]))).reshape((dim, dim))
    return inv(mat.T), channel_names


def compensate(fcs):
    """Compensate FCS measurement data if SPILL is available. In-place op.
    Args:
      fcs (FCMeasurement): Input FCS data
    """
    if 'SPILL' not in fcs.meta:
        _LG.info('Not compensating since `SPILL` is missing.')
        return
    _LG.info('Compensating samples')
    try:
        comp_mat, channels, = parse_compensation_matrix(fcs)
        vals = np.dot(fcs.data[channels], comp_mat.T)
        vals_ = {channel: vals[:, i] for i, channel in enumerate(channels)}
        fcs.data = fcs.data.assign(**vals_)
    except Exception:
        _LG.exception('Failed to compensate')
