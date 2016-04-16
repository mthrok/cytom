import logging

from cytom.logicle import logicle
from cytom.logicle import markertype

_LG = logging.getLogger(__name__)


def apply_logicle_transform(fcs, channels=None):
    """Apply logicle transformation to  FLUO measurements"""
    channel_info = fcs.meta['_channels_']
    _LG.debug('\n{}'.format(channel_info))
    if channels is None:
        channels = channel_info['$PnN']
    for name, range_ in zip(channels, map(int, channel_info['$PnR'])):
        markertype_ = markertype(name)
        if markertype_ == 'TIME':
            continue
        _LG.info('Processing: {}'.format(name))
        _LG.debug('Marker: {}, Type: {}, Range: {}'
                  .format(name, markertype_, range_))
        data = fcs.data[name].values  # numpy.ndarray
        if markertype_ == 'SCATTER':
            range_ = max(int(data.max()), range_)
            data = 4095.0 * data / range_
            data[4095 < data] = 4095
            data[data < 0] = 0
        else:
            _LG.info('Applying logicle transform.')
            range_ = max(int(data.max()), range_)
            data = 262144 * data / range_
            data = logicle(data, range_)
        fcs.data = fcs.data.assign(**{name: data})
