from FlowCytometryTools import FCMeasurement

import logging

from logicle import logicle
from logicle import markertype
from compensate import compensate

_LG = logging.getLogger(__name__)


def apply_logicle_transform(fcs):
    """Apply logicle transformation to  FLUO measurements"""
    channel_info = fcs.meta['_channels_']
    _LG.info('\n{}'.format(channel_info))
    names, ranges = channel_info['$PnN'], map(int, channel_info['$PnR'])
    for name, range_ in zip(names, ranges):
        markertype_ = markertype(name)
        _LG.debug('Marker: {}, Type: {}, Range: {}'
                  .format(name, markertype_, range_))
        if markertype_ == 'TIME':
            continue
        data = fcs.data[name].values  # numpy.ndarray
        if markertype_ == 'SCATTER':
            range_ = max(int(data.max()), range_)
            data = 4095.0 * data / range_
            data[4095 < data] = 4095
            data[data < 0] = 0
        else:
            range_ = max(int(data.max()), range_)
            data = 262144 * data / range_
            data = logicle(data, range_)
        fcs.data = fcs.data.assign(**{name: data})


def _process_fcs_file(filename):
    # TODO: Remove this script
    fcs = FCMeasurement(ID='Test Sample', datafile=filename)
    fcs_format = fcs.meta['__header__']['FCS format']
    if int(float(fcs_format[3:])) not in [2, 3]:
        raise RuntimeError('{} is not supported.'.format(fcs_format))
    if fcs.meta['$DATATYPE'] == 'I':
        _LG.info('Not processing integer type')
        return fcs
    compensate(fcs)
    apply_logicle_transform(fcs)


def _parse_command_line_arguments():
    from argparse import ArgumentParser
    ap = ArgumentParser('Apply logicle transformation on fcs file')
    ap.add_argument('input')
    ap.add_argument('--output')
    return ap.parse_args()


def init_logging(debug=True):
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


def main():
    args = _parse_command_line_arguments()
    init_logging()
    _process_fcs_file(args.input)


if __name__ == '__main__':
    main()
