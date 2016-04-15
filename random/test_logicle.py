from FlowCytometryTools import FCMeasurement

import logging

from cytom.compensate import compensate
from cytom.logicle import apply_logicle_transform

_LG = logging.getLogger(__name__)


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
