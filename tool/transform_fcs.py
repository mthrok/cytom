
import logging
from argparse import ArgumentParser as AP

from cytom.tool import load_fcs
from cytom.tool import apply_logicle_transform

_LG = logging.getLogger(__name__)


def parse_command_line_arguments():
    ap = AP('Apply logicle transformatoin to FCS file and dump to pickle obj')
    ap.add_argument('input', help='Input FCS file path.')
    ap.add_argument('output', help='Output FCS(Pickled) file path.')
    return ap.parse_args()


def main():
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line_arguments()

    _LG.info('Loading {}'.format(args.input))
    fcs = load_fcs(args.input)

    _LG.info('Applying logicle transform.')
    # Currently applying logicle transform takes a lot of time, and we want to
    # run this script inside CircleCI without timeout. So we put the
    # processing into subprocess and periodically print status.
    run_and_check(lambda: apply_logicle_transform(fcs,), 300)

    _LG.info('Saving the result to {}'.format(args.output))
    fcs.save(args.output)


if __name__ == '__main__':
    main()
