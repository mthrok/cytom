import unittest
from time import time
from multiprocessing import Process

import numpy as np

from cytom.compensate import compensate
from cytom.tool import apply_logicle_transform
from tests.fixture import load_fcs_data
from tests.fixture import load_fcstrans_data

import logging
_LG = logging.getLogger(__name__)

# Since the current loggicle function is so slow so as not to timeout in CircleCI
logging.basicConfig(level=logging.INFO)
logging.getLogger('cytom.logicle').setLevel(logging.DEBUG)


def run_and_check(task, interval):
    """Run task on subprocess"""
    t = Process(target=task)
    t.start()
    t0 = time()
    try:
        while t.is_alive():
            t.join(interval)
            _LG.info('Still processing... {} sec'.format(time() - t0))
    except BaseException:
        _LG.exception('')
        t.terminate()
        raise


def load_test_data(original_file, transed_file):
    fcs = load_fcs_data(original_file)
    fcs_trans = load_fcstrans_data(transed_file)
    return fcs, fcs_trans


class TestLogicleImplementation(unittest.TestCase):
    def base_logicle_test(self, original_file, transed_file):
        diff_threshold = 1.0
        fcs, fcs_trans = load_test_data(original_file, transed_file)
        compensate(fcs)
        # Currently applying logicle transform takes a lot of time, and
        # we want to run this script inside CircleCI without timeout.
        # So we put the processing into subprocess and periodically print
        # status.
        run_and_check(lambda: apply_logicle_transform(fcs,), 300)

        channel_info = fcs.meta['_channels_']
        result = {}
        fail = False
        for channel in channel_info['$PnN']:
            if channel == 'Time':
                continue
            val0 = fcs[channel].values
            val1 = fcs_trans[channel].values
            denom = (val0 + val1) / 2 + 1
            numer = np.abs(val0 - val1)
            diff = 100 * np.mean(numer/denom)
            result[channel] = diff
            if diff > diff_threshold:
                fail = True
        _LG.info(result)
        self.assertTrue(not fail)

    def test_logicle1(self):
        self.base_logicle_test('1-1.fcs', '1-1_transformed.txt')

    def test_logicle2(self):
        self.base_logicle_test('2-4.fcs', '2-4_transformed.txt')


if __name__ == '__main__':
    unittest.main()
