import os
import unittest

import numpy as np

from cytom.compensate import compensate
from cytom.logicle import apply_logicle_transform
from tests.fixture import load_fcs_data
from tests.fixture import load_fcstrans_data

import logging
logging.basicConfig(level=logging.WARN)
logging.getLogger('cytom.logicle').setLevel(logging.DEBUG)
# Since the current loggicle function is so slowso as not to timeout in CircleCI


def load_test_data():
    base_dir = '2015_12_16_full_staining'
    original = os.path.join(base_dir, '2015_12_16_full_staining.fcs')
    converted = os.path.join(base_dir, '2015_12_16_full_staining.txt.tar.gz')
    fcs = load_fcs_data(original)
    fcs_trans = load_fcstrans_data(converted)
    return fcs, fcs_trans


class TestLogicleImplementation(unittest.TestCase):
    def test_logicle(self):
        diff_threshold = 1.0
        fcs, fcs_trans = load_test_data()
        compensate(fcs)
        channel_info = fcs.meta['_channels_']
        for channel in channel_info['$PnN']:
            if channel == 'Time':
                continue
            apply_logicle_transform(fcs, channels=[channel])
            val0 = fcs[channel].values
            val1 = fcs_trans[channel].values
            denom = (val0 + val1) / 2 + 1
            numer = np.abs(val0 - val1)
            diff = 100 * np.mean(numer/denom)
            self.assertTrue(
                diff < diff_threshold,
                'Diff exceeded for {} by {}'.format(channel, diff))

if __name__ == '__main__':
    unittest.main()
