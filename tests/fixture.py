import os
import gzip
import pandas as pd

from FlowCytometryTools import FCMeasurement


_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def load_fcs_data(filename):
    filepath = os.path.join(_DATA_DIR, filename)
    return FCMeasurement(ID=filename, datafile=filepath)


def load_fcstrans_data(filename):
    file_ = os.path.join(_DATA_DIR, filename)
    if file_.endswith('tar.gz'):
        file_ = gzip.GzipFile(file_)
    return pd.read_csv(file_, sep='\t')
