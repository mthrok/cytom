from FlowCytometryTools import FCMeasurement


def load_fcs(filepath, ID=None):
    if ID is None:
        ID = filepath
    return FCMeasurement(ID=ID, datafile=filepath)
