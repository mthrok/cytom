from FlowCytometryTools import FCMeasurement


def extract(datafile):
    sample = FCMeasurement(ID='Test Sample', datafile=datafile)
    print sample.channel_names
    return (sample.data, sample.channel_names)


def parse_command_line_arguments():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument('input_file', help='fcs file')
    return ap.parse_args()


def main():
    args = parse_command_line_arguments()
    sample = extract(args.input_file)
    print sample.data
    print sample.channel_names

if __name__ == '__main__':
    main()
