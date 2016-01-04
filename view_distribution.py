#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from FlowCytometryTools import FCMeasurement


def analyze(datafile):
    sample = FCMeasurement(ID=datafile, datafile=datafile)

    hists = []
    channels = []
    bin_edges = []
    for channel, data in zip(sample.channel_names, sample.data.values.T):
        if channel == 'Time':
            continue
        print("{}: Min: {}, Max: {}".format(channel, np.min(data), np.max(data)))
        # Convert scale
        if channel.startswith('FSC') or channel.startswith('SSC'):
            pass
        else:
            data[data <= 0] = 1  # For convenience until biexponential is implemented
            data = np.log10(data)
        # modify range
        if channel == 'FSC-W':
            data = data[data < 1.5e5]
            pass
        hist, bin_edge = np.histogram(data.flatten(), bins=100)
        hists.append(hist)
        channels.append(channel)
        bin_edges.append(bin_edge)

    fig = plt.figure()
    i = 0
    for channel, hist, bin_edge in zip(channels, hists, bin_edges):
        i += 1
        ax = fig.add_subplot(len(sample.channel_names)/4, 4, i)
        ax.bar(bin_edge[:-1], hist/1000, width=bin_edge[1:]-bin_edge[:-1])
        ax.ticklabel_format(axis='x', style='sci', scilimits=(1, 4))
        ax.set_title(channel)
    plt.show()


def parse_command_line_arguments():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument('input_file', help='fcs file')
    return ap.parse_args()


def main():
    args = parse_command_line_arguments()
    analyze(args.input_file)

if __name__ == '__main__':
    main()
