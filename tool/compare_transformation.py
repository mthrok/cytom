"""Script for checking if our implementation of logicle function is closer to
the one from FCSTRANS"""

import csv
import logging

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from extract_data import extract
from logicle import logicle

logger_ = logging.getLogger(__name__)


def parse_command_line_arguments():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument('original_file', help='original fcs file')
    ap.add_argument('converted_file', help='output file from FCSTrans')
    ap.add_argument('--x', default='Pacific Blue-A')
    ap.add_argument('--y', default='SSC-A')
    ap.add_argument('--debug', action='store_true')
    return ap.parse_args()


def init_logging(debug=False):
    level = logging.INFO
    if debug:
        level = logging.DEBUG
    logging.basicConfig(level=level)


def get_data_from_fcs(filename, x_attr, y_attr):
    """Extract data from FCS file"""
    data, channel_names = extract(filename)
    x, y = data[x_attr].values, data[y_attr].values
    return x, y


def get_data_from_csv(filename, x_attr, y_attr):
    """Extract data from CSV file"""
    with open(filename, 'r') as f:
        x, y = [], []
        ind_x, ind_y = None, None
        for line in csv.reader(f, delimiter='\t'):
            if ind_x is None:
                ind_x = line.index(x_attr)
                ind_y = line.index(y_attr)
            else:
                x.append(float(line[ind_x]))
                y.append(float(line[ind_y]))
        return x, y


def heatmap2d(ax, x, y, bins=50, scale='log'):
    """Draw heatmap on the given axis"""
    heatmap, xedges, yedges = np.histogram2d(y, x, bins=50)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    if scale == 'log':
        heatmap = np.log10(heatmap)
        bins = 'log'
        cb_label = 'Count (log10)'
    else:
        bins = None
        cb_label = 'Count (Linear)'
    ax.contour(heatmap, extent=extent)
    hb = ax.hexbin(x, y, bins=bins, cmap=cm.Greys)
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label(cb_label)


def main():
    args = parse_command_line_arguments()
    init_logging(args.debug)
    x0, y0 = get_data_from_fcs(args.original_file, args.x, args.y)
    x1, y1 = get_data_from_csv(args.converted_file, args.x, args.y)
    x2 = logicle(x0)
    fig = plt.figure()
    ax0 = fig.add_subplot(2, 2, 1)
    ax1 = fig.add_subplot(2, 2, 2)
    ax2 = fig.add_subplot(2, 2, 3)
    heatmap2d(ax0, x0, y0)
    heatmap2d(ax1, x2, y0)
    heatmap2d(ax2, x1, y1)
    ax0.set_xlabel(args.x)
    ax0.set_ylabel(args.y)
    ax0.set_title('Original FCS')
    ax1.set_xlabel(args.x)
    ax1.set_ylabel(args.y)
    ax1.set_title('Processed with our implementation')
    ax2.set_xlabel(args.x)
    ax2.set_ylabel(args.y)
    ax2.set_title('Processed by FCSTrans')

    logger_.info('Plot ready')
    plt.show()


if __name__ == '__main__':
    main()
