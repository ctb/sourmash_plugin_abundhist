"""\
calculate abundance profiles from one or more abund sketches

Display histograms of k-mer/hash multiplicity in sourmash sketches created
with `-p abund`.

'abundhist' provides text output, as well as CSV output.
"""

usage="""
   sourmash scripts abundhist abund-sketch.sig.gz
"""

epilog="""
See https://github.com/ctb/sourmash_plugin_abundhist for more examples.

Need help? Have questions? Ask at http://github.com/sourmash/issues!
"""

import argparse
import sourmash

from sourmash import sourmash_args
from sourmash.cli.utils import add_ksize_arg, add_moltype_args
from sourmash.logging import debug_literal, set_quiet, notify
from sourmash.plugins import CommandLinePlugin

import numpy as np
import collections, csv
import termplotlib as tpl
import seaborn
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks

###

def find_rightmost_peak(abunds, counter):
    abunds = list(abunds)
    counts_dist = list(counter.items())
    sum_hist = sum(counter.values())

    sofar = 0
    for k, v in sorted(counts_dist):
        sofar += v
        if sofar >= 0.99*sum_hist:
            max_range = k
            break

    print(f"max_range: {max_range}, sum_hist: {sum_hist}, sofar: {sofar}")
    print(max_range, len(abunds), max(abunds))

    arr = np.array(abunds)
    arr = arr[(arr >= 5)].reshape(-1, 1)

    kde = KernelDensity(kernel="gaussian", bandwidth=20).fit(arr)

    x = np.linspace(0, max_range, max_range)
    y = kde.score_samples(x.reshape(-1, 1))

    yl = np.exp(y)

    rightmost = 0
    xx = find_peaks(yl, width=20)[0]
    print(xx)
    for i in xx:
        xval = x[i]
        if xval > max_range:
            print('XXX2', xval, max_range)
            break
        rightmost = xval

    print(f'right most peak: {rightmost}')
    return rightmost, yl[i]

#
# CLI plugin - supports 'sourmash scripts abundhist'
#

class Command_Abundhist(CommandLinePlugin):
    command = 'abundhist'
    description = __doc__
    usage = usage
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        subparser.add_argument('signature_file', nargs='+')

        subparser.add_argument(
            '--csv', metavar='FILE',
            help='output histogram to this file (in CSV format)'
        )
        subparser.add_argument('--figure', help='save figure to this file')
        subparser.add_argument(
            '--abundances-csv', metavar='FILE',
            help='output hashes and abundances to this file (in CSV format)')
        subparser.add_argument(
            '--md5', default=None,
            help='select signatures whose md5 contains this substring'
        )
        subparser.add_argument(
            '--name', default=None,
            help='select signatures whose name contains this substring'
        )
        subparser.add_argument(
            '-M', '--max', type=int, default=None,
            help='max value for histogram range (default none)')
        subparser.add_argument(
            '-m', '--min', type=int, default=None,
            help='min value for histogram range (default none)')
        subparser.add_argument(
            '--bins', type=int, default=10,
            help='number of bins (default 10)')
        subparser.add_argument('--ymax', type=int,
                               help='maximum Y value for histogram display')
        subparser.add_argument('--figure-title',
                               default=None, help="plot title")
        subparser.add_argument('-I', '--intersect',
                               help='plot only hashes that intersect with this signature')
        subparser.add_argument('--silent', action='store_true',
                               help="do not output histogram to text")
        add_ksize_arg(subparser, default=31)
        add_moltype_args(subparser)


    def main(self, args):
        """
        output abundance histogram and/or raw abundances.
        """

        set_quiet(args.quiet)
        moltype = sourmash_args.calculate_moltype(args)
        ksize = args.ksize

        outlist = []
        total_loaded = 0
        for filename in args.signature_file:
            siglist = sourmash.load_file_as_signatures(filename, ksize=ksize,
                                                       select_moltype=moltype)
            siglist = list(siglist)

            total_loaded += len(siglist)

            # select!
            if args.md5 is not None:
                siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
            if args.name is not None:
                siglist = [ ss for ss in siglist if args.name in ss.name() ]

        notify("loaded {} total that matched ksize & molecule type",
               total_loaded)
        if len(siglist) != total_loaded:
            notify("selected {} via name / md5 selectors".format(len(siglist)))
        notify('')

        # do we need to intersect hashes?
        intersect_hashes = None
        if args.intersect:
            notify(f"loading --intersect sketch from '{args.intersect}' k={ksize} moltype={moltype}")
            ss = sourmash.load_file_as_signatures(args.intersect,
                                                  ksize=ksize,
                                                  select_moltype=moltype)
            ss = list(ss)
            if len(ss) == 0:
                notify("ERROR: cannot find a sketch that matches ksize/moltype")
                sys.exit(-1)
            elif len(ss) > 1:
                notify("ERROR: find {len(ss)} sketches that match ksize/moltype")
                sys.exit(-1)

            ss = ss[0]

            intersect_hashes = set(ss.minhash.hashes)

        # track hashval abundanc distribution
        counts_d = collections.defaultdict(int)
        # track abundance distribution, too
        counter = collections.Counter()

        # collect across all input signatures
        for ss in siglist:
            hashes_d = ss.minhash.hashes
            hashvals = set(hashes_d)
            if intersect_hashes is not None:
                hashvals &= intersect_hashes

            for hashval in hashvals:
                abund = hashes_d[hashval]
                counts_d[hashval] += abund
                counter[abund] += 1

        all_counts = list(counts_d.values())
        counts_dist = list(counter.items())
        sum_hist = sum(counter.values())
        max_range = max(counter.values())

        # find count that covers 95% of distribution.
        sofar = 0
        for k, v in sorted(counts_dist):
            sofar += v
            if sofar >= 0.99*sum_hist:
                max_range = 2*k
                print(f'setting default max_range to {max_range} (2x 99% of counts)')
                break

        if args.max is not None:
            print(f'overriding default max_range with -M {args.max}')
            max_range = args.max

        min_range = 1
        if args.min is not None:
            min_range = args.min

        n_bins = args.bins
        print(f"set number of bins to {n_bins} (--bins)")
        if max_range - min_range + 1 < n_bins:
            n_bins = max_range - min_range + 1
            print(f"reducing to {n_bins} because of max/min range")

        ###
        rightmost_x, rightmost_y = find_rightmost_peak(counts_d.values(), counter)
        ###

        # make hist
        counts, bin_edges = np.histogram(all_counts,
                                            range=(min_range, max_range),
                                            bins=n_bins)
        counts_bin3 = counts[2]

        bin_edges = bin_edges.astype(int)

        # plot
        if not args.silent:
            fig = tpl.figure()
            f = fig.barh(counts, [ str(x) for x in bin_edges[1:] ], force_ascii=True)
            fig.show()

        # output histogram in csv?
        if args.csv:
            with sourmash_args.FileOutput(args.csv, 'wt') as fp:
                w = csv.writer(fp)
                w.writerow(['count', 'n_count'])
                for nc, c in zip(counts, bin_edges[1:]):
                    w.writerow([c, nc])

        # output raw counts tagged with hashval?
        if args.abundances_csv:
            with sourmash_args.FileOutput(args.abundances_csv, 'wt') as fp:
                w = csv.writer(fp)
                w.writerow(['hashval', 'count'])
                for hashval, count in counts_d.items():
                    w.writerow([hashval, count])

        # output figure?
        if args.figure:
            seaborn.histplot(all_counts,
                             binrange=(min_range, max_range),
                             bins=n_bins, kde=True)
            if args.ymax:
                plt.ylim(top=args.ymax)
            else:
                plt.ylim(top=counts_bin3)

            if args.figure_title:
                plt.title(args.figure_title)
            else:
                plt.title("K-mer abundance histogram")
            plt.xlim(min_range, max_range)

            _, max_y = plt.ylim()
            plt.plot([rightmost_x, rightmost_x,], [0, max_y], 'o--')
            plt.xlabel("k-mer abundance")
            plt.ylabel("N(k-mers at that abundance)")
            plt.savefig(args.figure)

            print('XYZ', rightmost_x, rightmost_y)
