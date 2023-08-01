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

import numpy, collections, csv
import termplotlib as tpl
import seaborn
import matplotlib.pyplot as plt


###

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
        add_ksize_arg(subparser, default=31)
        add_moltype_args(subparser)


    def main(self, args):
        """
        output abundance histogram and/or raw abundances.
        """

        set_quiet(args.quiet)
        moltype = sourmash_args.calculate_moltype(args)

        outlist = []
        total_loaded = 0
        for filename in args.signature_file:
            siglist = sourmash.load_file_as_signatures(filename, ksize=args.ksize,
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

        counts_d = collections.defaultdict(int)
        for ss in siglist:
            for hashval, abund in ss.minhash.hashes.items():
                counts_d[hashval] += abund

        all_counts = list(counts_d.values())

        min_range = 1
        if args.min is not None:
            min_range = args.min
        max_range = max(all_counts)
        if args.max is not None:
            max_range = args.max

        n_bins = args.bins
        if max_range - min_range + 1 < n_bins:
            n_bins = max_range - min_range + 1

        # make hist
        counts, bin_edges = numpy.histogram(all_counts,
                                            range=(min_range, max_range),
                                            bins=n_bins)
        bin_edges = bin_edges.astype(int)

        # plot
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
            seaborn.histplot(all_counts, binrange=(min_range, max_range),
                             bins=n_bins, kde=True)
            if args.ymax:
                plt.ylim(top=args.ymax)
            plt.savefig(args.figure)
