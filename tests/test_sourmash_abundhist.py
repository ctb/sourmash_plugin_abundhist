"""
Tests for sourmash_plugin_abundhist.
"""
import os
import pytest

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_run_sourmash(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('', fail_ok=True)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0                    # no args provided, ok ;)


def test_run_abundhist(runtmp):
    reads = utils.get_test_data('reads.sig.gz')

    runtmp.sourmash('scripts', 'abundhist', reads)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status == 0

    assert '65   [13596]  ****************************************' in runtmp.last_result.out


def test_run_abundhist_csv(runtmp):
    reads = utils.get_test_data('reads.sig.gz')

    csvfile = runtmp.output('zzz.csv')
    runtmp.sourmash('scripts', 'abundhist', reads, '--csv', csvfile)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status == 0

    assert '65   [13596]  ****************************************' in runtmp.last_result.out
    assert os.path.exists(csvfile)


def test_run_abundhist_fig(runtmp):
    reads = utils.get_test_data('reads.sig.gz')

    outfig = runtmp.output('fig.png')
    runtmp.sourmash('scripts', 'abundhist', reads,
                    '--max', '100', '--min', '1', '--bins', '100',
                    '--figure', outfig, '--ymax=200')

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status == 0

    assert '1    [10469]  ****************************************' in runtmp.last_result.out
    assert '35   [   25]  *' in runtmp.last_result.out

    assert os.path.exists(outfig)
