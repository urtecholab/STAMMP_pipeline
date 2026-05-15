#!/usr/bin/env python
# encoding: utf-8
# -*- mode: python; coding: utf-8 -*-
# vim: set fileencoding=utf-8:encoding=utf-8:syntax=python
import sys
import os
import errno
import argparse
from pathlib import Path

import pandas as pd

def main(argv=None):
    args = _cmdline(argv)
    samples = []
    for f in args.files:
        sample_id = Path(f).stem
        samples.append(pd.read_csv(f, sep='\t', index_col=0, names=['id', sample_id])[sample_id])
    tbl = pd.DataFrame(samples).fillna(0).T
    tbl = tbl.astype(int)
    tbl.index.name='taxon_id'
    tbl = tbl.reindex(index=tbl.index.sort_values(), columns=tbl.columns.sort_values())
    tbl.to_csv(sys.stdout)

    return 0  # success


def _cmdline(argv=None):
    """Parse commandline arguments"""
    p = argparse.ArgumentParser()
    p.add_argument(
        "-v", "--verbose", default=False, action="store_true", help="Verbose Logging"
    )
    p.add_argument("files", metavar="file", nargs="*", help="File(s) to process")


    return p.parse_args(argv)


def cli_wrap(cmd=main, argv=None):
    """
    Make sure to prevent exception tracebacks for common signals
    like SIGABORT (^C) or SIGPIPE (e.g. cmd | head)
    """
    try:
        status = cmd(argv)
    except KeyboardInterrupt:
        status = 1
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        status = 1
    except IOError as e:
        if e.errno == errno.EPIPE:
            status = 1
        else:
            raise
    return status

# ---------------
if __name__ == "__main__":
    raise SystemExit(cli_wrap(main))
