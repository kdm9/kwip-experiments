#!/usr/bin/env python3
from __future__ import print_function, absolute_import, division
import screed
import docopt

CLI = """
USAGE:
    fasplit <fasta> <prefix>
"""

opts = docopt.docopt(CLI)
prefix = opts['<prefix>']

with screed.open(opts['<fasta>']) as fh:
    for record in fh:
        fname = "{}{}.fasta".format(prefix, record.name)
        with open(fname, 'w') as ofh:
            print(">", record.name, sep='', file=ofh)
            print(str(record.sequence).translate({'-': ''}), file=ofh)

