#!/usr/bin/env python3
from __future__ import print_function, absolute_import, division
import screed
import docopt

def seq2fa(name, seq, linelen=80):
    lines = ['>{}'.format(name), ]

    for start in range(0, len(seq), linelen):
        lines.append(seq[start:start + linelen])

    return '\n'.join(lines) + '\n'

CLI = """
USAGE:
    splitfa <fasta> <prefix> [<name> ...]
"""

opts = docopt.docopt(CLI)
prefix = opts['<prefix>']
names = opts['<name>']

orig_names = []
with screed.open(opts['<fasta>']) as fh:
    for record in fh:
        orig_names.append(int(record.name))

newnames = {}
for orig, new in zip(sorted(orig_names), names):
    newnames[str(orig)] = new

with screed.open(opts['<fasta>']) as fh:
    for record in fh:
        name = newnames[record.name]
        fname = "{}{}.fasta".format(prefix, name)
        with open(fname, 'w') as ofh:
            seq = str(record.sequence).translate({'-': ''})
            print(seq2fa(name, seq), file=ofh)

