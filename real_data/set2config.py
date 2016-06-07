#!/usr/bin/env python3

from collections import defaultdict
import sys
import json
from os import path

if len(sys.argv) < 2:
    print("USAGE: set2config.py SETFILE ... >sets.json", file=sys.stderr)
    exit(1)

sets = defaultdict(dict)

for setfile in sys.argv[1:]:
    setname = path.splitext(path.basename(setfile))[0]
    project = path.basename(path.dirname(setfile))
    if project == 'sra':
        project = ''
    runs = []
    with open(setfile) as fh:
        for line in fh:
            run = line.strip()
            if run:
                runs.append(run)
    sets[project][setname] = runs
print(json.dumps(sets, indent=2))
