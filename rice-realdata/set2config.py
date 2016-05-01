#!/usr/bin/env python3

import sys
import json
from os import path


def set2dict(setfile):
    setname = path.splitext(path.basename(setfile))[0]
    runs = []
    with open(setfile) as fh:
        for line in fh:
            run = line.strip()
            if run:
                runs.append(run)
    return {setname: runs}


if __name__ == "__main__":
    sets = {}
    for sf in sys.argv[1:]:
        sets.update(set2dict(sf))
    print(json.dumps(sets, indent=2))
