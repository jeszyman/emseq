#!/usr/bin/env python3
"""Fraction of a control's fragments to draw so its depth per CpG matches a target table.

Depth per CpG is the number of rows of a JEEPERS merge table whose column 11
lists the CpG. The depth of a table is the median of that number over the
CpGs that at least one fragment covers. The fraction is the target table's
depth divided by the control table's depth.

Prints one tab-separated line: target depth, control depth, fraction.
"""
from __future__ import annotations

import argparse
import collections
import statistics


def cpg_depths(lines) -> collections.Counter:
    """{CpG position: number of rows listing it} from merge-table lines."""
    depth: collections.Counter = collections.Counter()
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) > 10 and f[10] not in ("-", ""):
            depth.update(f[10].split(","))
    return depth


def median_depth(path: str) -> float:
    with open(path) as fh:
        depth = cpg_depths(fh)
    return statistics.median(depth.values())


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--target-merge", required=True)
    p.add_argument("--control-merge", required=True)
    a = p.parse_args(argv)
    target, control = median_depth(a.target_merge), median_depth(a.control_merge)
    print(f"{target:g}\t{control:g}\t{target / control:.8f}")


if __name__ == "__main__":
    main()
