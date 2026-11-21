#!/usr/bin/env python3
"""Pool the draws of each spike-in control and compare mid-fragment methylation with the 3' bins.

A control subsampled to plasma depth holds too few fragments in one draw to
resolve a 5-point difference, so the fragments of all draws are pooled per
control and setting before the mid-fragment minus bin difference is taken.

Usage: emseq_jeepers_pool_controls.py <runs_dir>

<runs_dir> holds <setting>/<control>_draw<N>/t_JEEPERS_Sibling.txt for the
settings 0.9 and default and the controls puc19 and lambda. Prints one
tab-separated row per setting, control, method, and bin. The verdict column
is filled for pUC19 only: lambda has no bar.
"""
from __future__ import annotations

import glob
import os
import sys

import emseq_jeepers_score as score

HEADER = "setting\tcontrol\tdraws\tfragments\tmethod\tbin\tn\tmeth\tpct\tmid_minus_bin\tlower\tupper\tbin_verdict"


def pooled_rows(runs_dir: str):
    """Yield one output line per setting, control, method, and bin."""
    for setting in ("0.9", "default"):
        for control in ("puc19", "lambda"):
            dirs = sorted(glob.glob(os.path.join(runs_dir, setting, f"{control}_draw*")))
            fragments = []
            for d in dirs:
                fragments.extend(score.read_table(os.path.join(d, "t_JEEPERS_Sibling.txt")))
            if not fragments:
                continue
            for r in score.control_rows(fragments):
                pct = 100 * r["meth"] / r["n"] if r["n"] else float("nan")
                if r["interval"] is None:
                    tail = "\t\t\t"
                else:
                    verdict = r["verdict"] if control == "puc19" else ""
                    tail = "%.2f\t%.2f\t%.2f\t%s" % (*r["interval"], verdict)
                yield (f"{setting}\t{control}\t{len(dirs)}\t{len(fragments)}\t{r['method']}\t{r['bin']}"
                       f"\t{r['n']}\t{r['meth']}\t{pct:.2f}\t{tail}")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    print(HEADER)
    for line in pooled_rows(argv[0]):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
