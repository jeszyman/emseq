#!/usr/bin/env python3
"""Score JEEPERS's sibling-mode output against five bars fixed before the run.

Reads the run folders written by emseq_jeepers_run.sh, laid out as
<runs>/<setting>/<run>/, where setting is "0.9" or "default" and run is
plasma_all, plasma_without_a, plasma_without_b, puc19_draw<N>, or
lambda_draw<N>, and the pairs table of emseq_jeepers_holdout.py. Writes
gap.tsv, opposite_strand.tsv, controls.tsv, and values.json.

Table layout. <out>_Sibling.txt has 35 tab-separated columns and no header:
1 contig, 2 and 3 fragment start and end (1-based, inclusive), 4 block,
5 strand (OT or OB), 6 read name, 11 CpG positions (1-based position of the
C), 12 calls before correction, 25 calls after, 26 CpGs marked for
correction or "-", 27 methylated:total reference counts for the marked CpGs,
35 correction type. The merge table has the same first 24 columns.

Calls. M and m count as methylated, U and u as unmethylated; q and N are
left out of numerator and denominator.

Strand ends. On an OT fragment the 3' end is the fragment end. On an OB
fragment the called base is the G at position + 1 and the 3' end is the
fragment start.

Methods scored on the same fragments:
  uncorrected  calls before correction
  jeepers      calls after correction
  majority     calls before correction, with each CpG JEEPERS marked replaced
               by the call most other fragments hold at that CpG outside
               their own marked CpGs; a tie keeps the fragment's own call
  excl30/52    calls before correction more than 30 (52) bases from the 3' end
  drop_marked  calls before correction without the CpGs JEEPERS marked

Bars (each bin judged on a 95% interval, Wilson for a share and Newcombe for
a difference; a bin whose interval reaches more than 5 points to either side
of its estimate, or that has no calls, is unresolved; a bar is met when all
its bins are met, failed when any bin fails, unresolved otherwise):
  1  rows after correction equal rows before, block by block, and no
     Traceback or ERROR in the logs
  2  plasma, fragments of 160 bp or longer, all fragments: 5' minus 3'
     methylation per bin has an interval that contains 5 points or lies
     below; unresolved for JEEPERS if the uncorrected calls already meet it
  3  pUC19, each draw: mid-fragment minus bin (1-10, 11-20, 21-30 from the
     3' end) contains 5 points or lies below
  4  agreement with the opposite strand, JEEPERS minus uncorrected: above 0
     in bins 1-10 and 11-20, contains 0 or above in 21-30 and 31-50
  5  agreement with the opposite strand, JEEPERS minus the locus majority
     call over marked CpGs, and JEEPERS on all calls within 80 bases of the
     3' end minus each exclusion on the calls it keeps there: contains 0 or
     lies above
"""
from __future__ import annotations

import argparse
import collections
import dataclasses
import glob
import json
import os
import re

Z = 1.959963984540054
BINS = [(1, 10), (11, 20), (21, 30), (31, 50), (51, 80)]
MIN_GAP_LENGTH = 160
MID_MARGIN = 30
RESOLUTION = 5.0  # points an interval may reach to either side of its estimate
SETTINGS = ("0.9", "default")
EXCLUSIONS = {"excl30": 30, "excl52": 52}
METHODS = ("uncorrected", "jeepers", "majority", "excl30", "excl52", "drop_marked")


@dataclasses.dataclass
class Fragment:
    contig: str
    start: int
    end: int
    strand: str
    name: str
    pos: list
    before: list
    after: list
    marked: list
    support: dict
    kind: str
    block: str = ""


def _split(field: str) -> list:
    return [] if field in ("-", "") else field.split(",")


def parse_row(line: str) -> Fragment:
    """A Fragment from one row of a sibling table (35 columns) or a merge table (24)."""
    f = line.rstrip("\n").split("\t")
    pos = [int(p) for p in _split(f[10])]
    before = _split(f[11])
    if len(f) < 35:
        return Fragment(f[0], int(f[1]), int(f[2]), f[4], f[5], pos, before, list(before), [], {}, "", f[3])
    marked = [int(p) for p in _split(f[25])]
    counts = _split(f[26])
    support: dict = {}
    for i, p in enumerate(marked):
        entry = counts[i] if i < len(counts) else "-"
        support[p] = int(entry.split(":")[1]) if ":" in entry and entry.split(":")[1].isdigit() else None
    return Fragment(f[0], int(f[1]), int(f[2]), f[4], f[5], pos, before, _split(f[24]), marked, support, f[34], f[3])


def read_table(path: str, blocks=None):
    """Fragments of a sibling or merge table; with blocks, only those of the block numbers given."""
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("Chr\t") and "\tBlock\t" not in line[:60]:
                frag = parse_row(line)
                if blocks is None or float(frag.block) in blocks:
                    yield frag


def complete_blocks(run_dir: str):
    """(block numbers whose rows after correction equal the rows before, all block numbers).

    A block that JEEPERS did not finish holds only the fragments it wrote
    before stopping, which are those it left unchanged, so scoring is kept to
    the blocks it returned whole.
    """
    done, every = set(), set()
    with open(os.path.join(run_dir, "blocks.tsv")) as fh:
        next(fh)
        for line in fh:
            block, merge_rows, sibling_rows = line.split("\t")
            every.add(float(block))
            if int(merge_rows) == int(sibling_rows):
                done.add(float(block))
    return done, every


def distances(frag: Fragment, pos: int) -> tuple[int, int]:
    """(bases from the 3' end, bases from the 5' end) of the called base, counting the end base as 1."""
    if frag.strand == "OT":
        return frag.end - pos + 1, pos - frag.start + 1
    return pos + 2 - frag.start, frag.end - pos


def bin_of(d: int):
    for lo, hi in BINS:
        if lo <= d <= hi:
            return f"{lo}-{hi}"
    return None


def binary(symbol: str):
    if symbol in ("M", "m"):
        return 1
    if symbol in ("U", "u"):
        return 0
    return None


# ---- comparators ----------------------------------------------------------


def tally_unmarked(fragments, tally=None) -> dict:
    """{position: [unmethylated, methylated]} over before-correction calls at unmarked CpGs."""
    tally = collections.defaultdict(lambda: [0, 0]) if tally is None else tally
    for frag in fragments:
        marked = set(frag.marked)
        for p, s in zip(frag.pos, frag.before):
            b = binary(s)
            if b is not None and p not in marked:
                tally[p][b] += 1
    return tally


def majority_for(frag: Fragment, tally) -> dict:
    """{marked position: call} for one fragment: M or U by majority, own call on a tie or no other fragment."""
    own = dict(zip(frag.pos, frag.before))
    out = {}
    for p in frag.marked:
        u, m = tally.get(p, (0, 0))
        out[p] = "M" if m > u else "U" if u > m else own.get(p, "N")
    return out


def majority_calls(fragments) -> dict:
    """{read name: {marked position: call}} for each fragment with marked CpGs."""
    fragments = list(fragments)
    tally = tally_unmarked(fragments)
    return {f.name: majority_for(f, tally) for f in fragments if f.marked}


def drop_marked(frag: Fragment) -> Fragment:
    """The fragment without the calls at its marked CpGs."""
    marked = set(frag.marked)
    keep = [i for i, p in enumerate(frag.pos) if p not in marked]
    return dataclasses.replace(
        frag,
        pos=[frag.pos[i] for i in keep],
        before=[frag.before[i] for i in keep],
        after=[frag.after[i] for i in keep if i < len(frag.after)],
        marked=[],
        support={},
    )


def method_calls(frag: Fragment, method: str, tally=None) -> dict:
    """{position: 0 or 1} of the fragment's calls under one method; q and N are left out."""
    if method == "jeepers":
        symbols = dict(zip(frag.pos, frag.after))
    elif method == "drop_marked":
        dropped = drop_marked(frag)
        symbols = dict(zip(dropped.pos, dropped.before))
    else:
        symbols = dict(zip(frag.pos, frag.before))
        if method == "majority" and frag.marked:
            symbols.update(majority_for(frag, tally or {}))
        elif method in EXCLUSIONS:
            symbols = {p: s for p, s in symbols.items() if distances(frag, p)[0] > EXCLUSIONS[method]}
    return {p: b for p, b in ((p, binary(s)) for p, s in symbols.items()) if b is not None}


# ---- intervals and verdicts -----------------------------------------------


def _wilson_fraction(k, n):
    p = k / n
    mid = p + Z * Z / (2 * n)
    half = Z * ((p * (1 - p) / n + Z * Z / (4 * n * n)) ** 0.5)
    denom = 1 + Z * Z / n
    return p, max(0.0, (mid - half) / denom), min(1.0, (mid + half) / denom)


def wilson(k, n):
    """(estimate, lower, upper) in points for k of n by the Wilson score interval; None when n is 0."""
    if n == 0:
        return None
    return tuple(100 * x for x in _wilson_fraction(k, n))


def newcombe(k1, n1, k2, n2):
    """(estimate, lower, upper) in points for k1/n1 minus k2/n2 by Newcombe's hybrid score interval."""
    if n1 == 0 or n2 == 0:
        return None
    p1, l1, u1 = _wilson_fraction(k1, n1)
    p2, l2, u2 = _wilson_fraction(k2, n2)
    d = p1 - p2
    lo = d - ((p1 - l1) ** 2 + (u2 - p2) ** 2) ** 0.5
    hi = d + ((u1 - p1) ** 2 + (p2 - l2) ** 2) ** 0.5
    return 100 * d, 100 * lo, 100 * hi


def judge(interval, limit: float, mode: str) -> str:
    """met, failed, or unresolved for one interval.

    at_most: the interval contains the limit or lies below it.
    above: the interval lies above the limit.
    at_least: the interval contains the limit or lies above it.
    """
    if interval is None:
        return "unresolved"
    est, lo, hi = interval
    if hi - est > RESOLUTION or est - lo > RESOLUTION:
        return "unresolved"
    if mode == "at_most":
        ok = lo <= limit
    elif mode == "above":
        ok = lo > limit
    elif mode == "at_least":
        ok = hi >= limit
    else:
        raise ValueError(mode)
    return "met" if ok else "failed"


def bar_verdict(verdicts) -> str:
    verdicts = list(verdicts)
    if any(v == "failed" for v in verdicts):
        return "failed"
    if verdicts and all(v == "met" for v in verdicts):
        return "met"
    return "unresolved"


# ---- tallies --------------------------------------------------------------


def new_acc():
    return collections.defaultdict(lambda: [0, 0])


def add_gap(acc, frag: Fragment, tally=None, methods=("uncorrected", "jeepers", "majority", "drop_marked")):
    """Add one fragment's calls to acc[(method, "3p" or "5p", bin)] = [calls, methylated].

    Only fragments of MIN_GAP_LENGTH or longer count, so that a call within 80
    bases of one end is not within 80 of the other.
    """
    if frag.end - frag.start + 1 < MIN_GAP_LENGTH:
        return
    for method in methods:
        labels = [method] + ([f"{method}_sf"] if frag.kind == "Sf" and method in ("uncorrected", "jeepers") else [])
        for p, b in method_calls(frag, method, tally).items():
            d3, d5 = distances(frag, p)
            for end, d in (("3p", d3), ("5p", d5)):
                label = bin_of(d)
                if label:
                    for name in labels:
                        cell = acc[(name, end, label)]
                        cell[0] += 1
                        cell[1] += b


def gap_rows(acc, methods):
    """One dict per method and bin with the 5' minus 3' interval and its bin verdict."""
    rows = []
    for method in methods:
        for lo, hi in BINS:
            label = f"{lo}-{hi}"
            n3, m3 = acc.get((method, "3p", label), (0, 0))
            n5, m5 = acc.get((method, "5p", label), (0, 0))
            diff = newcombe(m5, n5, m3, n3)
            rows.append(
                dict(method=method, bin=label, n_3p=n3, meth_3p=m3, n_5p=n5, meth_5p=m5,
                     interval=diff, verdict=judge(diff, 5.0, "at_most"))
            )
    return rows


def mid_value(fragments):
    """(methylated, calls) per method over calls more than MID_MARGIN bases from both strand ends."""
    out = {"uncorrected": [0, 0], "jeepers": [0, 0]}
    for frag in fragments:
        for method in out:
            for p, b in method_calls(frag, method).items():
                d3, d5 = distances(frag, p)
                if d3 > MID_MARGIN and d5 > MID_MARGIN:
                    out[method][0] += b
                    out[method][1] += 1
    return {k: tuple(v) for k, v in out.items()}


def control_rows(fragments):
    """Per method: the mid-fragment value and, for the three bins nearest the 3' end, mid minus bin."""
    fragments = list(fragments)
    mid = mid_value(fragments)
    acc = new_acc()
    for frag in fragments:
        for method in mid:
            for p, b in method_calls(frag, method).items():
                label = bin_of(distances(frag, p)[0])
                if label:
                    acc[(method, label)][0] += 1
                    acc[(method, label)][1] += b
    rows = []
    for method, (mid_m, mid_n) in mid.items():
        rows.append(dict(method=method, bin="mid", n=mid_n, meth=mid_m, interval=None, verdict=""))
        for lo, hi in BINS[:3]:
            label = f"{lo}-{hi}"
            n, m = acc.get((method, label), (0, 0))
            diff = newcombe(mid_m, mid_n, m, n)
            rows.append(dict(method=method, bin=label, n=n, meth=m, interval=diff, verdict=judge(diff, 5.0, "at_most")))
    return rows


def support_class(n):
    if n is None:
        return "fallback"
    return "1" if n == 1 else "2" if n == 2 else "3+" if n >= 3 else "0"


def score_strand(acc, frag: Fragment, partner: dict, tally=None):
    """Add one strand's agreement with its partner strand to acc.

    partner is {position: 0 or 1}, the partner strand's calls before
    correction. A call counts only when it is methylated or unmethylated on
    both strands. Keys of acc:
      (method, scope) = [calls, calls equal to the partner's], with scope a
        bin from the 3' end, "1-80", "all", or "marked" (the CpGs JEEPERS
        marked on this strand)
      ("changed_to_methylated", support class) = [calls, partner unmethylated]
      ("dropped_marked", bin) = [calls, partner unmethylated]
    """
    marked = set(frag.marked)
    for method in METHODS:
        labels = [method] + ([f"{method}_sf"] if frag.kind == "Sf" and method in ("uncorrected", "jeepers") else [])
        for p, b in method_calls(frag, method, tally).items():
            pb = partner.get(p)
            if pb is None:
                continue
            d3 = distances(frag, p)[0]
            scopes = ["all"]
            label = bin_of(d3)
            if label:
                scopes += [label, "1-80"]
            if p in marked:
                scopes.append("marked")
            for name in labels:
                for scope in scopes:
                    cell = acc[(name, scope)]
                    cell[0] += 1
                    cell[1] += int(b == pb)
    before = {p: binary(s) for p, s in zip(frag.pos, frag.before)}
    after = {p: binary(s) for p, s in zip(frag.pos, frag.after)}
    for p in frag.pos:
        pb = partner.get(p)
        if pb is None:
            continue
        if before.get(p) == 0 and after.get(p) == 1:
            cell = acc[("changed_to_methylated", support_class(frag.support.get(p)))]
            cell[0] += 1
            cell[1] += int(pb == 0)
        if p in marked and before.get(p) is not None:
            d3 = distances(frag, p)[0]
            for scope in filter(None, ("all", bin_of(d3), "1-80" if d3 <= 80 else None)):
                cell = acc[("dropped_marked", scope)]
                cell[0] += 1
                cell[1] += int(pb == 0)


def difference_row(acc, bar, scope, method, comparator, mode, comparator_scope=None):
    n1, k1 = acc.get((method, scope), (0, 0))
    n2, k2 = acc.get((comparator, comparator_scope or scope), (0, 0))
    diff = newcombe(k1, n1, k2, n2)
    return dict(bar=bar, scope=scope, method=method, comparator=comparator, n=n1, agree=k1,
                comparator_n=n2, comparator_agree=k2, interval=diff, verdict=judge(diff, 0.0, mode))


def bar4_rows(acc):
    modes = {"1-10": "above", "11-20": "above", "21-30": "at_least", "31-50": "at_least"}
    return [difference_row(acc, 4, scope, "jeepers", "uncorrected", mode) for scope, mode in modes.items()]


def bar5_rows(acc):
    return [
        difference_row(acc, 5, "marked", "jeepers", "majority", "at_least"),
        difference_row(acc, 5, "1-80", "jeepers", "excl30", "at_least"),
        difference_row(acc, 5, "1-80", "jeepers", "excl52", "at_least"),
    ]


# ---- run folders ----------------------------------------------------------


def completeness(run_dir: str) -> dict:
    """Bar 1 for one run folder: block row counts, error text in logs, the runner's exit code."""
    out = dict(run=os.path.basename(run_dir), blocks=0, blocks_differing=None, merge_rows=0, sibling_rows=0,
               log_errors=[], exit_code=None, seconds={})
    blocks = os.path.join(run_dir, "blocks.tsv")
    if os.path.exists(blocks):
        out["blocks_differing"] = 0
        with open(blocks) as fh:
            next(fh)
            for line in fh:
                _, a, b = line.split("\t")
                out["blocks"] += 1
                out["merge_rows"] += int(a)
                out["sibling_rows"] += int(b)
                out["blocks_differing"] += int(int(a) != int(b))
    for log in sorted(glob.glob(os.path.join(run_dir, "*.log"))):
        with open(log, errors="replace") as fh:
            if re.search("Traceback|ERROR", fh.read()):
                out["log_errors"].append(os.path.basename(log))
    code = os.path.join(run_dir, "exit_code")
    if os.path.exists(code):
        out["exit_code"] = int(open(code).read().strip())
    run_tsv = os.path.join(run_dir, "run.tsv")
    if os.path.exists(run_tsv):
        with open(run_tsv) as fh:
            next(fh)
            for line in fh:
                step, seconds, _ = line.split("\t")
                out["seconds"][step] = int(seconds)
    ok = out["blocks_differing"] == 0 and out["blocks"] > 0 and not out["log_errors"] and out["exit_code"] in (0, None)
    out["verdict"] = "met" if ok else "failed"
    return out


def read_pairs(path: str):
    same, near = [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            kind, top, bottom, removed_in = line.rstrip("\n").split("\t")
            (same if kind == "same" else near).append((top, bottom, removed_in))
    return same, near


def scan(path: str, names=frozenset(), gap_acc=None, blocks=None):
    """One pass over a sibling table.

    Returns the unmarked-call tally, the fragments named in `names`, and
    counts of correction types and call symbols. With gap_acc, unmarked
    fragments are added to it at once and marked fragments are returned to be
    added once the tally is complete (their majority calls need it).
    """
    tally = collections.defaultdict(lambda: [0, 0])
    kept, deferred = {}, []
    counts = dict(kinds=collections.Counter(), after_symbols=collections.Counter(),
                  before_symbols=collections.Counter(), marked_fragments=0, marked_cpgs=0, fragments=0,
                  call_count_differs=0)
    for frag in read_table(path, blocks):
        tally_unmarked([frag], tally)
        counts["fragments"] += 1
        counts["kinds"][frag.kind] += 1
        counts["after_symbols"].update(frag.after)
        counts["before_symbols"].update(frag.before)
        # calls after correction are matched to positions by order, so a row with
        # a different number of calls than positions would be misread
        counts["call_count_differs"] += int(not len(frag.pos) == len(frag.before) == len(frag.after))
        if frag.marked:
            counts["marked_fragments"] += 1
            counts["marked_cpgs"] += len(frag.marked)
        if frag.name in names:
            kept[frag.name] = frag
        if gap_acc is not None:
            if frag.marked:
                deferred.append(frag)
            else:
                add_gap(gap_acc, frag)
    if gap_acc is not None:
        for frag in deferred:
            add_gap(gap_acc, frag, tally)
    return tally, kept, counts


def interval_cols(interval):
    if interval is None:
        return ["", "", ""]
    return [f"{x:.2f}" for x in interval]


def pct(k, n):
    return "" if n == 0 else f"{100 * k / n:.2f}"


def score_setting(runs: str, setting: str, same, near):
    base = os.path.join(runs, setting)
    out = dict(setting=setting, partial=[])
    names = {n for top, bottom, _ in same + near for n in (top, bottom)}

    # bar 1
    run_dirs = sorted(d for d in glob.glob(os.path.join(base, "*")) if os.path.isdir(d))
    out["completeness"] = [completeness(d) for d in run_dirs]
    out["bar1"] = bar_verdict(c["verdict"] for c in out["completeness"])

    # partner calls: before-correction calls from the table of all fragments
    partner = {}
    merge = os.path.join(base, "plasma_all", "t_Merge.txt")
    if os.path.exists(merge):
        for frag in read_table(merge):
            if frag.name in names:
                partner[frag.name] = method_calls(frag, "uncorrected")

    # bar 2 and the one-base-off pairs, from the run on all fragments
    gap_acc, near_acc = new_acc(), new_acc()
    sib = os.path.join(base, "plasma_all", "t_JEEPERS_Sibling.txt")
    if os.path.exists(sib):
        done, every = complete_blocks(os.path.dirname(sib))
        if done != every:
            out["partial"].append(f"plasma_all scored on blocks {sorted(done)} of {sorted(every)}")
        tally, kept, counts = scan(sib, names, gap_acc, done)
        out["plasma_all_counts"] = {k: (dict(v) if isinstance(v, collections.Counter) else v) for k, v in counts.items()}
        for top, bottom, _ in near:
            for own, other in ((top, bottom), (bottom, top)):
                if own in kept and other in partner:
                    score_strand(near_acc, kept[own], partner[other], tally)
    else:
        out["partial"].append("plasma_all missing")
    gap_methods = ("uncorrected", "jeepers", "majority", "drop_marked", "uncorrected_sf", "jeepers_sf")
    out["gap"] = gap_rows(gap_acc, gap_methods)
    uncorrected = bar_verdict(r["verdict"] for r in out["gap"] if r["method"] == "uncorrected")
    jeepers = bar_verdict(r["verdict"] for r in out["gap"] if r["method"] == "jeepers")
    out["bar2_uncorrected"] = uncorrected
    out["bar2_jeepers_as_scored"] = jeepers
    out["bar2"] = "unresolved" if uncorrected == "met" else jeepers

    # bars 4 and 5: each strand from the run in which its partner was removed
    same_acc = new_acc()
    scored = 0
    for removed in ("A", "B"):
        sib = os.path.join(base, f"plasma_without_{removed.lower()}", "t_JEEPERS_Sibling.txt")
        if not os.path.exists(sib):
            out["partial"].append(f"plasma_without_{removed.lower()} missing")
            continue
        done, every = complete_blocks(os.path.dirname(sib))
        if done != every:
            out["partial"].append(f"plasma_without_{removed.lower()} scored on blocks {sorted(done)} of {sorted(every)}")
        tally, kept, _ = scan(sib, names, blocks=done)
        for top, bottom, top_removed_in in same:
            own, other = (bottom, top) if top_removed_in == removed else (top, bottom)
            if own in kept and other in partner:
                score_strand(same_acc, kept[own], partner[other], tally)
                scored += 1
    out["strands_scored"] = scored
    out["same_acc"] = {"|".join(k): v for k, v in same_acc.items()}
    out["near_acc"] = {"|".join(k): v for k, v in near_acc.items()}
    out["bar4_rows"] = bar4_rows(same_acc)
    out["bar5_rows"] = bar5_rows(same_acc)
    out["bar4"] = bar_verdict(r["verdict"] for r in out["bar4_rows"])
    out["bar5"] = bar_verdict(r["verdict"] for r in out["bar5_rows"])
    out["near_bar4_rows"] = bar4_rows(near_acc)
    out["near_bar5_rows"] = bar5_rows(near_acc)

    # bar 3 and lambda
    out["controls"] = []
    for control in ("puc19", "lambda"):
        for d in sorted(glob.glob(os.path.join(base, f"{control}_*"))):
            sib = os.path.join(d, "t_JEEPERS_Sibling.txt")
            if not os.path.exists(sib):
                out["partial"].append(f"{os.path.basename(d)} missing")
                continue
            fragments = list(read_table(sib))
            kinds = collections.Counter(f.kind for f in fragments)
            for row in control_rows(fragments):
                out["controls"].append(dict(control=control, draw=os.path.basename(d), fragments=len(fragments),
                                            kinds=dict(kinds), **row))
    puc = [r["verdict"] for r in out["controls"] if r["control"] == "puc19" and r["method"] == "jeepers" and r["bin"] != "mid"]
    out["bar3"] = bar_verdict(puc)
    out["bar3_uncorrected"] = bar_verdict(
        r["verdict"] for r in out["controls"] if r["control"] == "puc19" and r["method"] == "uncorrected" and r["bin"] != "mid"
    )
    out["verdicts"] = {f"bar{i}": out[f"bar{i}"] for i in range(1, 6)}
    out["suitable"] = all(v == "met" for v in out["verdicts"].values())
    return out


def write_tables(results, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "gap.tsv"), "w") as fh:
        fh.write("setting\tmethod\tbin\tn_3p\tmeth_3p\tpct_3p\tn_5p\tmeth_5p\tpct_5p\t"
                 "diff_5p_minus_3p\tlower\tupper\tbin_verdict\n")
        for res in results:
            for r in res["gap"]:
                fh.write("\t".join(map(str, [
                    res["setting"], r["method"], r["bin"], r["n_3p"], r["meth_3p"], pct(r["meth_3p"], r["n_3p"]),
                    r["n_5p"], r["meth_5p"], pct(r["meth_5p"], r["n_5p"]), *interval_cols(r["interval"]), r["verdict"],
                ])) + "\n")
    with open(os.path.join(out_dir, "opposite_strand.tsv"), "w") as fh:
        fh.write("setting\tpairs\trow\tmethod\tscope\tn\tcount\tpct\tlower\tupper\t"
                 "bar\tcomparator\tcomparator_n\tcomparator_count\tdiff\tdiff_lower\tdiff_upper\tbin_verdict\n")
        for res in results:
            for pairs, acc_key, bars in (("same", "same_acc", ("bar4_rows", "bar5_rows")),
                                         ("near", "near_acc", ("near_bar4_rows", "near_bar5_rows"))):
                for key, (n, k) in sorted(res[acc_key].items()):
                    method, scope = key.split("|")
                    # count is calls equal to the partner's for agreement rows, and calls
                    # the partner reads unmethylated for the two partner_unmethylated rows
                    row = "partner_unmethylated" if method in ("changed_to_methylated", "dropped_marked") else "agreement"
                    fh.write("\t".join(map(str, [res["setting"], pairs, row, method, scope, n, k, pct(k, n),
                                                 *interval_cols(wilson(k, n))[1:], "", "", "", "", "", "", "", ""])) + "\n")
                for bar_rows in bars:
                    for r in res[bar_rows]:
                        fh.write("\t".join(map(str, [
                            res["setting"], pairs, "difference", r["method"], r["scope"], r["n"], r["agree"],
                            pct(r["agree"], r["n"]), "", "", r["bar"] if pairs == "same" else "", r["comparator"],
                            r["comparator_n"], r["comparator_agree"], *interval_cols(r["interval"]),
                            r["verdict"] if pairs == "same" else "",
                        ])) + "\n")
    with open(os.path.join(out_dir, "controls.tsv"), "w") as fh:
        fh.write("setting\tcontrol\tdraw\tfragments\tmethod\tbin\tn\tmeth\tpct\tmid_minus_bin\tlower\tupper\tbin_verdict\n")
        for res in results:
            for r in res["controls"]:
                verdict = r["verdict"] if r["control"] == "puc19" else ""
                fh.write("\t".join(map(str, [
                    res["setting"], r["control"], r["draw"], r["fragments"], r["method"], r["bin"], r["n"], r["meth"],
                    pct(r["meth"], r["n"]), *interval_cols(r["interval"]), verdict,
                ])) + "\n")
    with open(os.path.join(out_dir, "values.json"), "w") as fh:
        json.dump({"bins": BINS, "settings": {r["setting"]: r for r in results}}, fh, indent=1)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--runs", required=True, help="folder holding <setting>/<run>/")
    p.add_argument("--pairs", required=True, help="pairs table from emseq_jeepers_holdout.py")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--settings", nargs="+", default=list(SETTINGS))
    a = p.parse_args(argv)
    same, near = read_pairs(a.pairs)
    results = [score_setting(a.runs, s, same, near) for s in a.settings]
    write_tables(results, a.out_dir)
    for res in results:
        print(res["setting"], json.dumps(res["verdicts"]), "partial:", res["partial"])


if __name__ == "__main__":
    main()
