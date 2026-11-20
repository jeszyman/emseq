#!/usr/bin/env python3
"""Numbers and position table for the read-end handling validation.

Reads MethylDackel mbias tables and suggested-bounds lines for the lambda and
pUC19 controls of each library, and writes values.json and position.tsv.

A position counts as biased when the whole 95% Wilson interval of its pooled
methylation lies more than 5 percentage points from the pooled methylation of
the middle 60% of positions. The threshold and the reference are the ones
MethylDackel uses for its suggested bounds; the interval requirement keeps a
position with few calls from being called biased by chance.

Exit codes: 0 success, 1 runtime error, 2 usage error.
"""
import argparse
import json
import re
import sys
from pathlib import Path

# ==============================================================================
# Constants
# ==============================================================================

Z = 1.959963984540054  # two-sided 95%
THRESHOLD = 0.05  # MethylDackel's bias threshold
CLIP_WINDOW = 10  # bases a bisulfite-style clip removes at a read end

# Strand whose mbias position 1 sits at each read end (MethylDackel counts
# positions from the left end of the stored alignment); pinned by
# tests/unit/test_methyldackel_bounds.py.
END_STRAND = {"r1_start": "OT", "r1_end": "OB", "r2_start": "OB", "r2_end": "OT"}
END_READ = {"r1_start": 1, "r1_end": 1, "r2_start": 2, "r2_end": 2}
# Index of each end's left bound within "--OT a,b,c,d" / "--OB a,b,c,d"
END_BOUND_FIELD = {"r1_start": 0, "r1_end": 0, "r2_start": 2, "r2_end": 2}

CONTROL_CONTEXTS = {"puc19": ["CpG"], "unmeth_lambda": ["CpG", "CHG", "CHH"]}
OTHER_ENDS = ["r1_start", "r1_end", "r2_end"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--mbias-dir",
        required=True,
        type=Path,
        help="directory of {library}.{source}.{context}.txt and .bounds.txt files",
    )
    p.add_argument(
        "--control-dir",
        required=True,
        type=Path,
        help="directory of {library}.{source}.reads.txt and .tlen.tsv files",
    )
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument(
        "--library-index",
        required=True,
        nargs="+",
        metavar="LIBRARY=INDEX",
        help="index that stands for each library in the output files",
    )
    p.add_argument("--methyldackel-version", required=True)
    p.add_argument("--pipeline-commit", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    index = dict(item.rsplit("=", 1) for item in args.library_index)
    index = {lib: int(i) for lib, i in index.items()}
    try:
        counts, bounds = load(args.mbias_dir, index)
    except (ValueError, FileNotFoundError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_position_table(counts, args.out_dir / "position.tsv")
    values = build_values(
        counts, bounds, len(index), args.methyldackel_version, args.pipeline_commit
    )
    values.update(control_values(args.control_dir, index, values["read_length"]))
    (args.out_dir / "values.json").write_text(
        json.dumps(values, indent=2, sort_keys=True) + "\n"
    )
    return 0


# ==============================================================================
# Statistics
# ==============================================================================


def wilson(k, n):
    """95% Wilson interval for k of n, or None when n is 0."""
    if n == 0:
        return None
    p = k / n
    mid = p + Z * Z / (2 * n)
    half = Z * ((p * (1 - p) / n + Z * Z / (4 * n * n)) ** 0.5)
    denom = 1 + Z * Z / n
    return max(0.0, (mid - half) / denom), min(1.0, (mid + half) / denom)


def reference_middle(counts):
    """Pooled methylation of the middle 60% of positions."""
    n = len(counts)
    middle = counts[int(n * 0.2) : n - int(n * 0.2)]
    meth = sum(m for m, _ in middle)
    total = sum(m + u for m, u in middle)
    return meth / total if total else None


def is_biased(count, reference):
    """True when the whole interval lies more than THRESHOLD from the reference."""
    interval = wilson(count[0], count[0] + count[1])
    if interval is None or reference is None:
        return False
    lower, upper = interval
    return upper < reference - THRESHOLD or lower > reference + THRESHOLD


def biased_run(counts, reference):
    """Count of consecutive biased positions from the start of counts."""
    run = 0
    for count in counts:
        if not is_biased(count, reference):
            break
        run += 1
    return run


def biased_within(counts, reference, window):
    """Count of biased positions among the first `window` positions."""
    return sum(is_biased(count, reference) for count in counts[:window])


def judge(tables):
    """Judge predictions P1 to P3 on pooled tables.

    tables maps (source, context) to {end: [(n_methylated, n_unmethylated), ...]}
    with counts ordered from that read end inward.
    """

    def series(source, context, end):
        counts = tables[(source, context)][end]
        return counts, reference_middle(counts)

    counts, ref = series("puc19", "CpG", "r2_start")
    p1_run = biased_run(counts, ref)

    def window_counts(sources, ends):
        found = {}
        for source in sources:
            for context in CONTROL_CONTEXTS[source]:
                for end in ends:
                    counts, ref = series(source, context, end)
                    found[f"{source}.{context}.{end}"] = biased_within(
                        counts, ref, CLIP_WINDOW
                    )
        return found

    p2 = window_counts(["puc19", "unmeth_lambda"], OTHER_ENDS)
    p3 = window_counts(["unmeth_lambda"], list(END_STRAND))
    return {
        "P1": {"holds": p1_run >= 1, "biased_run": p1_run},
        "P2": {"holds": not any(p2.values()), "biased": p2},
        "P3": {"holds": not any(p3.values()), "biased": p3},
    }


# ==============================================================================
# MethylDackel output
# ==============================================================================


def parse_bounds(line, path):
    """Suggested bounds per strand from MethylDackel's suggestion line."""
    if "Suggested inclusion options" not in line:
        raise ValueError(f"no suggested bounds from MethylDackel in {path}")
    bounds = {}
    for strand in ("OT", "OB"):
        m = re.search(rf"--{strand} (\d+),(\d+),(\d+),(\d+)", line)
        bounds[strand] = tuple(int(x) for x in m.groups()) if m else None
    return bounds


def end_bases(bounds, end_strand):
    """Left bound MethylDackel suggests at each read end, or None without one."""
    result = {}
    for end, strand in end_strand.items():
        strand_bounds = bounds.get(strand)
        result[end] = (
            None if strand_bounds is None else strand_bounds[END_BOUND_FIELD[end]]
        )
    return result


def read_mbias(path):
    """{(strand, read): {position: (n_methylated, n_unmethylated)}} from one table."""
    table = {}
    lines = Path(path).read_text().splitlines()
    for line in lines[1:]:
        strand, read, position, meth, unmeth = line.split("\t")
        table.setdefault((strand, int(read)), {})[int(position)] = (
            int(meth),
            int(unmeth),
        )
    return table


def tlen_summary(paths, read_length):
    """Median fragment length and share at or under the read length, pooled over files."""
    lengths = {}
    for path in paths:
        for line in Path(path).read_text().splitlines():
            length, n = line.split("\t")
            lengths[int(length)] = lengths.get(int(length), 0) + int(n)
    pairs = sum(lengths.values())
    if pairs == 0:
        return {"median": None, "share_le_read_pct": None, "pairs": 0}
    ordered = sorted(lengths)
    cumulative, lower, upper = 0, None, None
    for length in ordered:
        cumulative += lengths[length]
        if lower is None and cumulative >= pairs / 2:
            lower = length
        if upper is None and cumulative >= pairs / 2 + 1:
            upper = length
    if pairs % 2 == 1 or upper is None:
        upper = lower
    short = sum(n for length, n in lengths.items() if length <= read_length)
    return {
        "median": (lower + upper) / 2,
        "share_le_read_pct": round(100 * short / pairs, 2),
        "pairs": pairs,
    }


def control_values(control_dir, index, read_length):
    """Usable read counts and fragment lengths of the controls, pooled over libraries."""
    values = {}
    for source in CONTROL_CONTEXTS:
        reads = sum(
            int((control_dir / f"{library}.{source}.reads.txt").read_text().split()[0])
            for library in index
        )
        summary = tlen_summary(
            [control_dir / f"{library}.{source}.tlen.tsv" for library in index],
            read_length,
        )
        values[f"{source}_usable_reads"] = reads
        values[f"{source}_usable_reads_text"] = f"{reads:,}"
        values[f"{source}_tlen_median"] = summary["median"]
        values[f"{source}_share_tlen_le_read_pct"] = summary["share_le_read_pct"]
    return values


def load(mbias_dir, index):
    """Counts and suggested bounds per library, source, and context."""
    counts, bounds = {}, {}
    for library, idx in index.items():
        for source, contexts in CONTROL_CONTEXTS.items():
            for context in contexts:
                stem = mbias_dir / f"{library}.{source}.{context}"
                counts[(source, context, idx)] = read_mbias(f"{stem}.txt")
                bounds_path = Path(f"{stem}.bounds.txt")
                bounds[(source, context, idx)] = parse_bounds(
                    bounds_path.read_text(), bounds_path.name
                )
    return counts, bounds


# ==============================================================================
# Outputs
# ==============================================================================


def pooled_series(counts, source, context, strand, read):
    """Counts per position, pooled over libraries, as a list from position 1."""
    per_library = [
        t.get((strand, read), {})
        for (s, c, _), t in counts.items()
        if (s, c) == (source, context)
    ]
    last = max((max(t) for t in per_library if t), default=0)
    series = []
    for position in range(1, last + 1):
        meth = sum(t.get(position, (0, 0))[0] for t in per_library)
        unmeth = sum(t.get(position, (0, 0))[1] for t in per_library)
        series.append((meth, unmeth))
    return series


def write_position_table(counts, path):
    """Per-library and pooled (library_index 0) counts with Wilson intervals."""
    header = [
        "source",
        "context",
        "strand",
        "read",
        "position",
        "library_index",
        "n_methylated",
        "n_unmethylated",
        "lower",
        "upper",
    ]
    rows = []
    for source, contexts in CONTROL_CONTEXTS.items():
        for context in contexts:
            for strand in ("OT", "OB"):
                for read in (1, 2):
                    pooled = pooled_series(counts, source, context, strand, read)
                    series = {0: pooled}
                    for (s, c, idx), table in counts.items():
                        if (s, c) == (source, context):
                            positions = table.get((strand, read), {})
                            series[idx] = [
                                positions.get(p, (0, 0))
                                for p in range(1, len(pooled) + 1)
                            ]
                    for idx in sorted(series):
                        for position, (meth, unmeth) in enumerate(series[idx], start=1):
                            interval = wilson(meth, meth + unmeth)
                            lower, upper = (
                                ("", "")
                                if interval is None
                                else (f"{interval[0]:.5f}", f"{interval[1]:.5f}")
                            )
                            rows.append(
                                [
                                    source,
                                    context,
                                    strand,
                                    read,
                                    position,
                                    idx,
                                    meth,
                                    unmeth,
                                    lower,
                                    upper,
                                ]
                            )
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")


def build_values(counts, bounds, n_libraries, methyldackel_version, pipeline_commit):
    """Flat dict of the numbers the report quotes."""
    tables = {}
    for source, contexts in CONTROL_CONTEXTS.items():
        for context in contexts:
            tables[(source, context)] = {
                end: pooled_series(
                    counts, source, context, END_STRAND[end], END_READ[end]
                )
                for end in END_STRAND
            }
    judged = judge(tables)
    values = {
        "methyldackel_version": methyldackel_version,
        "pipeline_commit": pipeline_commit,
        "n_libraries": n_libraries,
        "threshold_points": THRESHOLD * 100,
        "clip_window": CLIP_WINDOW,
        "read_length": max(
            len(series) for ends in tables.values() for series in ends.values()
        ),
        "P1_holds": judged["P1"]["holds"],
        "P2_holds": judged["P2"]["holds"],
        "P3_holds": judged["P3"]["holds"],
    }
    for (source, context), ends in tables.items():
        label = source if source == "puc19" else f"{source}_{context}"
        for end, series in ends.items():
            reference = reference_middle(series)
            key = f"{label}_{end}"
            values[f"{key}_biased_run"] = biased_run(series, reference)
            values[f"{key}_biased_in_window"] = biased_within(
                series, reference, CLIP_WINDOW
            )
            values[f"{key}_reference_pct"] = (
                None if reference is None else round(reference * 100, 2)
            )
            first = series[0] if series else (0, 0)
            total = first[0] + first[1]
            values[f"{key}_position1_pct"] = (
                round(first[0] / total * 100, 2) if total else None
            )
            values[f"{key}_position1_calls"] = total
            second = series[1] if len(series) > 1 else (0, 0)
            total2 = second[0] + second[1]
            values[f"{key}_position2_pct"] = (
                round(second[0] / total2 * 100, 2) if total2 else None
            )
            # mean over the first 30 positions that follow the biased run
            after = series[values[f"{key}_biased_run"] : 30]
            after_total = sum(m + u for m, u in after)
            values[f"{key}_mean_after_run_pct"] = (
                round(100 * sum(m for m, _ in after) / after_total, 2)
                if after_total
                else None
            )
            suggested = [
                end_bases(b, END_STRAND)[end]
                for (s, c, _), b in bounds.items()
                if (s, c) == (source, context)
            ]
            values[f"{key}_n_libraries_bound"] = sum(1 for x in suggested if x)
            values[f"{key}_suggested_bound_max"] = max(
                [x for x in suggested if x is not None], default=None
            )
    # totals per control, and the largest lambda value near a read end
    for source, contexts in CONTROL_CONTEXTS.items():
        meth = total = 0
        largest = 0.0
        for context in contexts:
            for strand, read in (("OT", 1), ("OT", 2), ("OB", 1), ("OB", 2)):
                series = pooled_series(counts, source, context, strand, read)
                meth += sum(m for m, _ in series)
                total += sum(m + u for m, u in series)
            for end, series in tables[(source, context)].items():
                for m, u in series[:30]:
                    if m + u:
                        largest = max(largest, 100 * m / (m + u))
        values[f"{source}_total_calls"] = total
        values[f"{source}_total_calls_text"] = f"{total:,}"
        values[f"{source}_overall_pct"] = (
            round(100 * meth / total, 2) if total else None
        )
        values[f"{source}_max_position_pct"] = round(largest, 2)
    return values


if __name__ == "__main__":
    sys.exit(main())
