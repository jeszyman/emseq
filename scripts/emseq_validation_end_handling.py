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
    p.add_argument(
        "--fragment-dir",
        type=Path,
        help="directory of {library}.{source}.ends.tsv files from the fragment-end tally",
    )
    p.add_argument(
        "--fragment-sources",
        nargs="*",
        default=[],
        help="sources tallied by distance from the strand ends",
    )
    p.add_argument(
        "--check-index",
        type=int,
        help="index of the library that checks the measures and is left out of the judged tallies",
    )
    p.add_argument(
        "--extract-dir",
        type=Path,
        help="directory of {library}.{setting}_CpG.bedGraph files from the calling settings",
    )
    p.add_argument(
        "--extract-settings",
        nargs="*",
        default=[],
        help="calling settings compared, with 'none' as the reference",
    )
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
    if args.fragment_sources:
        values.update(
            fragment_values(
                args.fragment_dir,
                args.fragment_sources,
                index,
                args.check_index,
                args.out_dir / "fragment_end.tsv",
            )
        )
    if args.extract_settings:
        values.update(
            extract_values(
                args.extract_dir, args.extract_settings, index, args.check_index
            )
        )
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


# ==============================================================================
# Plasma measures
# ==============================================================================

# Two-sided 95% critical values of Student's t by degrees of freedom
T_975 = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
    11: 2.201,
    12: 2.179,
    13: 2.160,
    14: 2.145,
    15: 2.131,
    20: 2.086,
    30: 2.042,
}
LOSS_DISTANCES = (1, 5, 10, 20, 30, 40, 50, 80)
LONG_MAX_DISTANCE = (
    80  # a base this near one end of a long fragment is nearer it than the other end
)
CHECK_DISTANCES = range(71, 81)


def proportion_difference(a, b):
    """a minus b as proportions of (methylated, unmethylated) counts, with a 95% interval."""
    na, nb = a[0] + a[1], b[0] + b[1]
    if na == 0 or nb == 0:
        return None
    pa, pb = a[0] / na, b[0] / nb
    half = Z * (pa * (1 - pa) / na + pb * (1 - pb) / nb) ** 0.5
    return pa - pb, pa - pb - half, pa - pb + half


def end_loss(five_prime, three_prime):
    """Per distance, 5' methylation minus 3' methylation with its 95% interval."""
    return [proportion_difference(a, b) for a, b in zip(five_prime, three_prime)]


def loss_run(losses, threshold=THRESHOLD):
    """Count of consecutive distances from the first whose whole interval exceeds the threshold."""
    run = 0
    for loss in losses:
        if loss is None or loss[1] <= threshold:
            break
        run += 1
    return run


def pool(counts):
    return sum(m for m, _ in counts), sum(u for _, u in counts)


def mate_difference(read1, read2):
    """Read 1 minus read 2 methylation, pooled over positions, with a 95% interval."""
    return proportion_difference(pool(read1), pool(read2))


def paired_interval(differences):
    """Mean of per-library differences with a 95% t interval."""
    n = len(differences)
    mean = sum(differences) / n
    if n < 2:
        return mean, None, None
    sd = (sum((d - mean) ** 2 for d in differences) / (n - 1)) ** 0.5
    t = T_975.get(n - 1, Z)
    half = t * sd / n**0.5
    return mean, mean - half, mean + half


def read_ends(path):
    """{fragments: {(end, mate, distance): (methylated, unmethylated, kept_neb, kept_all10)}}."""
    tallies = {}
    for line in Path(path).read_text().splitlines()[1:]:
        fragments, end, mate, distance, meth, unmeth, kept_neb, kept_all10 = line.split(
            "\t"
        )
        tallies.setdefault(fragments, {})[(end, int(mate), int(distance))] = (
            int(meth),
            int(unmeth),
            int(kept_neb),
            int(kept_all10),
        )
    return tallies


def summed(libraries, fragments, end, mate, distances):
    """Per distance, counts summed over libraries."""
    rows = []
    for distance in distances:
        cells = [
            t.get(fragments, {}).get((end, mate, distance), (0, 0, 0, 0))
            for t in libraries
        ]
        rows.append(tuple(sum(c[i] for c in cells) for i in range(4)))
    return rows


def fragment_values(fragment_dir, sources, index, check_index, table_path):
    """Loss by distance from the strand ends, judged on the libraries other than the check library."""
    values, table = {}, []
    long_distances = range(1, LONG_MAX_DISTANCE + 1)
    for source in sources:
        tallies = {
            idx: read_ends(fragment_dir / f"{library}.{source}.ends.tsv")
            for library, idx in index.items()
        }
        judged = [t for idx, t in tallies.items() if idx != check_index]
        five = summed(judged, "long", "5p", 0, long_distances)
        three = summed(judged, "long", "3p", 0, long_distances)
        losses = end_loss([c[:2] for c in five], [c[:2] for c in three])
        run = loss_run(losses)
        values[f"{source}_n_libraries_judged"] = len(judged)
        values[f"{source}_loss_run"] = run
        for distance in LOSS_DISTANCES:
            loss = losses[distance - 1]
            values[f"{source}_loss_points_d{distance}"] = (
                None if loss is None else round(100 * loss[0], 1)
            )
            for name, counts in (("p5", five), ("p3", three)):
                meth, unmeth = counts[distance - 1][:2]
                values[f"{source}_{name}_pct_d{distance}"] = (
                    round(100 * meth / (meth + unmeth), 1) if meth + unmeth else None
                )
        per_library = []
        for t in judged:
            loss = proportion_difference(
                summed([t], "long", "5p", 0, [10])[0][:2],
                summed([t], "long", "3p", 0, [10])[0][:2],
            )
            if loss is not None:
                per_library.append(round(100 * loss[0], 1))
        values[f"{source}_loss_points_d10_min"] = (
            min(per_library) if per_library else None
        )
        values[f"{source}_loss_points_d10_max"] = (
            max(per_library) if per_library else None
        )
        if check_index in tallies:
            check = [tallies[check_index]]
            far = proportion_difference(
                pool([c[:2] for c in summed(check, "long", "5p", 0, CHECK_DISTANCES)]),
                pool([c[:2] for c in summed(check, "long", "3p", 0, CHECK_DISTANCES)]),
            )
            values[f"{source}_check_difference_points"] = (
                None if far is None else round(100 * far[0], 1)
            )
            values[f"{source}_check_passes"] = (
                far is not None and abs(far[0]) <= THRESHOLD
            )
        if run:
            within = range(1, run + 1)
            mates = mate_difference(
                [c[:2] for c in summed(judged, "all", "3p", 1, within)],
                [c[:2] for c in summed(judged, "all", "3p", 2, within)],
            )
            if mates is not None:
                values[f"{source}_mate_difference_points"] = round(100 * mates[0], 1)
                values[f"{source}_mate_difference_lower_points"] = round(
                    100 * mates[1], 1
                )
                values[f"{source}_mate_difference_upper_points"] = round(
                    100 * mates[2], 1
                )
                values[f"{source}_P5_holds"] = (
                    -THRESHOLD <= mates[1] and mates[2] <= THRESHOLD
                )
            once = summed(judged, "all", "3p", 0, within)
            total = sum(c[0] + c[1] for c in once)
            values[f"{source}_kept_neb_pct"] = (
                round(100 * sum(c[2] for c in once) / total, 1) if total else None
            )
            values[f"{source}_kept_all10_pct"] = (
                round(100 * sum(c[3] for c in once) / total, 1) if total else None
            )
        pooled = {}
        for t in judged:
            for fragments, cells in t.items():
                for key, cell in cells.items():
                    old = pooled.setdefault((fragments,) + key, (0, 0, 0, 0))
                    pooled[(fragments,) + key] = tuple(a + b for a, b in zip(old, cell))
        for key, cell in sorted(pooled.items()):
            table.append((source,) + key + (0,) + cell)
        for idx, t in sorted(tallies.items()):
            for fragments, cells in t.items():
                for key, cell in sorted(cells.items()):
                    table.append((source, fragments) + key + (idx,) + cell)
    with open(table_path, "w") as out:
        out.write(
            "source\tfragments\tend\tmate\tdistance\tlibrary_index\tn_methylated\tn_unmethylated"
            "\tn_kept_neb\tn_kept_all10\n"
        )
        for row in table:
            out.write("\t".join(str(x) for x in row) + "\n")
    return values


def global_cpg(bedgraphs):
    """Per setting, percent methylated over CpGs called under every setting and over all its calls."""
    calls = {}
    for setting, path in bedgraphs.items():
        cells = {}
        for line in Path(path).read_text().splitlines():
            if line.startswith("track"):
                continue
            contig, start, _, _, meth, unmeth = line.split("\t")
            cells[(contig, start)] = (int(meth), int(unmeth))
        calls[setting] = cells
    shared = set.intersection(*(set(c) for c in calls.values()))
    result = {}
    for setting, cells in calls.items():
        meth_all, total_all = sum(m for m, _ in cells.values()), sum(
            m + u for m, u in cells.values()
        )
        meth_shared = sum(cells[k][0] for k in shared)
        total_shared = sum(cells[k][0] + cells[k][1] for k in shared)
        result[setting] = {
            "shared_pct": (
                round(100 * meth_shared / total_shared, 2) if total_shared else None
            ),
            "all_pct": round(100 * meth_all / total_all, 2) if total_all else None,
            "calls": total_all,
        }
    return result


def extract_values(extract_dir, settings, index, check_index):
    """Global CpG methylation under each calling setting, compared with 'none' per library."""
    per_library = {}
    for library, idx in index.items():
        if idx == check_index:
            continue
        per_library[idx] = global_cpg(
            {s: extract_dir / f"{library}.{s}_CpG.bedGraph" for s in settings}
        )
    values = {}
    for setting in settings:
        shared = [r[setting]["shared_pct"] for r in per_library.values()]
        values[f"host_global_{setting}_pct"] = round(sum(shared) / len(shared), 2)
        every = [r[setting]["all_pct"] for r in per_library.values()]
        values[f"host_global_{setting}_all_calls_pct"] = round(
            sum(every) / len(every), 2
        )
        if setting == "none":
            continue
        mean, lower, upper = paired_interval(
            [r[setting]["all_pct"] - r["none"]["all_pct"] for r in per_library.values()]
        )
        values[f"host_global_all_calls_diff_{setting}_points"] = round(mean, 2)
        values[f"host_global_all_calls_diff_{setting}_lower_points"] = (
            None if lower is None else round(lower, 2)
        )
        values[f"host_global_all_calls_diff_{setting}_upper_points"] = (
            None if upper is None else round(upper, 2)
        )
        differences = [
            r[setting]["shared_pct"] - r["none"]["shared_pct"]
            for r in per_library.values()
        ]
        mean, lower, upper = paired_interval(differences)
        values[f"host_global_diff_{setting}_points"] = round(mean, 2)
        values[f"host_global_diff_{setting}_lower_points"] = (
            None if lower is None else round(lower, 2)
        )
        values[f"host_global_diff_{setting}_upper_points"] = (
            None if upper is None else round(upper, 2)
        )
        dropped = [
            100 * (1 - r[setting]["calls"] / r["none"]["calls"])
            for r in per_library.values()
        ]
        values[f"host_calls_dropped_{setting}_pct"] = round(
            sum(dropped) / len(dropped), 2
        )
    return values


if __name__ == "__main__":
    sys.exit(main())
