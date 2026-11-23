#!/usr/bin/env python3
"""Numbers for the comparison of read-end settings on a two-group cohort.

Reads the four tables written by the settings comparison (global.tsv,
global_diff.tsv, tiles.tsv, tiles_pair.tsv) from one results folder and writes
values.json beside them, for the val macro of the method text.

Settings: none (all read positions counted), neb (the NEB pipeline's read-2
bounds), all10 (10 bases excluded at both ends of both mates).

Usage: emseq_validation_settings_values.py <results folder>
Exit codes: 0 success, 1 runtime error, 2 usage error.
"""
import csv
import json
import sys
from pathlib import Path

SETTINGS = ["none", "neb", "all10"]
MASKS = ["neb", "all10"]


def read_tsv(path: Path) -> list:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by_key(rows: list, key: str) -> dict:
    return {r[key]: r for r in rows}


def interval(values: dict, prefix: str, row: dict, digits: int = 2) -> None:
    """Write estimate, lower and upper of one comparison row under prefix."""
    for name, col in (("", "estimate"), ("_lower", "lower"), ("_upper", "upper")):
        values[f"{prefix}{name}"] = round(float(row[col]), digits)


def settings_values(folder: Path) -> dict:
    glob = read_tsv(folder / "global.tsv")
    diff = by_key(read_tsv(folder / "global_diff.tsv"), "comparison")
    tiles = by_key(read_tsv(folder / "tiles.tsv"), "setting")
    pairs = read_tsv(folder / "tiles_pair.tsv")

    libraries = {r["library_index"]: r["group"] for r in glob}
    values: dict[str, float] = {
        "n_libraries": len(libraries),
        "n_patients": sum(g == "patient" for g in libraries.values()),
        "n_healthy": sum(g == "healthy" for g in libraries.values()),
    }
    for s in SETTINGS:
        pct = [float(r["pct"]) for r in glob if r["setting"] == s]
        if len(pct) != len(libraries):
            raise ValueError(f"setting {s}: {len(pct)} rows for {len(libraries)} libraries")
        values[f"global_mean_{s}_pct"] = round(sum(pct) / len(pct), 2)
        interval(values, f"group_diff_{s}_points", diff[f"patient minus healthy, {s}"])
        values[f"tiles_tested_{s}"] = int(tiles[s]["tiles_tested"])
        values[f"called_q01_d25_{s}"] = int(tiles[s]["called_q01_d25"])
        values[f"called_q05_d10_{s}"] = int(tiles[s]["called_q05_d10"])
        for key in (f"tiles_tested_{s}", f"called_q05_d10_{s}"):
            values[f"{key}_text"] = f"{int(values[key]):,}"
    for m in MASKS:
        interval(values, f"shift_{m}_points", diff[f"{m} minus none"])
        interval(values, f"group_shift_{m}_points", diff[f"patient minus healthy, of {m} minus none"])
    sign_differs = 0
    for r in pairs:
        first, second = r["comparison"].split(" against ")
        if second != "none":
            continue
        th = "q01_d25" if float(r["q"]) == 0.01 else "q05_d10"
        values[f"r_meth_diff_{first}"] = round(float(r["r_meth_diff"]), 3)
        values[f"median_abs_change_{first}_points"] = round(float(r["median_abs_change"]), 2)
        values[f"called_both_{th}_{first}"] = int(r["called_both"])
        values[f"jaccard_{th}_{first}_pct"] = round(100 * float(r["jaccard"]), 1)
        sign_differs += int(r["sign_differs"])
    values["sign_differs_vs_none"] = sign_differs
    return values


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    folder = Path(argv[0])
    try:
        values = settings_values(folder)
    except (OSError, KeyError, ValueError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    (folder / "values.json").write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
