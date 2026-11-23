"""Tests for scripts/emseq_validation_settings_values.py on constructed tables."""
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "emseq_validation_settings_values.py"
spec = importlib.util.spec_from_file_location("settings_values", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def write(path: Path, header: list, rows: list) -> None:
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


@pytest.fixture
def folder(tmp_path):
    glob = []
    for lib, group, base in ((1, "healthy", 78.0), (2, "patient", 74.0), (3, "patient", 76.0)):
        for s, add in (("none", 0.0), ("neb", 1.0), ("all10", 1.5)):
            glob.append([lib, s, 100, 1, 1, base + add, group])
    write(tmp_path / "global.tsv", ["library_index", "setting", "positions", "meth", "unmeth", "pct", "group"], glob)
    diff = [["neb minus none", 3, 1.0, 0.9, 1.1], ["all10 minus none", 3, 1.5, 1.4, 1.6],
            ["all10 minus neb", 3, 0.5, 0.4, 0.6]]
    diff += [[f"patient minus healthy, {s}", 3, -3.0, -4.0, -2.0] for s in ("none", "neb", "all10")]
    diff += [["patient minus healthy, of neb minus none", 3, -0.123, -0.2, -0.05],
             ["patient minus healthy, of all10 minus none", 3, -0.15, -0.25, -0.04]]
    write(tmp_path / "global_diff.tsv", ["comparison", "n", "estimate", "lower", "upper"], diff)
    write(tmp_path / "tiles.tsv", ["setting", "tiles_tested", "called_q01_d25", "called_q05_d10"],
          [["none", 1000, 10, 100], ["neb", 1000, 11, 105], ["all10", 999, 14, 120]])
    header = ["comparison", "q", "min_diff", "tiles_shared", "tiles_only_first", "tiles_only_second", "r_meth_diff",
              "median_abs_change", "max_abs_change", "called_first", "called_second", "called_both", "jaccard",
              "sign_differs"]
    pairs = [["neb against none", 0.01, 25, 1000, 0, 0, 0.9978, 0.08, 20, 11, 10, 9, 0.75, 0],
             ["neb against none", 0.05, 10, 1000, 0, 0, 0.9978, 0.08, 20, 105, 100, 95, 0.8636, 1],
             ["all10 against none", 0.01, 25, 999, 0, 1, 0.988, 0.15, 30, 14, 10, 7, 0.4118, 0],
             ["all10 against none", 0.05, 10, 999, 0, 1, 0.988, 0.15, 30, 120, 100, 90, 0.6923, 2],
             ["all10 against neb", 0.01, 25, 999, 0, 1, 0.99, 0.1, 30, 14, 11, 8, 0.47, 5]]
    write(tmp_path / "tiles_pair.tsv", header, pairs)
    return tmp_path


def test_counts_and_global_means(folder):
    v = mod.settings_values(folder)
    assert (v["n_libraries"], v["n_patients"], v["n_healthy"]) == (3, 2, 1)
    assert v["global_mean_none_pct"] == 76.0
    assert v["global_mean_all10_pct"] == 77.5


def test_intervals_and_group_shift(folder):
    v = mod.settings_values(folder)
    assert (v["shift_neb_points"], v["shift_neb_points_lower"], v["shift_neb_points_upper"]) == (1.0, 0.9, 1.1)
    assert v["group_shift_neb_points"] == -0.12
    assert v["group_diff_all10_points_lower"] == -4.0


def test_tiles_against_none_only(folder):
    v = mod.settings_values(folder)
    assert v["called_q01_d25_all10"] == 14
    assert v["tiles_tested_none_text"] == "1,000"
    assert v["called_both_q05_d10_neb"] == 95
    assert v["jaccard_q01_d25_all10_pct"] == 41.2
    assert v["r_meth_diff_neb"] == 0.998
    # the all10-against-neb row (5 sign changes) is not a comparison against none
    assert v["sign_differs_vs_none"] == 3


def test_missing_setting_row_is_an_error(folder):
    rows = (folder / "global.tsv").read_text().splitlines()
    (folder / "global.tsv").write_text("\n".join(r for r in rows if not r.startswith("1\tneb")) + "\n")
    with pytest.raises(ValueError, match="setting neb"):
        mod.settings_values(folder)


def test_main_writes_values_json(folder):
    assert mod.main([str(folder)]) == 0
    assert json.loads((folder / "values.json").read_text())["tiles_tested_none"] == 1000
    assert mod.main([]) == 2
