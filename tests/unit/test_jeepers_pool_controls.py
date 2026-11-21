"""Pooling of control draws by emseq_jeepers_pool_controls.py, on hand-built rows."""
import emseq_jeepers_pool_controls as pool
from test_jeepers_score import row


def write_draw(base, setting, name, rows):
    d = base / setting / name
    d.mkdir(parents=True)
    (d / "t_JEEPERS_Sibling.txt").write_text("".join(rows))


def test_draws_are_pooled_before_the_difference_is_taken(tmp_path):
    # each draw alone holds 1,000 fragments; the 3' call is unmethylated in 100 of the second draw only
    for n, unmeth in ((1, 0), (2, 100)):
        rows = [
            row(10000001, 10000200, "OT", f"d{n}r{i}", "10000100,10000196", "M," + ("U" if i < unmeth else "M"),
                contig="chrpUC19")
            for i in range(1000)
        ]
        write_draw(tmp_path, "0.9", f"puc19_draw{n}", rows)
    lines = [line.split("\t") for line in pool.pooled_rows(str(tmp_path))]
    by_bin = {(f[4], f[5]): f for f in lines}
    assert {f[0] for f in lines} == {"0.9"} and {f[1] for f in lines} == {"puc19"}
    f = by_bin[("uncorrected", "1-10")]
    assert (f[2], f[3], f[6], f[7]) == ("2", "2000", "2000", "1900")
    assert f[9] == "5.00" and f[12] == "met"
    assert by_bin[("uncorrected", "mid")][6:8] == ["2000", "2000"]


def test_lambda_rows_carry_no_verdict(tmp_path):
    rows = [row(10000001, 10000200, "OT", f"r{i}", "10000100,10000196", "U,U", contig="chrLambda") for i in range(50)]
    write_draw(tmp_path, "default", "lambda_draw1", rows)
    lines = [line.split("\t") for line in pool.pooled_rows(str(tmp_path))]
    assert lines and all(f[1] == "lambda" and f[12] == "" for f in lines if f[5] != "mid")
