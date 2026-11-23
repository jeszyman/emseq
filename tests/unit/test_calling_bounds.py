"""The pipeline's MethylDackel calls exclude the NEB read-2 bases.

Both calling rules of workflows/emseq.smk (host and spike-in) pass
EMSEQ_CALL_BOUNDS to MethylDackel extract as a rule parameter, and the
constant equals the "neb" setting compared in workflows/validation.smk, whose
effect on constructed read pairs tests/unit/test_methyldackel_bounds.py pins.
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SMK = (REPO / "workflows" / "emseq.smk").read_text()


def rule_body(name: str) -> str:
    match = re.search(rf"^rule {name}:\n(.*?)(?=^rule |\Z)", SMK, re.S | re.M)
    assert match, f"rule {name} not found"
    return match.group(1)


def test_constant_equals_validation_neb_setting():
    match = re.search(r'^EMSEQ_CALL_BOUNDS = "([^"]*)"', SMK, re.M)
    assert match, "EMSEQ_CALL_BOUNDS not defined"
    validation = (REPO / "workflows" / "validation.smk").read_text()
    neb = re.search(r'"neb": "([^"]*)"', validation)
    assert neb
    assert match.group(1) == neb.group(1) == "--nOT 0,0,0,5 --nOB 0,0,5,0"


def test_both_calling_rules_pass_the_bounds():
    for rule in ("emseq_methyldackel", "emseq_methyldackel_spike"):
        body = rule_body(rule)
        assert re.search(r"^\s+bounds\s*=\s*EMSEQ_CALL_BOUNDS,", body, re.M), rule
        extract = body[body.index("MethylDackel extract"):]
        assert "{params.bounds}" in extract.split('"{input.fasta}"')[0], rule
