# tests/unit/test_emseq_config_contract.py
import emseq_config_contract as ecc


def test_slice_preamble_stops_at_first_rule_and_drops_configfile():
    text = (
        'configfile: "config/test.yaml"\n'
        "x = config['a']\n"
        "rule all:\n"
        "    input: x\n"
        "y = config['should_not_appear']\n"
    )
    pre = ecc.slice_preamble(text)
    assert "configfile:" not in pre
    assert "config['a']" in pre
    assert "should_not_appear" not in pre


def test_configpath_is_hashable_and_value_equal():
    a = ecc.ConfigPath(path=("envs", "emseq"), requiredness="mandatory")
    b = ecc.ConfigPath(path=("envs", "emseq"), requiredness="mandatory")
    assert a == b
    assert len({a, b}) == 1
