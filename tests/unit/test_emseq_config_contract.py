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


def test_extract_mandatory_flat_and_nested():
    src = (
        "a = config['main-data-dir']\n"
        "b = config['envs']['emseq']\n"
        "c = config[\"library-ids\"]\n"
    )
    contract = ecc.extract_from_preamble(src, "w.smk")
    paths = {p.path: p for p in contract.paths}
    assert paths[("main-data-dir",)].requiredness == "mandatory"
    assert paths[("envs", "emseq")].requiredness == "mandatory"
    assert ("library-ids",) in paths


def test_extract_optional_and_conditional_get():
    src = (
        "m = config.get('mosdepth-quant-levels', '1,5,10,20')\n"
        "n = config.get('emseq-mincov', 2)\n"
        "f = config.get('fastp', {}).get('extra', '')\n"
    )
    contract = ecc.extract_from_preamble(src, "w.smk")
    paths = {p.path: p for p in contract.paths}
    assert paths[("mosdepth-quant-levels",)].requiredness == "optional"
    assert paths[("mosdepth-quant-levels",)].default == "1,5,10,20"
    assert paths[("emseq-mincov",)].default == 2
    fx = paths[("fastp", "extra")]
    assert fx.requiredness == "conditional"
    assert fx.condition == ("fastp",)
    assert fx.default == ""


def test_extract_aliases_and_baked_in():
    src = (
        "meth_map = config['meth-map']\n"
        "ENV = config['envs']['emseq']\n"
        "emseq_ref_names = ['chr22']\n"
        "emseq_align_methods = ['bwa_meth', 'biscuit']\n"
    )
    contract, aliases = ecc.extract_from_preamble(src, "w.smk", return_aliases=True)
    assert aliases["meth_map"] == ("meth-map",)
    assert aliases["ENV"] == ("envs", "emseq")
    assert contract.baked_in["emseq_ref_names"] == ["chr22"]
    assert contract.baked_in["emseq_align_methods"] == ["bwa_meth", "biscuit"]


def test_dynamic_key_marks_incomplete():
    src = "name = 'x'\nv = config[name]\n"
    contract = ecc.extract_from_preamble(src, "w.smk")
    assert contract.incomplete is True
    assert any("dynamic config key" in w for w in contract.warnings)


def test_extract_comprehension_wildcard_subkey():
    src = "REF = {k: v['input'] for k, v in config['emseq_ref_assemblies'].items()}\n"
    contract = ecc.extract_from_preamble(src, "w.smk")
    paths = {p.path for p in contract.paths}
    assert ("emseq_ref_assemblies",) in paths
    assert ("emseq_ref_assemblies", "*", "input") in paths


def test_extract_comprehension_values_branch():
    src = "{v['k']: v['w'] for v in config['sec'].values()}\n"
    contract = ecc.extract_from_preamble(src, "w.smk")
    paths = {p.path for p in contract.paths}
    assert ("sec", "*", "k") in paths
    assert ("sec", "*", "w") in paths
