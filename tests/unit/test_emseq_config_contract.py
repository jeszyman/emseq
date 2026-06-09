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


import pathlib


def test_find_includes_resolves_relative_and_warns_on_conditional(tmp_path):
    wrapper = tmp_path / "w.smk"
    wrapper.write_text(
        'include: "emseq.smk"\n'
        "if cond:\n"
        '    include: "maybe.smk"\n'
    )
    (tmp_path / "emseq.smk").write_text("# module\n")
    includes, warnings = ecc.find_includes(wrapper.read_text(), wrapper.parent)
    assert (tmp_path / "emseq.smk").resolve() in includes
    assert any("conditional include" in w for w in warnings)


REPO = pathlib.Path(__file__).resolve().parents[2]


def _paths(contract):
    return {p.path for p in contract.paths}


def test_build_contract_test_smk_snapshot():
    c = ecc.build_contract(REPO / "workflows" / "test.smk")
    p = _paths(c)
    for expect in [("main-data-dir",), ("library-ids",), ("keep-bed",),
                   ("exclude-bed",), ("meth-map",), ("envs", "emseq"),
                   ("envs", "methylkit"), ("repos", "emseq"),
                   ("emseq_ref_assemblies",), ("emseq_ref_assemblies", "*", "input")]:
        assert expect in p, f"missing {expect}"
    byp = {x.path: x for x in c.paths}
    assert byp[("mosdepth-quant-levels",)].requiredness == "optional"
    assert byp[("mosdepth-quant-levels",)].default == "1,5,10,20"
    assert byp[("emseq-mincov",)].default == 2
    # test.smk does NOT include the analysis module: no meth-map inner schema
    assert not any(x[:2] == ("meth-map", "*") for x in p)
    # analysis-only keys absent
    assert ("haplotype", "cpg-ref") not in p
    assert ("deconv", "atlas") not in p


def test_build_contract_analysis_superset_and_methmap_schema():
    c = ecc.build_contract(REPO / "workflows" / "test-analysis.smk")
    p = _paths(c)
    for leaf in ["libs", "tx", "mincov", "mingroup", "chunksize", "win_size",
                 "emseq_ref_name", "align_method"]:
        assert ("meth-map", "*", leaf) in p, f"missing meth-map.*.{leaf}"
    for expect in [("haplotype", "cpg-ref"), ("haplotype", "mhb-bed"),
                   ("haplotype", "metrics"), ("deconv", "genome-name"),
                   ("deconv", "atlas"), ("repos", "mhaptools"),
                   ("envs", "haplotype"), ("envs", "deconv")]:
        assert expect in p, f"missing {expect}"


def test_test_smk_is_subset_of_analysis():
    c_test = _paths(ecc.build_contract(REPO / "workflows" / "test.smk"))
    c_an = _paths(ecc.build_contract(REPO / "workflows" / "test-analysis.smk"))
    assert c_test <= c_an


def test_alias_subkey_word_boundary_and_multibracket():
    aliases = {"meth_map": ("meth-map",)}
    text = (
        'a = meth_map[wc.e]["tx"]\n'
        'b = big_meth_map[wc.e]["libs"]\n'          # substring of alias: must NOT match
        'c = meth_map[wc.e][wc.a]["mincov"]\n'      # two dynamic brackets: must match
    )
    got = {p.path for p in ecc.extract_alias_subkeys(text, aliases, "m.smk")}
    assert ("meth-map", "*", "tx") in got
    assert ("meth-map", "*", "mincov") in got
    assert ("meth-map", "*", "libs") not in got


def test_alias_subkey_scan():
    aliases = {"meth_map": ("meth-map",)}
    module_text = (
        '        library_id=meth_map[wc.experiment]["libs"],\n'
        "        tx = lambda wc: meth_map[wc.experiment]['tx'],\n"
        "        other = some_other[wc.x]['nope'],\n"
    )
    paths = ecc.extract_alias_subkeys(module_text, aliases, "mod.smk")
    got = {p.path for p in paths}
    assert ("meth-map", "*", "libs") in got
    assert ("meth-map", "*", "tx") in got
    assert ("meth-map", "*", "nope") not in got
