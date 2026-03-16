# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/emseq/emseq.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-16 11:33:40
# ============================================================

# test-analysis.smk — Test wrapper for emseq_analysis.smk
# Tangled from emseq.org; do not edit directly.
import os

configfile: "config/test.yaml"

def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

# --- Environments ---
ENV_EMSEQ = config['envs']['emseq']
ENV_METHYLKIT = config['envs']['methylkit']

# --- Repositories ---
R_EMSEQ = config['repos']['emseq']

# --- Data directories ---
D_DATA = config['main-data-dir']
D_EMSEQ = f"{D_DATA}/emseq"
D_REF = f"{D_DATA}/ref"
D_LOGS = f"{D_DATA}/logs"
D_BENCHMARK = f"{D_DATA}/benchmark"
D_INPUTS = f"{D_DATA}/inputs"

# --- Parameters ---
MOSDEPTH_QUANT_LEVELS = config.get("mosdepth-quant-levels", "1,5,10,20")
EMSEQ_MINCOV = config.get("emseq-mincov", 2)
FASTP_EXTRA = config.get("fastp", {}).get("extra", "")
EMSEQ_REF_INPUTS = {k: v['input'] for k, v in config['emseq_ref_assemblies'].items()}

# --- Samples and references ---
emseq_library_ids = config["library-ids"]
spike_builds = ["puc19", "unmeth_lambda"]
emseq_ref_names = ["chr22"]
KEEP_BED = config["keep-bed"]
EXCL_BED = config["exclude-bed"]
meth_map = config["meth-map"]

# --- Rule all: analysis targets ---
rule all:
    input:
        # DMR - per-base
        expand(
            f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
            experiment=meth_map.keys(),
        ),
        expand(
            f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
            experiment=meth_map.keys(),
        ),
        # DMR - tiled
        expand(
            f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz",
            experiment=meth_map.keys(),
        ),
        # DMR - meth extract
        expand(
            f"{D_EMSEQ}/dmr/diff/{{experiment}}_pos_meth.tsv",
            experiment=meth_map.keys(),
        ),
        # DMR - annotation
        expand(
            f"{D_EMSEQ}/dmr/annotation/{{experiment}}_annotated.tsv",
            experiment=meth_map.keys(),
        ),
        # Haplotype
        expand(
            f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_mhl.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth"],
        ),
        # Deconvolution
        f"{D_EMSEQ}/deconv/uxm_results.csv",

shell.prefix("set -e; ")

# Include both modules
include: "emseq.smk"
include: "emseq_analysis.smk"
