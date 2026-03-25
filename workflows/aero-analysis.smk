# aero-analysis.smk — AERO project wrapper (core + analysis)
# Hand-crafted for local run; not tangled from org.
import os

configfile: "config/aero.yaml"

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
emseq_ref_names = ["ncbi_decoy_hg38"]
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

# --- Input symlinks ---
rule symlink_input_fastqs:
    message: "Create symlinks for raw input FASTQs into workflow directory"
    input:
        r1 = f"{D_INPUTS}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_INPUTS}/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_symlink_input_fastqs.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_symlink_input_fastqs.tsv"
    output:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[symlink-fastqs] $(date) lib={wildcards.library_id}"
        mkdir -p "$(dirname "{output.r1}")"
        ln -sfr "{input.r1}" "{output.r1}"
        ln -sfr "{input.r2}" "{output.r2}"
        """

# Include both modules
include: "emseq.smk"
include: "emseq_analysis.smk"
