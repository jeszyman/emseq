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
ENV_HAPLOTYPE = config['envs']['haplotype']
ENV_DECONV = config['envs']['deconv']

# --- Repositories ---
R_EMSEQ = config['repos']['emseq']
R_MHAPTOOLS = config['repos']['mhaptools']
R_WGBSTOOLS = config['repos']['wgbs_tools']
R_UXM = config['repos']['uxm_deconv']

# --- Data directories ---
D_DATA = config['main-data-dir']
D_EMSEQ = f"{D_DATA}/emseq"
D_REF = f"{D_DATA}/ref"
D_LOGS = f"{D_DATA}/logs"
D_BENCHMARK = f"{D_DATA}/benchmark"
D_INPUTS = f"{D_DATA}/inputs"

# --- Tool parameters ---
MOSDEPTH_QUANT_LEVELS = config.get("mosdepth-quant-levels", "1,5,10,20")
EMSEQ_MINCOV = config.get("emseq-mincov", 2)
FASTP_EXTRA = config.get("fastp", {}).get("extra", "")
EMSEQ_REF_INPUTS = {k: v['input'] for k, v in config['emseq_ref_assemblies'].items()}

# --- Samples and references ---
emseq_library_ids = config["library-ids"]
emseq_ref_names = ["chr22"]
emseq_align_methods = ["bwa_meth", "biscuit"]
spike_builds = ["puc19", "unmeth_lambda"]
KEEP_BED = config["keep-bed"]
EXCL_BED = config["exclude-bed"]
meth_map = config["meth-map"]

# --- Analysis-specific variables ---
CPG_REF = config["haplotype"]["cpg-ref"]
MHB_BED = config["haplotype"]["mhb-bed"]
HAP_METRICS = " ".join(config["haplotype"]["metrics"])
DECONV_GENOME = config["deconv"]["genome-name"]
DECONV_ATLAS = config["deconv"]["atlas"]

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
        # Haplotype — excluded from test: mHapTools requires C++ compilation
        # from source with bundled htslib; too fragile for CI. Validated on
        # production data.
        # expand(
        #     f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_mhl.txt",
        #     library_id=emseq_library_ids,
        #     emseq_ref_name=emseq_ref_names,
        #     align_method=["bwa_meth"],
        # ),
        # Deconvolution — excluded from test: UXM requires more CpG coverage
        # than chr22 test data provides. Validated on production data.
        # f"{D_EMSEQ}/deconv/uxm_results.csv",

shell.prefix("set -e; ")

# Test alias: duplicate lib003 as lib003b for unique methylKit sample IDs
rule alias_lib003b:
    input:
        r1 = f"{D_INPUTS}/lib003.raw_R1.fastq.gz",
        r2 = f"{D_INPUTS}/lib003.raw_R2.fastq.gz",
    output:
        r1 = f"{D_INPUTS}/lib003b.raw_R1.fastq.gz",
        r2 = f"{D_INPUTS}/lib003b.raw_R2.fastq.gz",
    shell:
        """
        ln -sfr "{input.r1}" "{output.r1}"
        ln -sfr "{input.r2}" "{output.r2}"
        """

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
