# -----------------------------
# Imports
# -----------------------------
import os

# -----------------------------
# Path expansion for strings in config (~, $VARS)
# -----------------------------
def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

# -----------------------------
# Environments
# -----------------------------
ENV_EMSEQ = config['envs']['emseq']
ENV_METHYLKIT = config['envs']['methylkit']

# -----------------------------
# Repositories
# -----------------------------
R_EMSEQ = config['repos']['emseq']

# -----------------------------
# Data directories (derived from main-data-dir)
# -----------------------------
D_DATA = f"{config['main-data-dir']}"
D_EMSEQ = f"{D_DATA}/emseq"
D_REF = f"{D_DATA}/ref"
D_LOGS = f"{D_DATA}/logs"
D_BENCHMARK = f"{D_DATA}/benchmark"
D_INPUTS = f"{D_DATA}/inputs"

# -----------------------------
# Tool/global params
# -----------------------------
mosdepth_quant_levels = config["mosdepth-quant-levels"]
emseq_mincov = 2  # used by per-sample methylKit object rule

# -----------------------------
# Sample set (required by emseq.smk)
# -----------------------------
emseq_library_ids = config["library-ids"]

# -----------------------------
# Reference selections (kept explicit)
# -----------------------------
spike_builds = ["puc19", "unmeth_lambda"]
emseq_ref_names = ["chr22"]

# -----------------------------
# Region filtering inputs
# -----------------------------
KEEP_BED = config["keep-bed"]
EXCL_BED = config["exclude-bed"]

# -----------------------------
# Experiments map from YAML
# (differential methylation experiments for methylKit)
# -----------------------------
meth_map = config["meth-map"]

rule all:
    input:
        # FASTQs
        expand(f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_{{read}}.fastq.gz",
               library_id=emseq_library_ids,
               read=["R1","R2"]),
        #
        # Alignment and methylation calling

        expand(f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
               library_id = emseq_library_ids,
               emseq_ref_name = emseq_ref_names,
               align_method = ["bwa_meth", "biscuit"]),
        #
        # Spike workflow
        expand(f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
               library_id = emseq_library_ids,
               emseq_ref_name = spike_builds,
               align_method = "bwa_meth"),

        # QC
        expand(f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
               library_id = emseq_library_ids,
               processing = ["raw","trimmed"],
               read = ["R1","R2"]),

        expand(f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.summary.txt",
               library_id = emseq_library_ids,
               emseq_ref_name = emseq_ref_names,
               align_method = ["bwa_meth", "biscuit"]),

        expand(f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.txt",
               library_id = emseq_library_ids,
               emseq_ref_name = emseq_ref_names,
               align_method = ["bwa_meth", "biscuit"]),

        f"{D_EMSEQ}/qc/multiqc.html",
        # NOTE: MultiQC rule DOES NOT prompt inputs to run, need individual QC output calls
        #
        # MethylKit
        expand(f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
               experiment = meth_map.keys()),

        expand(f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
               experiment = meth_map.keys()),

        expand(f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}_tiled.txt.bgz",
               experiment = meth_map.keys()),


rule symlink_input_fastqs:
    input:
        r1=f"{D_DATA}/inputs/{{library_id}}.raw_R1.fastq.gz",
        r2=f"{D_DATA}/inputs/{{library_id}}.raw_R2.fastq.gz",
    output:
        r1=f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2=f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    shell:
        r'''
        mkdir -p "$(dirname "{output.r1}")"
        ln -sfr "{input.r1}" "{output.r1}"
        ln -sfr "{input.r2}" "{output.r2}"
        '''

include: f"{R_EMSEQ}/workflows/emseq.smk"
