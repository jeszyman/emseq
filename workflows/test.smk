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
D_DATA = config['main-data-dir']
D_EMSEQ = f"{D_DATA}/emseq"
D_REF = f"{D_DATA}/ref"
D_LOGS = f"{D_DATA}/logs"
D_BENCHMARK = f"{D_DATA}/benchmark"
D_INPUTS = f"{D_DATA}/inputs"

# -----------------------------
# Tool/global params (UPPERCASE for module consumption)
# -----------------------------
MOSDEPTH_QUANT_LEVELS = config.get("mosdepth-quant-levels", "1,5,10,20")
EMSEQ_MINCOV = config.get("emseq-mincov", 2)
FASTP_EXTRA = config.get("fastp", {}).get("extra", "")

# -----------------------------
# Reference assembly lookup (for index rules)
# -----------------------------
EMSEQ_REF_INPUTS = {k: v['input'] for k, v in config['emseq_ref_assemblies'].items()}

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

# -----------------------------
# Rule all
# -----------------------------
rule all:
    input:
        # FASTQs
        expand(
            f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_{{read}}.fastq.gz",
            library_id=emseq_library_ids,
            read=["R1", "R2"],
        ),
        # Alignment and methylation calling
        expand(
            f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth", "biscuit"],
        ),
        # Spike workflow
        expand(
            f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
            library_id=emseq_library_ids,
            emseq_ref_name=spike_builds,
            align_method="bwa_meth",
        ),
        # QC - FastQC
        expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
            library_id=emseq_library_ids,
            processing=["raw", "trimmed"],
            read=["R1", "R2"],
        ),
        # QC - mosdepth
        expand(
            f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.summary.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth", "biscuit"],
        ),
        # QC - M-bias
        expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth", "biscuit"],
        ),
        # QC - MultiQC (does NOT prompt inputs to run)
        f"{D_EMSEQ}/qc/multiqc.html",
        # MethylKit - per-base
        expand(
            f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
            experiment=meth_map.keys(),
        ),
        expand(
            f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
            experiment=meth_map.keys(),
        ),
        # MethylKit - tiled
        expand(
            f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz",
            experiment=meth_map.keys(),
        ),

# -----------------------------
# Input symlinks
# -----------------------------
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

# -----------------------------
# Include module
# -----------------------------
include: "emseq.smk"
