# aero.smk — AERO project wrapper (core processing)
# Hand-crafted for local run; not tangled from org.
import os
import re

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

# --- Resource helpers ---
# Defaults mirror the original hardcoded values so emseq.smk works without config
_RESOURCE_DEFAULTS = {
    "align":             {"threads": 48, "concurrency": 50},
    "align-spike":       {"threads": 48, "concurrency": 50},
    "dedup":             {"threads": 8,  "concurrency": 25},
    "filter-bam":        {"threads": 16, "concurrency": 25},
    "fastp":             {"threads": 8,  "concurrency": 50},
    "fastqc":            {"threads": 4,  "concurrency": 25},
    "mosdepth":          {"threads": 8,  "concurrency": 20},
    "mbias":             {"threads": 10, "concurrency": 50},
    "methyldackel":      {"threads": 20, "concurrency": 25},
    "methyldackel-spike":{"threads": 8,  "concurrency": 50},
    "samtools-stats":    {"threads": 8,  "concurrency": 40},
    "methylkit":         {"threads": 1,  "concurrency": 50},
    "multiqc":           {"threads": 4,  "concurrency": 20},
}
_RES_CFG = config.get("resources", {})

def rule_threads(name):
    """Return min(configured threads, workflow.cores) for a rule class."""
    t = _RES_CFG.get(name, {}).get("threads", _RESOURCE_DEFAULTS[name]["threads"])
    return min(t, workflow.cores)

def rule_concurrency(name):
    """Return configured concurrency for a rule class."""
    return _RES_CFG.get(name, {}).get("concurrency", _RESOURCE_DEFAULTS[name]["concurrency"])

# --- Samples and references ---
emseq_library_ids = config["library-ids"]
spike_builds = ["puc19", "unmeth_lambda"]
emseq_ref_names = ["ncbi_decoy_hg38"]
KEEP_BED = config["keep-bed"]
EXCL_BED = config["exclude-bed"]
meth_map = config["meth-map"]

# --- Reference download lookup ---
# Map input filename → URL for each reference assembly
_REF_URL_MAP = {v['input']: v['url'] for v in config['emseq_ref_assemblies'].values()}

# Blacklist URL (same source as get_test_data.sh)
_BLACKLIST_URL = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"

# Autosome contig pattern for keep bed generation
_AUTOSOME_PATTERN = r'^chr[0-9]+\t'

# --- Rule all: core processing targets ---
rule all:
    input:
        # Trimmed FASTQs
        expand(
            f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_{{read}}.fastq.gz",
            library_id=emseq_library_ids,
            read=["R1", "R2"],
        ),
        # Alignment and methylation calling (bwa_meth only)
        expand(
            f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth"],
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
            align_method=["bwa_meth"],
        ),
        # QC - M-bias
        expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth"],
        ),
        # QC - MultiQC (skipped — emseq_multiqc hardcodes biscuit align_method)
        # f"{D_EMSEQ}/qc/multiqc.html",

shell.prefix("set -e; ")

# -----------------------------
# Reference downloads
# -----------------------------
rule download_ref_fasta:
    message: "Download reference FASTA from upstream URL"
    wildcard_constraints:
        ref_input = "|".join(re.escape(k) for k in _REF_URL_MAP),
    output:
        f"{D_INPUTS}/{{ref_input}}"
    params:
        url = lambda wc: _REF_URL_MAP[wc.ref_input],
    log:
        cmd = f"{D_LOGS}/download_{{ref_input}}.log",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[download-ref] $(date) file={wildcards.ref_input} url={params.url}"
        mkdir -p "$(dirname "{output}")"
        TMP="{output}.tmp"
        curl -fsSL "{params.url}" -o "$TMP"
        # If source is plain text FASTA, gzip it; otherwise keep as-is
        if file -b "$TMP" | grep -qi gzip; then
            mv "$TMP" "{output}"
        else
            gzip -c "$TMP" > "{output}"
            rm "$TMP"
        fi
        """

rule download_blacklist:
    message: "Download ENCODE hg38 blacklist BED"
    output:
        EXCL_BED,
    params:
        url = _BLACKLIST_URL,
    log:
        cmd = f"{D_LOGS}/download_blacklist.log",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[download-blacklist] $(date) url={params.url}"
        mkdir -p "$(dirname "{output}")"
        curl -fsSL "{params.url}" -o "{output}"
        """

rule build_keep_bed:
    message: "Generate autosomes keep BED from reference FASTA index"
    input:
        fai = f"{D_REF}/bwa_meth/ncbi_decoy_hg38/ncbi_decoy_hg38.fa.fai",
    output:
        KEEP_BED,
    log:
        cmd = f"{D_LOGS}/build_keep_bed.log",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[build-keep-bed] $(date)"
        mkdir -p "$(dirname "{output}")"
        awk 'BEGIN{{OFS="\\t"}} /^chr[0-9]+\\t/{{print $1,0,$2}}' "{input.fai}" \
          | sort -k1,1V -k2,2n > "{output}"
        """

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

# Include core processing module
include: "emseq.smk"
