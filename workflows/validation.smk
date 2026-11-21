# validation.smk: wrapper for validations on data held outside the repository
# Tangled from emseq.org; do not edit directly.
#
# The config is given on the command line (--configfile). A cohort's config
# stays in that cohort's project repository.
import os

def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

# Concurrency budget: defaults from the available-concurrency config key (100 if unset);
# override per run with --resources concurrency=N (the CLI value wins over config).
workflow.global_resources.setdefault("concurrency", int(config.get("available-concurrency", 100)))

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

# --- Tool parameters ---
MOSDEPTH_QUANT_LEVELS = config.get("mosdepth-quant-levels", "1,5,10,20")
EMSEQ_MINCOV = config.get("emseq-mincov", 2)
FASTP_EXTRA = config.get("fastp", {}).get("extra", "")
EMSEQ_REF_INPUTS = {k: v['input'] for k, v in config['emseq_ref_assemblies'].items()}

# --- Samples and references ---
emseq_library_ids = config["library-ids"]
emseq_ref_names = [config["validation-host-ref"]]
emseq_align_methods = ["bwa_meth"]
spike_builds = ["puc19", "unmeth_lambda"]
KEEP_BED = config["keep-bed"]
EXCL_BED = config["exclude-bed"]
meth_map = config.get("meth-map", {})

# --- Validation ---
VALIDATION_DATASET = config["validation-dataset"]
VALIDATION_INDEX = config["validation-library-index"]
VALIDATION_HOST_REF = config["validation-host-ref"]
VALIDATION_FRAGMENT_SOURCES = config.get("validation-fragment-end-sources", [])
VALIDATION_CHECK_INDEX = config.get("validation-check-index")
# Read-end settings compared at calling, with flags that
# tests/unit/test_methyldackel_bounds.py pins for MethylDackel 0.6.1
VALIDATION_END_SETTINGS = {
    "none": "",
    "neb": "--nOT 0,0,0,5 --nOB 0,0,5,0",
    "all10": "--nOT 10,11,10,11 --nOB 10,11,10,11",
} if config.get("validation-compare-settings", False) else {}
# pUC19 is methylated at CpG only; lambda is unmethylated in all three contexts
VALIDATION_CONTROL_CONTEXTS = {"puc19": ["CpG"], "unmeth_lambda": ["CpG", "CHG", "CHH"]}

_unindexed = [lib for lib in emseq_library_ids if lib not in VALIDATION_INDEX]
if _unindexed:
    raise ValueError(f"validation-library-index has no entry for: {', '.join(_unindexed)}")

# Panels of the strand-end figure: this data set's sources, then panels read
# from the committed tables of other data sets
VALIDATION_FRAGMENT_PANELS = [
    (panel.get("dataset", VALIDATION_DATASET), panel["source"], panel["label"], panel.get("exclude-index", 0))
    for panel in config.get("validation-fragment-figure", [])
]

# --- Rule all ---
rule all:
    input:
        f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}/values.json",
        f"{R_EMSEQ}/resources/figures/validation_end_handling_{VALIDATION_DATASET}.png",
        [f"{R_EMSEQ}/resources/figures/validation_fragment_end_{VALIDATION_DATASET}.png"] if VALIDATION_FRAGMENT_PANELS else [],

# -----------------------------
# Include modules
# -----------------------------
include: "emseq.smk"
include: "emseq_validation.smk"
