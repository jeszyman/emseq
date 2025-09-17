import pandas as pd
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

ENV_EMSEQ = config['ENV-EMSEQ']

D_EMSEQ = f"{config['D-DATA']}/emseq"
D_REF = f"{config['D-DATA']}/ref"
D_LOGS = f"{config['D-DATA']}/logs"
D_BENCHMARK = f"{config['D-DATA']}/benchmark"
D_INPUTS = f"{config['D-DATA']}/inputs"

R_EMSEQ = config['R-EMSEQ']

emseq_mincov = 2

emseq_build = ["ncbi_decoy_hg38"]

spike_builds = ["puc19", "unmeth_lambda"]
library_ids = ["AERO_24"]

rule all:
    input:
        # FASTQs
        expand(f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
               library_id = library_ids,
               processing = ["trimmed"],
               read = ["R1","R2"]),

        expand(f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.methyldackel.txt.bgz",
               library_id = library_ids,
               ref_name = emseq_build,
               align_method = "bwa_meth"),

        expand(f"{D_EMSEQ}/spike/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
               library_id = library_ids,
               ref_name = spike_builds,
               align_method = "bwa_meth"),

include: f"{R_EMSEQ}/workflows/emseq.smk"
