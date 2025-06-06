# Development snakemake to test
mosdepth_quant_levels = config["mosdepth-quant-levels"]
repo = "~/repos/emseq"
data_dir=config["data_dir"]
emseq_script_dir = "~/repos/emseq/scripts"
log_dir = f"{data_dir}/logs"

emseq_mincov = 2
emseq_build = "hg38"

# Explicitly select which references to build
index_targets = ["unmeth_lambda", "puc19"]

library_ids = ["NH_15_L3", "PRO_6_L2"]

spike_ref_names = ["unmeth_lambda"]

ref_names = ["ncbi_decoy_hg38"]

align_methods = ["bwa_meth"]

emseq_library_ids = library_ids

meth_map = {
    "test": {
        "build": "hg38",
        "mincov": "5",
        "libs": library_ids,
        "tx": "0,1"
}
}

rule all:
    input:
        # FASTQs
        expand(f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_{{processing}}_{{read}}.fastq.gz",
               library_id = library_ids,
               processing = ["raw","trimmed"],
               read = ["R1","R2"]),

        expand(f"{data_dir}/qc/{{library_id}}_{{processing}}_{{read}}_fastqc.{{suffix}}",
               library_id = library_ids,
               processing = ["raw","trimmed"],
               read = ["R1","R2"],
               suffix = ["zip","html"]),

        # Spike-ins
        expand(f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam",
               library_id = library_ids,
               ref_name = spike_ref_names,
               align_method = "bwa_meth"),

        expand(f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
               library_id=library_ids,
               ref_name=spike_ref_names,
               align_method= "bwa_meth"),

        # Quick bwa-meth alignment

        expand(f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam",
               library_id = library_ids,
               ref_name = ref_names,
               align_method = "bwa_meth"),

        # Deduplicate
        expand(f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
               library_id = library_ids,
               ref_name = ref_names,
               align_method = "bwa_meth"),


        # expand(f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
        #        library_id = library_ids,
        #        ref_name = ref_names,
        #        align_method = "bwa_meth"),

        # expand(f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.txt",
        #        library_id = library_ids,
        #        ref_name = ref_names,
        #        align_method = "bwa_meth"),



include: "./dev.smk"
