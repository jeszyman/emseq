# Development snakemake to test
mosdepth_quant_levels = config["mosdepth-quant-levels"]
repo = "~/repos/emseq"
data_dir=config["data_dir"]
emseq_script_dir = "~/repos/emseq/scripts"
log_dir = f"{data_dir}/logs"

emseq_mincov = 2
emseq_build = "hg38"

threads = 80
# We specify em-seq bam directory directly to allow for workflows that merge at the bam level:
emseq_bam_dir = f"{data_dir}/analysis/emseq/bams"


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


mosdepth_map = {
    "tests": {
        "library_ids": ["NH_15_L3", "PRO_6_L2"],
        "ref_name": "ncbi_decoy_hg38",
        "align_method": "bwa_meth"
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
        expand(f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam",
               library_id = library_ids,
               ref_name = spike_ref_names,
               align_method = "bwa_meth"),

        expand(f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
               library_id=library_ids,
               ref_name=spike_ref_names,
               align_method= "bwa_meth"),

        # Biscuit 1 steps
        ## Index
        expand(f"{data_dir}/ref/biscuit/{{name}}/{{name}}.fa.fai",
               name = "ncbi_decoy_hg38"),

        ## Align
        expand(f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.biscuit.coorsort.bam",
               library_id = library_ids,
               ref_name = ref_names),

        # BWA-meth
        # ## Index
        # expand(f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t",
        #        name = ref_names),

        # ## Align

        # expand(f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam",
        #        library_id = library_ids,
        #        ref_name = ref_names,
        #        align_method = "bwa_meth"),

        ## COMMON DEDUP HERE

        ## Pileup
        # expand(f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}.biscuit_pileup.{{suffix}}",
        #        library_id = library_ids,
        #        ref_name = ref_names,
        #        suffix = ["vcf.gz","vcf_meth_average.tsv"]),

        # ## Make per-library methylkit objects
        # expand(f"{data_dir}/analysis/emseq/post-biscuit/{{library_id}}.{{ref_name}}_biscuit.{{suffix}}",
        #        library_id = library_ids,
        #        ref_name = ref_names,
        #        suffix = ["txt", "txt.bgz", "txt.bgz.tbi"]),

        # Common post-alignment per-library steps
        ## Deduplicate
        expand(f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
               library_id = library_ids,
               ref_name = ref_names,
               align_method = ["biscuit","bwa_meth"]),

        ## Depth
        expand(f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.summary.txt",
               library_id = library_ids,
               ref_name = ref_names,
               align_method = ["biscuit","bwa_meth"]),

        ## Call methylation
        expand(f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
               library_id = library_ids,
               ref_name = ref_names,
               align_method = ["biscuit", "bwa_meth"]),

        ## Create per-library methylkit objects
        expand(f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.txt",
               library_id = library_ids,
               ref_name = ref_names,
               align_method = "bwa_meth"),

        expand(f"{data_dir}/qc/{{experiment}}.emseq_mosdepth_agg_plot.pdf",
               experiment = mosdepth_map.keys())


include: "./dev.smk"
