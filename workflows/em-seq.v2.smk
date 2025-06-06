rule emseq_fastp:
    input:
        r1 = f"{emseq_dir}/fastqs/{{library_id}}_raw_R1.fastq.gz",
        r2 = f"{emseq_dir}/fastqs/{{library_id}}_raw_R2.fastq.gz",
    log:
        cmd = f"{log_dir}/{{library_id}}_emseq_fastp.log",
        json = f"{log_dir}/{{library_id}}_emseq_fastp.json",
        html = f"{log_dir}/{{library_id}}_emseq_fastp.html",
    output:
        r1 = f"{emseq_dir}/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{emseq_dir}/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        failed = f"{emseq_fastq_dir}/{{library_id}}_failed.fastq.gz",
    params:
        script = f"{emseq_script_dir}/fastp_emseq_wrapper.sh",
        threads = 16,
    shell:
        """
        fastp \
        --detect_adapter_for_pe \
        --disable_quality_filtering \
        --failed_out {output.failed} \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --json {log.json} \
        --html {log.html} \
        --out1 {output.r1} \
        --out2 {output.r2} \
        --thread {params.threads} \
        """
# biscuit_ref = config["biscuit_reference"]["ensembl_hg38"]

# biscuit_reference:
#   ensembl_hg38:
#     url: https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#     name: ensembl_hg38
#     fa_name: Homo_sapiens.GRCh38.dna.primary_assembly.fa (because dupsifter d/n recognize .fna)

rule biscuit_index_ensembl_hg38:
    output:
        temp(gz = f"ref/biscuit/{ref['name']}/ref['fa_name']Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
        fa = f"ref/biscuit/{ref['name']}/{ref['fa_name']}",
        fai = f"ref/biscuit-{ref['name']}/{ref['fa_name']}.fai",
        index_done = f"ref/biscuit-{ref['name']}/biscuit_index.done",
    shell:
        r"""
        mkdir -p ref/biscuit/{ref['name']}
        wget -O {temp.gz} {ref['url']}
        gunzip -c {temp.gz} > {output.fa}
        samtools faidx {output.fa}
        biscuit index {output.fa}
        touch {output.index_done}
        """

ref = config["biscuit_reference"]["ensembl_hg38"]
ref_dir = f"ref/biscuit-{ref['name']}"

rule biscuit_index_ensembl_hg38:
    output:
        fa = f"{ref_dir}/{ref['fa_name']}",
        fai = f"{ref_dir}/{ref['fa_name']}.fai",
        index_done = f"{ref_dir}/biscuit_index.done",
    temp:
        gz = f"{ref_dir}/{ref['download_name']}",
    shell:
        r"""
        mkdir -p {ref_dir}
        wget -O {output.gz} {ref['url']}
        gunzip -c {output.gz} > {output.fa}
        samtools faidx {output.fa}
        biscuit index {output.fa}
        touch {output.index_done}
        """
rule emseq_biscuit_align:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        r1 = f"{emseq_dir}/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{emseq_dir}/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        cmd = f"{log_dir}/{{library_id}}_emseq_biscuit_align.log",
    output:
        bam = f"{emseq_dir}/biscuit-bams-raw/{{library_id}}.bam",
    params:
        threads = config["threads"],
    resources:
        concurrency=100
    shell:
        """
        mkdir -p {data_dir}/tmp && \
        biscuit align \
        -@ {params.threads} \
        -biscuit-ref {input.fasta} \
        {input.r1} {input.r2} \
        | samtools sort -n \
        -@ 8 \
        -m 2G \
        -T {data_dir}/tmp/{wildcards.library_id}_sorttmp \
        -o {output.bam} &>> {log}
        """
rule emseq_dedup:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.bam",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        f"{log_dir}/{{library_id}}_emseq_dedup.log",
    output:
        bam = f"{emseq_bam_dir}/{{library_id}}_deduped.bam",
        index = f"{emseq_bam_dir}/{{library_id}}_deduped.bam.bai",
    shell:
        r"""
        rm -f {output.bam}.tmp.*.bam
        samtools view -h -f 0x2 {input.bam} \
        | samtools sort -n -@ 4 -O BAM -o /dev/stdout \
        | dupsifter \
            --add-mate-tags \
            --stats-output {log} \
            {input.fasta} - \
        | samtools sort -@ 8 -o {output.bam}
        samtools index -@ 8 {output.bam}
        """
rule emseq_pileup:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}_deduped.bam",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        f"{log_dir}/{{library_id}}_emseq_pileup.log",
    output:
        vcf = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf.gz",
        tsv = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf_meth_average.tsv",
    params:
        out_base = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf",
    shell:
        """
        biscuit pileup \
	-@ 8 \
	-o {params.out_base} \
        {input.fasta} {input.bam} \
        && bgzip -@ 8 {params.out_base}
        """
rule emseq_post_pileup:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        vcf = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf.gz",
    log:
        f"{log_dir}/{{library_id}}_emseq_post_pileup.log",
    output:
        tbi = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf.gz.tbi",
        bed = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.bed",
        bismark = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_bismark_cov.bed",
    shell:
        """
        tabix -p vcf {input.vcf} \
        && biscuit vcf2bed \
	-t cg {input.vcf} > {output.bed} \
        && biscuit vcf2bed -c {input.vcf} > {output.bismark}
        """
rule make_single_methylkit_obj:
    conda:
        "../config/methylkit-conda-env.yaml",
    input:
        bismark = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_bismark_cov.bed",
    log:
        f"{log_dir}/methylkit_{{library_id}}.log",
    output:
        txt = f"{emseq_dir}/dmr/tabix/{{library_id}}.txt",
        bgz = f"{emseq_dir}/dmr/tabix/{{library_id}}.txt.bgz",
        tbi = f"{emseq_dir}/dmr/tabix/{{library_id}}.txt.bgz.tbi",
    params:
        Rscript = f"{emseq_script_dir}/make_single_methylkit_obj.R",
        out_dir = f"{emseq_dir}/dmr/tabix",
        mincov = emseq_mincov,
        build = emseq_build,
        treatment = 1,
    shell:
        """
        Rscript {params.Rscript} \
          --bismark_cov_bed {input.bismark} \
          --library_id {wildcards.library_id} \
          --mincov {params.mincov} \
          --out_dir {params.out_dir} \
          --treatment {params.treatment} \
          --build {params.build} \
          &>> {log}
        """
rule make_methylkit_diff_db:
    input:
        mkit_lib_db = lambda wildcards: expand(
            f"{emseq_dir}/dmr/tabix/{{library_id}}.txt.bgz",
            library_id = meth_map[wildcards.experiment]['libs']
        ),
    log:
        f"{log_dir}/{{experiment}}_make_methylkit_diff_db.log",
    output:
        unite = f"{emseq_dir}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
        diff = f"{emseq_dir}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
    params:
        library_id = lambda wildcards: " ".join(meth_map[wildcards.experiment]['libs']),
        treatment_list = lambda wildcards: meth_map[wildcards.experiment]['tx'],
        out_dir = f"{emseq_dir}/dmr/diff",
        script = f"{emseq_script_dir}/make_methylkit_diff_db.R",
    shell:
        """
        Rscript {params.script} \
        --lib_db_list "{input.mkit_lib_db}" \
        --lib_id_list "{params.library_id}" \
        --treatment_list "{params.treatment_list}" \
        --cores 32 \
        --out_dir {params.out_dir} \
        --suffix {wildcards.experiment} > {log} 2>&1
        """
rule all_experiment_methylation:
    input:
        f"{emseq_dir}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
    log:
        f"{log_dir}/all_experiment_methylation_{{experiment}}.log",
    output:
        f"{emseq_dir}/dmr/diff/{{experiment}}_pos_meth.tsv",
    params:
        script = f"{emseq_script_dir}/all_experiment_methylation.R",
    shell:
        """
        Rscript {params.script} \
        --db_file {input} \
        --out_file {output} > {log} 2>&1
        """
rule make_methylkit_diff_db_tiled:
    input:
        mkit_lib_db = lambda wildcards: expand(
            f"{emseq_dir}/dmr/tabix/{{library_id}}.txt.bgz",
            library_id = meth_map[wildcards.experiment]['libs']
        ),
    log:
        f"{log_dir}/{{experiment}}_make_methylkit_diff_tiled_db.log",
    output:
        unite = f"{emseq_dir}/dmr/diff/methylBase_{{experiment}}_tiled.txt.bgz",
        diff = f"{emseq_dir}/dmr/diff/methylDiff_{{experiment}}_tiled.txt.bgz",
    params:
        library_id = lambda wildcards: " ".join(meth_map[wildcards.experiment]['libs']),
        treatment_list = lambda wildcards: meth_map[wildcards.experiment]['tx'],
        out_dir = f"{emseq_dir}/dmr/diff",
        script = f"{emseq_script_dir}/make_methylkit_diff_tiled_db.R",
        win_size = 1000000,
        step_size= 1000000,
    shell:
        """
        Rscript {params.script} \
        --lib_db_list "{input.mkit_lib_db}" \
        --lib_id_list "{params.library_id}" \
        --treatment_list "{params.treatment_list}" \
        --cores 32 \
        --out_dir {params.out_dir} \
        --win_size {params.win_size} \
        --step_size {params.step_size} \
        --suffix {wildcards.experiment} \
        > {log} 2>&1
        """
rule all_experiment_tiled_methylation:
    input:
        f"{emseq_dir}/dmr/diff/methylBase_{{experiment}}_tiled.txt.bgz",
    log:
        f"{log_dir}/all_experiment_methylation_{{experiment}}_tiled.log",
    output:
        f"{emseq_dir}/dmr/diff/{{experiment}}_tiled_meth.tsv",
    params:
        script = f"{emseq_script_dir}/all_experiment_methylation.R",
    shell:
        """
        Rscript {params.script} \
        --db_file {input} \
        --out_file {output} > {log} 2>&1
        """
rule emseq_fastqc:
    input:
        f"{emseq_fastq_dir}/{{library_id}}_{{processing}}_{{read}}.fastq.gz",
    log:
        f"{log_dir}/{{library_id}}_{{processing}}_{{read}}_fastqc.log",
    output:
        f"{qc_dir}/{{library_id}}_{{processing}}_{{read}}_fastqc.html",
        f"{qc_dir}/{{library_id}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = qc_dir,
        threads = 2,
    resources:
        concurrency=20
    shell:
        """
        fastqc \
        --outdir {params.outdir} \
        --quiet \
        --svg \
        --threads {params.threads} \
        {input} &> {log}
        """
# Will follow symlinks
# rule emseq_index_bam_check:
#     input:
#         bam = ancient(f"{emseq_bam_dir}/{{library_id}}_deduped.bam"),
#     output:
#         bai = f"{emseq_bam_dir}/{{library_id}}_deduped.bam.bai",
#     shell:
#         """
#         samtools index -@ 8 {input.bam} {output.bai}
#         """

rule emseq_mosdepth:
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}_deduped.bam",
        index = f"{emseq_bam_dir}/{{library_id}}_deduped.bam.bai",
    output:
        summary = f"{qc_dir}/mosdepth_{{library_id}}.mosdepth.summary.txt",
        global_dist = f"{qc_dir}/mosdepth_{{library_id}}.mosdepth.global.dist.txt",
        region_dist = f"{qc_dir}/mosdepth_{{library_id}}.mosdepth.region.dist.txt",
        regions = f"{qc_dir}/mosdepth_{{library_id}}.regions.bed.gz",
        regions_idx = f"{qc_dir}/mosdepth_{{library_id}}.regions.bed.gz.csi",
        quantized = f"{qc_dir}/mosdepth_{{library_id}}.quantized.bed.gz",
        quantized_idx = f"{qc_dir}/mosdepth_{{library_id}}.quantized.bed.gz.csi",
        thresholds = f"{qc_dir}/mosdepth_{{library_id}}.thresholds.bed.gz",
        thresholds_idx = f"{qc_dir}/mosdepth_{{library_id}}.thresholds.bed.gz.csi",
    params:
        script = f"{emseq_script_dir}/emseq_mosdepth.sh",
        quant_levels = mosdepth_quant_levels,
        out_dir = qc_dir,
    threads: 8
    shell:
        """
        {params.script} \
        {input.bam} \
        {params.out_dir} \
        {wildcards.library_id} \
        "{params.quant_levels}" \
        {threads}
        """
print("emseq_library_ids:", emseq_library_ids)
print("type of first item:", type(emseq_library_ids[0]))

rule emseq_mosdepth_agg_plot:
    input:
        thresholds = expand(f"{qc_dir}/mosdepth_{{library_id}}.thresholds.bed.gz", library_id=emseq_library_ids),
        regions = expand(f"{qc_dir}/mosdepth_{{library_id}}.regions.bed.gz", library_id=emseq_library_ids),
    output:
        pdf = f"{qc_dir}/mosdepth_agg_plot.pdf",
    params:
        script = f"{emseq_script_dir}/emseq_mosdepth_agg_plot.R",
        library_list = " ".join(emseq_library_ids),
        threshold_list = lambda wildcards, input: " ".join(input.thresholds),
        regions_list = lambda wildcards, input: " ".join(input.regions),
    shell:
        """
        Rscript {params.script} \
        --threshold_list "{params.threshold_list}" \
        --regions_list "{params.regions_list}" \
        --library_list "{params.library_list}" \
        --output_pdf {output.pdf}
        """
