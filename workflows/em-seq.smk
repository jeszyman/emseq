rule emseq_fastp:
    input:
        r1 = f"{emseq_fastq_dir}/{{library_id}}_raw_R1.fastq.gz",
        r2 = f"{emseq_fastq_dir}/{{library_id}}_raw_R2.fastq.gz",
    log:
        cmd = f"{log_dir}/{{library_id}}-emseq-fastp.log",
        json = f"{log_dir}/{{library_id}}-emseq-fastp.json",
        html = f"{log_dir}/{{library_id}}-emseq-fastp.html",
    output:
        r1 = f"{emseq_fastq_dir}/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{emseq_fastq_dir}/{{library_id}}_trimmed_R2.fastq.gz",
        failed = f"{emseq_fastq_dir}/{{library_id}}_failed.fastq.gz",
    params:
        script = f"{emseq_script_dir}/fastp-emseq-wrapper.sh",
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
rule emseq_biscuit_align:
    input:
        r1 = f"{emseq_fastq_dir}/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{emseq_fastq_dir}/{{library_id}}_trimmed_R2.fastq.gz",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        cmd = f"{log_dir}/{{library_id}}_emseq_biscuit_align.log",
    output:
        bam = f"{emseq_bam_dir}/{{library_id}}.bam",
    resources:
        concurrency=100
    shell:
        """
        mkdir -p {data_dir}/tmp && \
        biscuit align \
        -@ 82 \
        -biscuit-ref {input.fasta} \
        {input.r1} {input.r2} \
        | samtools sort \
        -@ 8 \
        -m 2G \
        -T {data_dir}/tmp/{wildcards.library_id}_sorttmp \
        -o {output.bam} &>> {log}
        """
rule emseq_dedup:
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.bam",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        f"{log_dir}/{{library_id}}_emseq_dedup.log",
    output:
        bam = f"{emseq_bam_dir}/{{library_id}}_deduped.bam",
        index = f"{emseq_bam_dir}/{{library_id}}_deduped.bam.bai",
    shell:
        """
        dupsifter \
        --add-mate-tags \
        --stats-output {log} \
        {input.fasta} \
        {input.bam} \
        | samtools sort \
	-o {output.bam} \
	-@ 8 && samtools index -@ 8 {output.bam}
        """
rule emseq_pileup:
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
        threads = threads,
    shell:
        """
        fastqc \
        --outdir {params.outdir} \
        --quiet \
        --svg \
        --threads {params.threads} \
        {input} &> {log}
        """
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
        "{params.script}" \
        "{input.bam}" \
        "{params.out_dir}" \
        "{wildcards.library_id}" \
        "{params.quant_levels}" \
        {threads}
        """
rule make_single_methylkit_obj:
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
        --cores 8 \
        --out_dir {params.out_dir} \
        --suffix {wildcards.experiment} > {log} 2>&1
        """
