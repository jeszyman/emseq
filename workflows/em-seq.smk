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
rule methylkit_dmr_obj:
    input:
        bismark_cov lambda wildcards: expand(f"{emseq_dir}/pileup/{{library_id}}_bismark_cov.bed",
                                             library = emseq_map[wildcards.experiment]['libs']),
    log:
    output:
        f"{}
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
rule emseq_multiqc:
    input:
        expand(f"{qc_dir}/{{library_id}}_{{processing}}_{{read}}_fastqc.zip",
               library_id = library_ids,
               processing = ["raw","trimmed"],
               read = ["R1", "R2"]),
    log:
        f"{log_dir}/emseq_multiqc.log",
    output:
        f"{qc_dir}/emseq_multiqc/emseq_multiqc.html",
    params:
        out_dir = f"{qc_dir}/emseq_multiqc",
        out_name = "emseq_multiqc",
    shell:
        """
        multiqc \
        {input} \
        --force \
        --outdir {params.out_dir} \
        --filename {params.out_name}
        """
