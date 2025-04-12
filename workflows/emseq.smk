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
