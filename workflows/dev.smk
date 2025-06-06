rule emseq_fastp:
    input:
        r1 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_raw_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_raw_R2.fastq.gz",
    log:
        cmd = f"{log_dir}/{{library_id}}-emseq-fastp.log",
        json = f"{log_dir}/{{library_id}}-emseq-fastp.json",
        html = f"{log_dir}/{{library_id}}-emseq-fastp.html",
    output:
        r1 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        failed = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_failed.fastq.gz",
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
rule emseq_fastqc:
    input:
        f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_{{processing}}_{{read}}.fastq.gz",
    log:
        f"{data_dir}/logs/{{library_id}}_{{processing}}_{{read}}_fastqc.log",
    output:
        f"{data_dir}/qc/{{library_id}}_{{processing}}_{{read}}_fastqc.html",
        f"{data_dir}/qc/{{library_id}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = f"{data_dir}/qc",
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
        bam = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    output:
        summary = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.summary.txt",
        global_dist = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.global.dist.txt",
        region_dist = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.region.dist.txt",
        regions = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.regions.bed.gz",
        regions_idx = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.regions.bed.gz.csi",
        quantized = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.quantized.bed.gz",
        quantized_idx = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.quantized.bed.gz.csi",
        thresholds = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.thresholds.bed.gz",
        thresholds_idx = f"{data_dir}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.thresholds.bed.gz.csi",
    params:
        script = f"{emseq_script_dir}/emseq_mosdepth.sh",
        quant_levels = mosdepth_quant_levels,
        out_dir = f"{data_dir}/qc",
    threads: 8
    shell:
        """
        {params.script} \
        {input.bam} \
        {params.out_dir} \
        {wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method} \
        "{params.quant_levels}" \
        {threads}
        """
print("emseq_library_ids:", emseq_library_ids)
print("type of first item:", type(emseq_library_ids[0]))

rule emseq_mosdepth_agg_plot:
    input:
        thresholds = expand(f"{data_dir}/qc/mosdepth_{{library_id}}.thresholds.bed.gz", library_id=emseq_library_ids),
        regions = expand(f"{data_dir}/qc/mosdepth_{{library_id}}.regions.bed.gz", library_id=emseq_library_ids),
    output:
        pdf = f"{data_dir}/qc/mosdepth_agg_plot.pdf",
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
rule bwa_meth_index:
    conda:
        "../config/emseq-conda-env.yaml"
    input:
        lambda wildcards: f"{data_dir}/inputs/{config['emseq_ref_assemblies'][wildcards.name]['input']}"
    output:
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t",
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.0123",
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.amb",
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.ann",
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.bwt.2bit.64",
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.pac",
        f"{data_dir}/ref/bwa_meth/{{name}}/{{name}}.fa",
    params:
        fasta_target = lambda wildcards: f"{data_dir}/ref/bwa_meth/{wildcards.name}/{wildcards.name}.fa"
    log:
        f"{log_dir}/{{name}}_bwa_meth_index.log"
    shell:
        """
        mkdir -p $(dirname {params.fasta_target}) && \
        zcat {input} > {params.fasta_target} && \
        samtools faidx -@ 8 {params.fasta_target} && \
        bwameth.py index-mem2 {params.fasta_target} > {log} 2>&1
        """
rule emseq_align_bwameth:
    conda:
        "../config/emseq-conda-env.yaml"
    input:
        r1 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        ref = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
        c2t = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t",
    output:
        bam = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam"
    threads: 16
    params:
        temp_prefix = lambda wildcards: f"{data_dir}/tmp/{wildcards.library_id}.{wildcards.ref_name}",
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.bwameth.log"
    shell:
        """
        bwameth.py --threads {threads} \
            --reference {input.ref} \
            {input.r1} {input.r2} 2> {log} | \
        samtools view -u - | \
        samtools sort -@ {threads} -T {params.temp_prefix} -o {output.bam}
        """
rule emseq_dedup_bwameth:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam",
        fasta = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    log:
        f"{data_dir}/logs/{{library_id}}.{{ref_name}}.emseq_dedup_bwameth.log",
    output:
        bam = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.deduped.bam",
        index = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.deduped.bam.bai",
    params:
        temp_prefix = lambda wildcards: f"{data_dir}/tmp/{wildcards.library_id}.{wildcards.ref_name}.namesort",
    shell:
        r"""
        set -eo pipefail

        # Extract only properly paired reads (-f 0x2), include header (-h)
        # Name-sort BAM using samtools; output to stdout
        # Deduplicate using dupsifter, using reference FASTA and streaming input from stdin
        # Coordinate-sort deduplicated BAM for downstream tools
        # Index final BAM

        samtools view -h -f 0x2 {input.bam} \
        | samtools sort -n -@ 4 -O BAM -T {params.temp_prefix} -o - \
        | dupsifter \
            --add-mate-tags \
            --stats-output {log} \
            {input.fasta} - \
        | samtools sort -@ 8 -o {output.bam}
        samtools index -@ 8 {output.bam}
        """
rule emseq_methyldackel_dedup:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{data_dir}/analysis/emseq/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam",
        fasta = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    output:
        bed = f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        f"{log_dir}/{{library_id}}_{{ref_name}}_{{align_method}}_emseq_methyldackel_dedup.log",
    params:
        out_prefix = f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel",
        threads = 36,
    shell:
        """
        MethylDackel extract \
        -@ {params.threads} \
        --methylKit \
        {input.fasta} \
        {input.bam} \
        -o {params.out_prefix} > {log} 2>&1
        """
rule make_single_methylkit_amp_obj:
    conda:
        "../config/methylkit-conda-env.yaml",
    input:
        amp = f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.{{align_method}}_single_methylkit_amp.log",
    output:
        txt = f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.txt",
        bgz = f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.txt.bgz",
        tbi = f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.txt.bgz.tbi",
    params:
        Rscript = f"{emseq_script_dir}/make_single_amp_methylkit_obj.R",
        out_dir = f"{data_dir}/analysis/emseq/dmr/tabix",
        mincov = emseq_mincov,
        build = emseq_build,
        treatment = 1,
    shell:
        """
        Rscript {params.Rscript} \
          --amp_file {input.amp} \
          --library_id "{wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method}" \
          --mincov {params.mincov} \
          --out_dir {params.out_dir} \
          --treatment {params.treatment} \
          --build {params.build} \
          &>> {log}
        """
rule emseq_align_bwameth_spike:
    conda:
        "../config/emseq-conda-env.yaml"
    input:
        r1 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        ref = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa"
    output:
        bam = f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam"
    threads: 48
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.emseq_align_bwameth_spike.log"
    params:
        temp_prefix = lambda wildcards: f"{data_dir}/tmp/{wildcards.library_id}.{wildcards.ref_name}"
    shell:
        """
        bwameth.py --threads {threads} \
            --reference {input.ref} \
            {input.r1} {input.r2} 2> {log} | \
        samtools view -u -F 4 - | \
        samtools sort -@ {threads} -T {params.temp_prefix} -o {output.bam}
        """
rule emseq_methyldackel_spike:
    conda:
        "../config/emseq-conda-env.yaml"
    input:
        bam = f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsorted.bam",
        fasta = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    output:
        bed = f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        f"{log_dir}/{{library_id}}_{{ref_name}}_{{align_method}}_methyldackel.log",
    params:
        out_prefix = f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel",
    shell:
        """
        MethylDackel extract \
        --methylKit \
        {input.fasta} \
        {input.bam} \
        -o {params.out_prefix} > {log} 2>&1
        """
