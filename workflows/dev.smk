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
        threads = 8,
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
        concurrency = 25,
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
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
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
    threads: 8,
    resources:
        concurrency = 20,
    shell:
        """
        {params.script} \
        {input.bam} \
        {params.out_dir} \
        {wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method} \
        '{params.quant_levels}' \
        {threads}
        """
rule emseq_mosdepth_agg_plot:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        thresholds = lambda wildcards: expand(
            f"{data_dir}/qc/mosdepth_{{library_id}}.{mosdepth_map[wildcards.experiment]['ref_name']}."
            f"{mosdepth_map[wildcards.experiment]['align_method']}.thresholds.bed.gz",
            library_id=mosdepth_map[wildcards.experiment]['library_ids']
        ),
        regions = lambda wildcards: expand(
            f"{data_dir}/qc/mosdepth_{{library_id}}.{mosdepth_map[wildcards.experiment]['ref_name']}."
            f"{mosdepth_map[wildcards.experiment]['align_method']}.regions.bed.gz",
            library_id=mosdepth_map[wildcards.experiment]['library_ids']
        )
    output:
        pdf = f"{data_dir}/qc/{{experiment}}.emseq_mosdepth_agg_plot.pdf",
        tsv = f"{data_dir}/qc/{{experiment}}.emseq_mosdepth_agg.tsv",
    params:
        script = f"{emseq_script_dir}/emseq_mosdepth_agg_plot.R",
        library_list = lambda wildcards: " ".join(mosdepth_map[wildcards.experiment]['library_ids']),
        threshold_list = lambda wildcards, input: " ".join(input.thresholds),
        regions_list = lambda wildcards, input: " ".join(input.regions),
    shell:
        """
        Rscript {params.script} \
        --threshold_list "{params.threshold_list}" \
        --regions_list "{params.regions_list}" \
        --library_list "{params.library_list}" \
        --output_pdf {output.pdf} \
        --output_tsv {output.tsv}
        """
rule emseq_mbias:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsorted.deduped.bam",
        fasta = f"{data_dir}/ref/{{align_method}}/{{ref_name}}/{{ref_name}}.fa",
    log:
        f"{data_dir}/logs/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mbias.log",
    output:
        f"{data_dir}/qc/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mbias.txt",
    shell:
        r"""
        MethylDackel mbias \
        -@ 10 \
        --noSVG \
        {input.fasta} {input.bam} > {ouput}
        """
rule emseq_dedup:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{data_dir}/ref/{{align_method}}/{{ref_name}}/{{ref_name}}.fa",
    log:
        f"{data_dir}/logs/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_dedup.log",
    output:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    params:
        temp_prefix = lambda wildcards: f"{data_dir}/tmp/{wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method}.coorsort",
    resources:
        concurrency = 25,
    shell:
        r"""
        # Clean up any existing temp BAM chunks from previous failed sort
        rm -f {params.temp_prefix}.*

        # Extract only properly paired reads (-f 0x2), include header (-h)
        # Name-sort BAM using samtools; output to stdout
        # Deduplicate using dupsifter, using reference FASTA and streaming input from stdin
        # Coordinate-sort deduplicated BAM for downstream tools
        # Index final BAM

        samtools view -h -f 0x2 {input.bam} \
        | samtools sort -n -@ 8 -O BAM -T {params.temp_prefix} -o - \
        | dupsifter \
            --add-mate-tags \
            --stats-output {log} \
            {input.fasta} - \
        | samtools sort -@ 8 -o {output.bam}
        samtools index -@ 8 {output.bam}
        """
rule emseq_align_bwameth_spike:
    conda:
        "../config/emseq-conda-env.yaml"
    input:
        r1 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        ref = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa"
    output:
        bam = f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam"
    threads: 48
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.emseq_align_bwameth_spike.log"
    params:
        temp_prefix = lambda wildcards: f"{data_dir}/tmp/{wildcards.library_id}.{wildcards.ref_name}"
    shell:
        """
        mkdir -p $(dirname {params.temp_prefix})
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
        bam = f"{data_dir}/analysis/emseq/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam",
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
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam"
    threads: 16
    params:
        temp_prefix = lambda wildcards: f"{data_dir}/tmp/{wildcards.library_id}.{wildcards.ref_name}",
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.bwameth.log"
    shell:
        """
        mkdir -p $(dirname {params.temp_prefix})
        bwameth.py --threads {threads} \
            --reference {input.ref} \
            {input.r1} {input.r2} 2> {log} | \
        samtools view -u - | \
        samtools sort -@ {threads} -T {params.temp_prefix} -o {output.bam}
        """
rule emseq_methyldackel:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        fasta = f"{data_dir}/ref/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    output:
        f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_methyldackel_dedup.log",
    params:
        out_prefix = f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel",
        threads = 20,
    resources:
        concurrency = 25,
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
        f"{data_dir}/analysis/emseq/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        f"{log_dir}/{{library_id}}.{{ref_name}}.{{align_method}}_single_methylkit_amp.log",
    output:
        txt = f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.methyldackel.txt",
        bgz = f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.methyldackel.txt.bgz",
        tbi = f"{data_dir}/analysis/emseq/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.methyldackel.txt.bgz.tbi",
    params:
        Rscript = f"{emseq_script_dir}/make_single_amp_methylkit_obj.R",
        out_dir = f"{data_dir}/analysis/emseq/dmr/tabix",
        mincov = emseq_mincov,
        build = emseq_build,
        treatment = 1,
    shell:
        """
        Rscript {params.Rscript} \
          --amp_file {input} \
          --library_id "{wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method}.methyldackel" \
          --mincov {params.mincov} \
          --out_dir {params.out_dir} \
          --treatment {params.treatment} \
          --build {params.build} \
          &>> {log}
        """
rule emseq_biscuit_index:
    conda:
        "../config/emseq-conda-env.yaml"
    input:
        lambda wildcards: f"{data_dir}/inputs/{config['emseq_ref_assemblies'][wildcards.name]['input']}"
    output:
        fasta = f"{data_dir}/ref/biscuit/{{name}}/{{name}}.fa",
        fai = f"{data_dir}/ref/biscuit/{{name}}/{{name}}.fa.fai",
        biscuit_index_done = f"{data_dir}/ref/biscuit/{{name}}/{{name}}.fa.biscuit.index.done"
    log:
        f"{log_dir}/{{name}}_biscuit_index.log"
    shell:
        """
        mkdir -p $(dirname {output.fasta}) && \
        zcat {input} > {output.fasta} && \
        samtools faidx {output.fasta} && \
        biscuit index {output.fasta} > {log} 2>&1 && \
        touch {output.biscuit_index_done}
        """
rule emseq_align_biscuit:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        r1 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/emseq/fastqs/{{library_id}}_trimmed_R2.fastq.gz",
        fasta = f"{data_dir}/ref/biscuit/{{ref_name}}/{{ref_name}}.fa",
        index = f"{data_dir}/ref/biscuit/{{ref_name}}/{{ref_name}}.fa.par.sa",
    log:
        cmd = f"{data_dir}/logs/{{library_id}}.{{ref_name}}.biscuit.coorsorted.bam",
    output:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.biscuit.coorsorted.bam",
    params:
        threads = 80,
    resources:
        concurrency=100
    shell:
        """
        mkdir -p {data_dir}/tmp && \
        biscuit align \
        -@ {params.threads} \
        -biscuit-ref {input.fasta} \
        {input.r1} {input.r2} \
        | samtools sort \
        -@ 8 \
        -m 2G \
        -T {data_dir}/tmp/{wildcards.library_id}_sorttmp \
        -o {output.bam} &>> {log}
        """
rule emseq_biscuit_pileup:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.{{ref_name}}.biscuit.coorsort.deduped.bam",
        fasta = f"{data_dir}/ref/biscuit/{{ref_name}}/{{ref_name}}.fa",
    log:
        f"{data_dir}/logs/{{library_id}}.{{ref_name}}.biscuit_emseq_pileup.log",
    output:
        vcf = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}.biscuit_pileup.vcf.gz",
        tsv = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}.biscuit_pileup.vcf_meth_average.tsv",
    params:
        out_base = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}.biscuit_pileup.vcf",
    shell:
        """
        biscuit pileup \
        -@ 20 \
        -o {params.out_base} \
        {input.fasta} {input.bam} \
        && bgzip -@ 8 {params.out_base}
        """
rule emseq_biscuit_post_pileup:
    conda:
        "../config/emseq-conda-env.yaml",
    input:
        vcf = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}.biscuit_pileup.vcf.gz",
    log:
        f"{data_dir}/logs/{{library_id}}.{{ref_name}}_emseq_biscuit_post_pileup.log",
    output:
        tbi = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}_pileup.vcf.gz.tbi",
        bed = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}_pileup.bed",
        bismark = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}_bismark_cov.bed",
    shell:
        """
        tabix -p vcf {input.vcf} \
        && biscuit vcf2bed \
	-t cg {input.vcf} > {output.bed} \
        && biscuit vcf2bed -c {input.vcf} > {output.bismark} &> {log}
        """
rule make_single_biscuit_methylkit_obj:
    conda:
        "../config/methylkit-conda-env.yaml",
    input:
        bismark = f"{data_dir}/analysis/emseq/pileup/{{library_id}}.{{ref_name}}_bismark_cov.bed",
    log:
        f"{data_dir}/logs/{{library_id}}.{{ref_name}}_make_single_biscuit_methylkit_obj.log",
    output:
        txt = f"{data_dir}/analysis/emseq/post-biscuit/{{library_id}}.{{ref_name}}_biscuit.txt",
        bgz = f"{data_dir}/analysis/emseq/post-biscuit/{{library_id}}.{{ref_name}}_biscuit.txt.bgz",
        tbi = f"{data_dir}/analysis/emseq/post-biscuit/{{library_id}}.{{ref_name}}_biscuit.txt.bgz.tbi",
    params:
        Rscript = f"{emseq_script_dir}/make_single_biscuit_methylkit_obj.R",
        out_dir = f"{data_dir}/analysis/emseq/post-biscuit",
    shell:
        """
        Rscript {params.Rscript} \
          --bismark_cov_bed {input.bismark} \
          --library_id {wildcards.library_id} \
          --out_dir {params.out_dir} \
          &> {log}
        """
