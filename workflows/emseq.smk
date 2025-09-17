############################
###   EM-Seq Snakefile   ###
############################

#########1#########2#########3#########4#########5#########6#########7#########8
#
# A snakefile for basic processing of EM-seq sequencing data

# Fallbacks for lint-only execution
try:
    ENV_EMSEQ
except NameError:
    ENV_EMSEQ   = "envs/emseq.yml"
    D_EMSEQ     = "emseq"
    D_REF       = "ref"
    D_LOGS      = "logs"
    D_BENCHMARK = "benchmark"

emseq_mincov = 2
emseq_build = "hg38"

rule emseq_fastp:
    message: "EM-seq fastp for {wildcards.library_id}"
    conda: ENV_EMSEQ
    threads: 8
    params:
        extra = config.get("fastp", {}).get("extra", ""),
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    output:
        failed = f"{D_EMSEQ}/fastqs/{{library_id}}.failed.fastq.gz",
        html = f"{D_EMSEQ}/qc/{{library_id}}_emseq_fastp.html",
        json = f"{D_EMSEQ}/qc/{{library_id}}_emseq_fastp.json",
        r1     = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2     = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
    log:
        cmd  = f"{D_LOGS}/{{library_id}}_emseq_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_emseq_fastp.tsv"
    shell:
        r"""

        # Logging and console output
        exec &>> "{log.cmd}"
        echo "[fastp] $(date) lib={wildcards.library_id} threads={threads}"

        # Main
        fastp \
          --detect_adapter_for_pe \
          --disable_quality_filtering \
          --in1 "{input.r1}" --in2 "{input.r2}" \
          --out1 "{output.r1}" --out2 "{output.r2}" \
          --failed_out "{output.failed}" \
          --json "{output.json}" --html "{output.html}" \
          --thread {threads} \
          {params.extra}
        """

rule emseq_fastqc:
    message: "EM-seq FastQC for {wildcards.library_id} {wildcards.processing} {wildcards.read}"
    conda: ENV_EMSEQ
    threads: 4
    resources:
        concurrency = 25
    input:
        fq = f"{D_EMSEQ}/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
    output:
        html = f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.html",
        zip  = f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{processing}}_{{read}}_emseq_fastqc.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{processing}}_{{read}}_emseq_fastqc.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[fastqc] $(date) lib={wildcards.library_id} proc={wildcards.processing} read={wildcards.read} threads={threads}"

        fastqc \
          --outdir "$(dirname "{output.html}")" \
          --quiet \
          --threads {threads} \
          "{input.fq}"
        """

rule emseq_mosdepth:
    message: "EM-seq mosdepth for {wildcards.library_id} {wildcards.ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 8
    resources:
        concurrency = 20
    params:
        script       = f"{D_EMSEQ}/scripts/emseq_mosdepth.sh",
        quant_levels = config.get("mosdepth-quant-levels", ""),
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    output:
        summary       = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.summary.txt",
        global_dist   = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.global.dist.txt",
        region_dist   = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.region.dist.txt",
        regions       = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.regions.bed.gz",
        regions_idx   = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.regions.bed.gz.csi",
        quantized     = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.quantized.bed.gz",
        quantized_idx = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.quantized.bed.gz.csi",
        thresholds    = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.thresholds.bed.gz",
        thresholds_idx= f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.thresholds.bed.gz.csi",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mosdepth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mosdepth.tsv"
    shell:
        r"""

        exec &>> "{log.cmd}"
        echo "[mosdepth] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} aln={wildcards.align_method} threads={threads}"

        "{params.script}" \
          "{input.bam}" \
          "$(dirname "{output.summary}")" \
          "{wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method}" \
          '{params.quant_levels}' \
          {threads}
        """

rule emseq_mbias:
    message: "EM-seq MethylDackel mbias for {wildcards.library_id} {wildcards.ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 10
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        fasta = f"{D_REF}/{{align_method}}/{{ref_name}}/{{ref_name}}.fa",
    output:
        txt = f"{D_EMSEQ}/qc/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mbias.txt",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mbias.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_mbias.tsv"
    shell:
        r"""

        exec &>> "{log.cmd}"
        echo "[mbias] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} aln={wildcards.align_method} threads={threads}"

        MethylDackel mbias \
          -@ {threads} \
          --noSVG \
          "{input.fasta}" "{input.bam}" > "{output.txt}"
        """

rule emseq_dedup:
    message: "EM-seq dedup (dupsifter) for {wildcards.library_id} {wildcards.ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 8
    resources:
        concurrency = 25
    params:
        temp_prefix = lambda wc: f"{D_REF}/tmp/{wc.library_id}.{wc.ref_name}.{wc.align_method}.coorsort"
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{D_REF}/{{align_method}}/{{ref_name}}/{{ref_name}}.fa",
    output:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_dedup.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_dedup.tsv"
    shell:
        r"""

        mkdir -p "$(dirname "{params.temp_prefix}")"

        # append all stdout+stderr to the command log
        exec >> "{log.cmd}" 2>&1    # use this if you want POSIX; otherwise: exec &>> "{log.cmd}"

        echo "[dedup] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} aln={wildcards.align_method} threads={threads}"

        # Clean previous temp chunks (unquoted glob on purpose)
        rm -f {params.temp_prefix}.tmp.*

        # Proper pairs -> name-sort (BAM) -> dupsifter -> coord-sort -> index
        samtools view -bh -f 0x2 "{input.bam}" \
        | samtools sort -n -@ {threads} -O BAM -T "{params.temp_prefix}.tmp" -o - \
        | dupsifter \
        --add-mate-tags \
        --stats-output "{log.cmd}.dupsifter.tsv" \
        "{input.fasta}" - \
        | samtools sort -@ {threads} -O BAM -T "{params.temp_prefix}.tmp" -o "{output.bam}"

        samtools index -@ {threads} "{output.bam}"
        """

rule emseq_align_bwameth_spike:
    message: "EM-seq bwameth (spike) for {wildcards.library_id} {wildcards.ref_name}"
    conda: ENV_EMSEQ
    threads: 48
    params:
        temp_prefix = lambda wc: f"{D_EMSEQ}/tmp/{wc.library_id}.{wc.ref_name}"
    input:
        r1  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        ref = f"{D_REF}/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    output:
        bam = f"{D_EMSEQ}/spike/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}_emseq_align_bwameth_spike.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}_emseq_align_bwameth_spike.tsv"
    shell:
        r"""

        exec &>> "{log.cmd}"
        echo "[bwameth-spike] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        bwameth.py --threads {threads} \
          --reference "{input.ref}" \
          "{input.r1}" "{input.r2}" \
        | samtools view -u -F 4 - \
        | samtools sort -@ {threads} -T "{params.temp_prefix}" -o "{output.bam}"
        """

rule emseq_methyldackel_spike:
    message: "EM-seq MethylDackel (spike) for {wildcards.library_id} {wildcards.ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 8
    input:
        bam   = f"{D_EMSEQ}/spike/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{D_REF}/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    output:
        bed = f"{D_EMSEQ}/spike/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_methyldackel_spike.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_methyldackel_spike.tsv"
    params:
        out_prefix = lambda wc, input, output: output.bed.rsplit("_CpG.methylKit", 1)[0]
    shell:
        r"""
        mkdir -p "$(dirname "{params.out_prefix}")"


        exec &>> "{log.cmd}"
        echo "[methyldackel-spike] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} aln={wildcards.align_method} threads={threads}"

        MethylDackel extract \
          -@ {threads} \
          --methylKit \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """

# -------- Index (bwa-meth) --------
rule bwa_meth_index:
    message: "bwa-meth index for {wildcards.name}"
    conda: ENV_EMSEQ
    threads: 8
    input:
        lambda wc: f"{D_INPUTS}/{config['emseq_ref_assemblies'][wc.name]['input']}"
    output:
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t",
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.0123",
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.amb",
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.ann",
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.bwt.2bit.64",
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa.bwameth.c2t.pac",
        f"{D_REF}/bwa_meth/{{name}}/{{name}}.fa",
    params:
        fasta_target = lambda wc: f"{D_REF}/bwa_meth/{wc.name}/{wc.name}.fa"
    log:
        cmd = f"{D_LOGS}/{{name}}_bwa_meth_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{name}}_bwa_meth_index.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[bwa-meth index] $(date) name={wildcards.name} threads={threads}"

        zcat "{input}" > "{params.fasta_target}"
        samtools faidx -@ {threads} "{params.fasta_target}"
        bwameth.py index-mem2 "{params.fasta_target}"
        """

# -------- Align (bwa-meth) --------
rule emseq_align_bwameth:
    message: "EM-seq bwameth for {wildcards.library_id} {wildcards.ref_name}"
    conda: ENV_EMSEQ
    threads: 16
    params:
        temp_prefix = lambda wc: f"{D_EMSEQ}/tmp/{wc.library_id}.{wc.ref_name}"
    input:
        r1  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        ref = f"{D_REF}/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
        c2t = f"{D_REF}/bwa_meth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t",
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.bwa_meth.coorsort.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}_bwameth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}_bwameth.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[bwameth] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        mkdir -p "$(dirname "{params.temp_prefix}")"

        bwameth.py --threads {threads} \
          --reference "{input.ref}" \
          "{input.r1}" "{input.r2}" \
        | samtools view -u - \
        | samtools sort -@ {threads} -T "{params.temp_prefix}" -o "{output.bam}"
        """

# -------- Call methylation (MethylDackel) --------
rule emseq_methyldackel:
    message: "EM-seq MethylDackel for {wildcards.library_id} {wildcards.ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 20
    resources:
        concurrency = 25
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.deduped.bam",
        fasta = f"{D_REF}/bwa_meth/{{ref_name}}/{{ref_name}}.fa",
    output:
        bed = f"{D_EMSEQ}/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    params:
        out_prefix = lambda wc, input, output: output.bed.rsplit("_CpG.methylKit", 1)[0]
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_methyldackel_dedup.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}.{{align_method}}_emseq_methyldackel.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[methyldackel] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} aln={wildcards.align_method} threads={threads}"

        MethylDackel extract \
          -@ {threads} \
          --methylKit \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """

# -------- Build single methylKit (R) --------
rule make_single_methylkit_amp_obj:
    message: "Build tabix-backed methylKit object for {wildcards.library_id} {wildcards.ref_name} {wildcards.align_method}"
    conda: "../config/methylkit-conda-env.yaml"
    threads: 1
    input:
        f"{D_EMSEQ}/meth/{{library_id}}.{{ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    output:
        bgz = f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.methyldackel.txt.bgz",
        tbi = f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{ref_name}}.{{align_method}}.methyldackel.txt.bgz.tbi",
    params:
        Rscript   = f"{R_EMSEQ}/scripts/make_single_amp_methylkit_obj.R",
        mincov    = emseq_mincov,
        build     = emseq_build,
        treatment = 1,
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}.{{align_method}}_single_methylkit_amp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}.{{align_method}}_single_methylkit_amp.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[methylKit-amp] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} aln={wildcards.align_method}"

        Rscript "{params.Rscript}" \
          --amp_file "{input}" \
          --library_id "{wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method}.methyldackel" \
          --mincov {params.mincov} \
          --out_dir "$(dirname "{output.bgz}")" \
          --treatment {params.treatment} \
          --build {params.build}
        """
