############################
###   EM-Seq Snakefile   ###
############################

#########1#########2#########3#########4#########5#########6#########7#########8
#
# A snakefile for basic processing of EM-seq sequencing data
rule emseq_align_bwameth_spike:
    message: "EM-seq bwameth (spike) for {wildcards.library_id} {wildcards.emseq_ref_name}"
    conda: ENV_EMSEQ
    threads: 48
    params:
        temp_prefix = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.emseq_ref_name}"
    resources:
        concurrency=50
    input:
        r1  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        ref = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    output:
        bam = f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.bwa_meth.coorsort.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}_emseq_align_bwameth_spike.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}_emseq_align_bwameth_spike.tsv"
    shell:
        r"""

        exec &>> "{log.cmd}"
        echo "[bwameth-spike] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} threads={threads}"

        bwameth.py --threads {threads} \
          --reference "{input.ref}" \
          "{input.r1}" "{input.r2}" \
        | samtools view -@ 8 -u -F 4 - \
        | samtools sort -@ 8 -T "{params.temp_prefix}" -o "{output.bam}"
        """
rule emseq_methyldackel_spike:
    message: "EM-seq MethylDackel (spike) for {wildcards.library_id} {wildcards.emseq_ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 8
    input:
        bam   = f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    output:
        bed = f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel_spike.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel_spike.tsv"
    params:
        out_prefix = lambda wc, input, output: output.bed.rsplit("_CpG.methylKit", 1)[0]
    shell:
        r"""
        mkdir -p "$(dirname "{params.out_prefix}")"


        exec &>> "{log.cmd}"
        echo "[methyldackel-spike] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"

        MethylDackel extract \
          -@ {threads} \
          --methylKit \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """
# -------- Index (bwa-meth) --------
rule bwa_meth_index:
    message: "bwa-meth index for {wildcards.emseq_ref_name}"
    conda: ENV_EMSEQ
    input:
        lambda wc: f"{D_INPUTS}/{config['emseq_ref_assemblies'][wc. emseq_ref_name]['input']}"
    output:
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.0123",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.amb",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.ann",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.bwt.2bit.64",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.pac",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
    params:
        fasta_target = lambda wc: f"{D_REF}/bwa_meth/{wc.emseq_ref_name}/{wc.emseq_ref_name}.fa"
    log:
        cmd = f"{D_LOGS}/{{emseq_ref_name}}_bwa_meth_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{emseq_ref_name}}_bwa_meth_index.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[bwa-meth index] $(date) emseq_ref_name={wildcards.emseq_ref_name} threads={threads}"


        if file -b "{input}" | grep -qi gzip; then
            zcat "{input}" > "{params.fasta_target}"
        else
            cat "{input}" > "{params.fasta_target}"
        fi

        samtools faidx "{params.fasta_target}"
        bwameth.py index-mem2 "{params.fasta_target}"
        """
# -------- Align (bwa-meth) --------
rule emseq_align_bwameth:
    message: "EM-seq bwameth for {wildcards.library_id} {wildcards.emseq_ref_name}"
    conda: ENV_EMSEQ
    threads: min(workflow.cores, 48)
    params:
        temp_prefix = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.emseq_ref_name}",
    resources:
        concurrency = 50,
    input:
        r1  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        ref = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        c2t = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t",
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.bwa_meth.coorsort.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}_bwameth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}_bwameth.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[bwameth] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} threads={threads}"

        mkdir -p "$(dirname "{params.temp_prefix}")"

        bwameth.py --threads {threads} \
          --reference "{input.ref}" \
          "{input.r1}" "{input.r2}" \
        | samtools view -u - \
        | samtools sort -@ 8 -T "{params.temp_prefix}" -o "{output.bam}"
        """
rule biscuit_index:
    message: "biscuit index for {wildcards.emseq_ref_name}"
    conda: ENV_EMSEQ
    input:
        lambda wc: f"{D_INPUTS}/{config['emseq_ref_assemblies'][wc.emseq_ref_name]['input']}"
    output:
        fasta = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        fai = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
        index = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.par.sa",
        biscuit_index_done = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.biscuit.index.done"
    log:
        cmd = f"{D_LOGS}/{{emseq_ref_name}}_biscuit_index.log"
    shell:
        r"""
        exec &>> "{log.cmd}"
        mkdir -p "$(dirname "{output.fasta}")"
        zcat "{input}" > "{output.fasta}"
        samtools faidx "{output.fasta}"
        biscuit index "{output.fasta}"
        touch "{output.biscuit_index_done}"
        """
rule emseq_align_biscuit:
    message: "Biscuit alignment for {wildcards.emseq_ref_name}"
    conda: ENV_EMSEQ
    threads: min(workflow.cores, 48)
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        fasta = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        index = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.par.sa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.biscuit.coorsort.log",
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.biscuit.coorsort.bam",
    params:
        tmp_dir = D_DATA,
    resources:
        concurrency = 100
    shell:
        r"""
        exec &>> "{log.cmd}"
        mkdir -p "{params.tmp_dir}/tmp"
        biscuit align -@ {threads} "{input.fasta}" "{input.r1}" "{input.r2}" \
        | samtools sort -@ {threads} -m 2G \
            -T "{params.tmp_dir}/tmp/{wildcards.library_id}_sorttmp" \
            -o "{output.bam}"
        """
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
    message: "EM-seq mosdepth for {wildcards.library_id} {wildcards.emseq_ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 8
    resources:
        concurrency = 20
    params:
        script       = f"{R_EMSEQ}/scripts/emseq_mosdepth.sh",
        quant_levels = config.get("mosdepth-quant-levels", ""),
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    output:
        summary       = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.summary.txt",
        global_dist   = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.global.dist.txt",
        region_dist   = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.region.dist.txt",
        regions       = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.regions.bed.gz",
        regions_idx   = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.regions.bed.gz.csi",
        quantized     = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.quantized.bed.gz",
        quantized_idx = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.quantized.bed.gz.csi",
        thresholds    = f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.thresholds.bed.gz",
        thresholds_idx= f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.thresholds.bed.gz.csi",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mosdepth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mosdepth.tsv"
    shell:
        r"""

        exec &>> "{log.cmd}"
        echo "[mosdepth] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"

        "{params.script}" \
        "{input.bam}" \
        "$(dirname "{output.summary}")" \
        "{wildcards.library_id}.{wildcards.emseq_ref_name}.{wildcards.align_method}" \
        '{params.quant_levels}' \
        {threads}
        """
rule emseq_mbias:
    message: "EM-seq MethylDackel mbias for {wildcards.library_id} {wildcards.emseq_ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 10
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        fasta = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    output:
        txt = f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.txt",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.tsv"
    shell:
        r"""

        exec &>> "{log.cmd}"
        echo "[mbias] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"

        MethylDackel mbias \
          -@ {threads} \
          --noSVG \
          "{input.fasta}" "{input.bam}" > "{output.txt}"
        """
rule emseq_dedup:
    message: "EM-seq dedup (dupsifter) for {wildcards.library_id} {wildcards.emseq_ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 8
    resources:
        concurrency = 25
    params:
        temp_prefix = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.emseq_ref_name}.{wc.align_method}.coorsort"
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    output:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_dedup.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_dedup.tsv"
    shell:
        r"""

        mkdir -p "$(dirname "{params.temp_prefix}")"

        # append all stdout+stderr to the command log
        exec >> "{log.cmd}" 2>&1    # use this if you want POSIX; otherwise: exec &>> "{log.cmd}"

        echo "[dedup] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"

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
rule emseq_filter_bam:
    message: "Filtering deduped BAM for {wildcards.library_id}"
    conda: ENV_EMSEQ
    threads: 16
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        bai   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
        keep_bed    = KEEP_BED,     # from preamble, bedtools-sorted with this FAI
        exclude_bed = EXCL_BED,     # from preamble, bedtools-sorted with this FAI
        fai   = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",  # produced by your index rule
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
        bai = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam.bai",
    shell: r"""
        samtools view -@ {threads} -u -f 2 -q 30 -F 2816 "{input.bam}" \
          | bedtools intersect -sorted -g "{input.fai}" -a stdin -b "{input.exclude_bed}" -v -ubam \
          | bedtools intersect -sorted -g "{input.fai}" -a stdin -b "{input.keep_bed}"    -ubam \
          > "{output.bam}"
        samtools index -@ {threads} "{output.bam}" "{output.bai}"
    """
rule emseq_samtools_stats:
    message: "Samtools stats + flagstat for {wildcards.library_id}.{wildcards.emseq_ref_name}"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
    log:
        f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.stats.log",
    output:
        stats   = f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.stats.txt",
        flagstat = f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.flagstat.txt",
    params:
        threads = 8,
    threads: 8
    resources:
        concurrency=40
    shell:
        r"""
        samtools stats -@ {threads} {input.bam} > {output.stats} 2>> {log}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} 2>> {log}
        """
# -------- Call methylation (MethylDackel) --------
rule emseq_methyldackel:
    message: "EM-seq MethylDackel for {wildcards.library_id} {wildcards.emseq_ref_name} {wildcards.align_method}"
    conda: ENV_EMSEQ
    threads: 20
    resources:
        concurrency = 25
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
        fasta = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    output:
        bed = f"{D_EMSEQ}/meth/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    params:
        out_prefix = lambda wc, input, output: output.bed.rsplit("_CpG.methylKit", 1)[0]
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel_dedup.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[methyldackel] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"

        MethylDackel extract \
          -@ {threads} \
          --methylKit \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """
# -------- Build single methylKit (R) --------
rule make_single_methylkit_amp_obj:
    message: "Build tabix-backed methylKit object for {wildcards.library_id} {wildcards.emseq_ref_name} {wildcards.align_method}"
    conda: ENV_METHYLKIT
    threads: 1
    input:
        f"{D_EMSEQ}/meth/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    output:
        bgz = f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
        tbi = f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz.tbi",
    params:
        Rscript   = f"{R_EMSEQ}/scripts/make_single_amp_methylkit_obj.R",
        mincov    = emseq_mincov,
        build     = lambda wc: {wc.emseq_ref_name},
        treatment = 1,
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_single_methylkit_amp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_single_methylkit_amp.tsv"
    shell:
        r"""
        exec &>> "{log.cmd}"
        echo "[methylKit-amp] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method}"

        Rscript "{params.Rscript}" \
          --amp_file "{input}" \
          --library_id "{wildcards.library_id}.{wildcards.emseq_ref_name}.{wildcards.align_method}.methyldackel" \
          --mincov {params.mincov} \
          --out_dir "$(dirname "{output.bgz}")" \
          --treatment {params.treatment} \
          --build {params.build}
        """
rule emseq_multiqc:
    message: f"MultiQC report for EM-seq (n={len(emseq_library_ids)} libraries)"
    conda: ENV_EMSEQ
    input:
        fastqc = expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
            library_id=emseq_library_ids,
            processing=["raw","trimmed"],
            read=["R1","R2"],
            allow_missing=True,
        ),
        fastp_html = expand(
            f"{D_EMSEQ}/qc/{{library_id}}_emseq_fastp.html",
            library_id=emseq_library_ids,
            allow_missing=True,
        ),
        fastp_json = expand(
            f"{D_EMSEQ}/qc/{{library_id}}_emseq_fastp.json",
            library_id=emseq_library_ids,
            allow_missing=True,
        ),
        mbias = expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["biscuit","bwa_meth"],
            allow_missing=True,
        ),
        mosdepth_summary = expand(
            f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.summary.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["biscuit","bwa_meth"],
            allow_missing=True,
        ),
        mosdepth_dists = expand(
            f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.{{dist}}.dist.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["biscuit","bwa_meth"],
            dist=["global","region"],
            allow_missing=True,
        ),
        samstats = expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.{{stat}}.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["biscuit","bwa_meth"],
            stat = ["stats","flagstat"]),
    log:
        f"{D_LOGS}/emseq_multiqc.log"
    output:
        html = f"{D_EMSEQ}/qc/multiqc.html",
        data = directory(f"{D_EMSEQ}/qc/multiqc_data"),  # <-- match MultiQC
    params:
        extra = "--force",
    threads: 4
    resources:
        concurrency=20
    shell:
        r"""
        mkdir -p $(dirname {output.html})

        multiqc \
            {input} \
            {params.extra} \
            --outdir $(dirname {output.html}) \
            --filename $(basename {output.html}) \
            &> {log}
        """
rule make_methylkit_unite_db:
    conda: ENV_METHYLKIT
    threads: 32
    input:
        mkit_lib_db=lambda wc: expand(
            f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
            library_id=meth_map[wc.experiment]["libs"],
            emseq_ref_name=meth_map[wc.experiment]["emseq_ref_name"],
            align_method=meth_map[wc.experiment]["align_method"],
        ),
    output:
        mbase=f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
    params:
        library_id=lambda wc: " ".join(meth_map[wc.experiment]["libs"]),
        treatment_list=lambda wc: " ".join(map(str, meth_map[wc.experiment]["tx"])),
        script=f"{R_EMSEQ}/scripts/make_methylkit_unite_db.R",
        mincov=lambda wc: meth_map[wc.experiment]["mincov"],
        min_per_group=lambda wc: meth_map[wc.experiment]["mingroup"],
        chunk_size=lambda wc: meth_map[wc.experiment]["chunksize"],
    log:
        f"{D_LOGS}/{{experiment}}_make_methylkit_unite_db.log",
    shell:
        r"""
        # destructive pre-clean (script also unlinks)
        rm -f {output.mbase}*

        Rscript {params.script} \
          --lib_db_list "{input.mkit_lib_db}" \
          --lib_id_list "{params.library_id}" \
          --treatment_list "{params.treatment_list}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.mbase}")" \
          --suffix {wildcards.experiment} \
          --assembly "hg38" \
          --mincov {params.mincov} \
          --min_per_group {params.min_per_group} \
          --chunk_size {params.chunk_size} \
          > {log} 2>&1
        """
rule make_methylkit_diff_db:
    conda: ENV_METHYLKIT
    threads: 32
    input:
        mbase = f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
    output:
        mdiff = f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
    params:
        script     = f"{R_EMSEQ}/scripts/make_methylkit_diff_db.R",
        chunk_size = lambda wc: meth_map[wc.experiment]["chunksize"],
    log:
        f"{D_LOGS}/{{experiment}}_make_methylkit_diff_db.log",
    shell:
        r"""
        # destructive pre-clean
        rm -f {output.mdiff} {output.mdiff}.tbi

        Rscript {params.script} \
          --mbase "{input.mbase}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.mdiff}")" \
          --suffix {wildcards.experiment} \
          --chunk_size {params.chunk_size} \
          > {log} 2>&1
        """
rule make_methylkit_diff_db_tiled:
    conda: ENV_METHYLKIT
    input:
        mkit_lib_db = lambda wildcards: expand(f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
                                               library_id = meth_map[wildcards.experiment]['libs'],
                                               emseq_ref_name = meth_map[wildcards.experiment]['emseq_ref_name'],
                                               align_method = meth_map[wildcards.experiment]['align_method']),
    log:
        f"{D_LOGS}/{{experiment}}_make_methylkit_diff_tiled_db.log",
    output:
        unite = f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}_tiled.txt.bgz",
        diff = f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}_tiled.txt.bgz",
    params:
        library_id = lambda wildcards: " ".join(meth_map[wildcards.experiment]['libs']),
        treatment_list = lambda wc: " ".join(map(str, meth_map[wc.experiment]['tx'])),
        script = f"{R_EMSEQ}/scripts/make_methylkit_diff_tiled_db.R",
        mincov=lambda wc: meth_map[wc.experiment]["mincov"],
        chunk_size=lambda wc: meth_map[wc.experiment]["chunksize"],
        min_per_group=lambda wc: meth_map[wc.experiment]["mingroup"],
        win_size = lambda wc: meth_map[wc.experiment]["win_size"],
    shell:
        """
        # remove old outputs if they exist
        rm -f {output.unite}* {output.diff}*

        Rscript {params.script} \
        --lib_db_list "{input.mkit_lib_db}" \
        --lib_id_list "{params.library_id}" \
        --treatment_list "{params.treatment_list}" \
        --cores {threads} \
        --out_dir "$(dirname "{output.unite}")" \
        --suffix {wildcards.experiment} \
        --assembly "hg38" \
        --mincov {params.mincov} \
        --win_size {params.win_size} \
        --min_per_group {params.min_per_group} \
        --chunk_size {params.chunk_size} \
        > {log} 2>&1
        """
