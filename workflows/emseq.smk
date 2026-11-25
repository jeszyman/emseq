############################
###   EM-Seq Snakefile   ###
############################

#########1#########2#########3#########4#########5#########6#########7#########8
#
# A snakefile for basic processing of EM-seq sequencing data
#
# Required wrapper-provided variables (must be in scope before include):
#   ENV_EMSEQ, ENV_METHYLKIT       — conda environment YAML paths
#   R_EMSEQ                        — repository root path
#   D_DATA, D_EMSEQ, D_REF, D_LOGS, D_BENCHMARK, D_INPUTS — directory layout
#   emseq_library_ids              — list of sample IDs
#   emseq_ref_names                — list of reference genome names to align against
#   spike_builds                   — list of spike-in reference names
#   spike_host_ref                 : sample reference that the controls are aligned together with
#   EMSEQ_REF_INPUTS               — dict mapping ref name to input FASTA filename
#   KEEP_BED, EXCL_BED             — region filter BED file paths
#   meth_map                       — experiment map for differential methylation
#   MOSDEPTH_QUANT_LEVELS          — coverage quantization thresholds
#   emseq_align_methods            — list of alignment methods (e.g. ["bwa_meth"] or ["bwa_meth","biscuit"])
#   EMSEQ_MINCOV                   — minimum coverage for methylKit
#   FASTP_EXTRA                    — additional fastp arguments

# Read pairs kept for dedup and filtering, applied with `samtools view -f 1 -F 12`:
# mates on the same chromosome, facing inward, fragment <= 1000 bp. This
# replaces the proper-pair flag (0x2), which bwa-meth sets from an insert-size
# distribution inferred per library; on cfDNA it stops at about 270-280 bp and
# a `-f 2` filter drops every dinucleosome-length fragment. In samtools
# expressions flag.reverse is 16 and flag.mreverse is 32 when set, so the
# strands are compared after negation, which turns each into 0 or 1.
EMSEQ_PAIR_FILTER = (
    "mrname == rname && !flag.reverse != !flag.mreverse"
    " && ((!flag.reverse && tlen > 0) || (flag.reverse && tlen < 0))"
    " && tlen >= -1000 && tlen <= 1000"
)

# Bases excluded from methylation calls (MethylDackel --nOT/--nOB): the first
# five sequenced bases of read 2, which start at the 3' end of the original
# strand that end repair fills in with unmethylated C. These are the NEB EM-seq
# pipeline's options; tests/unit/test_methyldackel_bounds.py pins what they
# exclude on constructed read pairs.
EMSEQ_CALL_BOUNDS = "--nOT 0,0,0,5 --nOB 0,0,5,0"

# Reference that holds the sample genome and every control
EMSEQ_SPIKE_REF = f"{spike_host_ref}_spike"
# ── Spike-in Controls ────────────────────────────────────────────────────────
# Align to pUC19 (methylated) and Lambda (unmethylated) for conversion QC.
# Duplicates are retained — every read contributes to conversion rate estimation.
rule emseq_align_bwameth_spike:
    message: "Align to the sample reference and the controls together, keeping the reads on a control"
    conda: ENV_EMSEQ
    input:
        r1   = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2   = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        ref  = f"{D_REF}/bwa_meth/{EMSEQ_SPIKE_REF}/{EMSEQ_SPIKE_REF}.fa",
        c2t  = f"{D_REF}/bwa_meth/{EMSEQ_SPIKE_REF}/{EMSEQ_SPIKE_REF}.fa.bwameth.c2t",
        fais = expand(f"{D_REF}/bwa_meth/{{spike}}/{{spike}}.fa.fai", spike=spike_builds),
    log:
        cmd = f"{D_LOGS}/{{library_id}}_emseq_align_bwameth_spike.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_emseq_align_bwameth_spike.tsv"
    params:
        temp_prefix = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{EMSEQ_SPIKE_REF}",
    threads: workflow.cores
    resources:
        concurrency = 100
    output:
        bam = temp(f"{D_EMSEQ}/spike/{{library_id}}.{EMSEQ_SPIKE_REF}.bwa_meth.coorsort.bam"),
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[bwameth-spike] $(date) lib={wildcards.library_id} threads={threads}"
        mkdir -p "$(dirname "{params.temp_prefix}")"
        # rname == "a" || rname == "b" over the control contigs
        on_control=$(cut -f1 {input.fais} | awk '{{printf "%srname == \\"%s\\"", (NR > 1 ? " || " : ""), $1}}')
        # -F 2308: primary alignments only. A pair whose mates lie on different
        # references is dropped, since a mate on the sample leaves the pair's
        # control origin open.
        bwameth.py --threads {threads} \
          --reference "{input.ref}" \
          "{input.r1}" "{input.r2}" \
        | samtools view -@ 8 -u -F 2308 -e "mrname == rname && ($on_control)" - \
        | samtools sort -@ 8 -T "{params.temp_prefix}" -o "{output.bam}"
        """
rule emseq_split_spike:
    message: "Write the alignments of one control, with only its contigs in the header"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/spike/{{library_id}}.{EMSEQ_SPIKE_REF}.bwa_meth.coorsort.bam",
        fai = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}_emseq_split_spike.log",
    wildcard_constraints:
        emseq_ref_name = "|".join(spike_builds),
    threads: 2
    output:
        bam = f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.bwa_meth.coorsort.bam",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[split-spike] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name}"
        # MethylDackel reads the control's own FASTA, so the header must hold no
        # other contig
        samtools view -h "{input.bam}" \
          | awk 'BEGIN {{FS = "\t"}} FNR == NR {{keep[$1]; next}} /^@SQ/ {{if (substr($2, 4) in keep) print; next}} /^@/ {{print; next}} $3 in keep' "{input.fai}" - \
          | samtools view -b -o "{output.bam}" -
        """
rule emseq_methyldackel_spike:
    message: "MethylDackel CpG extraction on spike-in alignments, allowing duplicates"
    conda: ENV_EMSEQ
    input:
        bam   = f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel_spike.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel_spike.tsv"
    params:
        out_prefix = lambda wc, input, output: output.bed.rsplit("_CpG.methylKit", 1)[0],
        bounds     = EMSEQ_CALL_BOUNDS,
    threads: 8
    output:
        bed = f"{D_EMSEQ}/spike/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methyldackel-spike] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        mkdir -p "$(dirname "{params.out_prefix}")"
        # --keepDiscordant: count pairs without the proper-pair flag, and
        # exclude the same read-2 bases, as the host calling does
        MethylDackel extract \
          -@ {threads} \
          --methylKit \
          --keepDiscordant \
          {params.bounds} \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """
# ── Reference Indexing ───────────────────────────────────────────────────────
def _emseq_index_fastas(wildcards):
    """FASTA of one reference, or of the sample reference and every control for EMSEQ_SPIKE_REF."""
    name = wildcards.emseq_ref_name
    names = [spike_host_ref, *spike_builds] if name == EMSEQ_SPIKE_REF else [name]
    return [f"{D_INPUTS}/{EMSEQ_REF_INPUTS[n]}" for n in names]

rule emseq_bwa_meth_index:
    message: "Build bwa-meth bisulfite index from reference FASTA"
    conda: ENV_EMSEQ
    input:
        fastas = _emseq_index_fastas,
    log:
        cmd = f"{D_LOGS}/{{emseq_ref_name}}_emseq_bwa_meth_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{emseq_ref_name}}_emseq_bwa_meth_index.tsv"
    params:
        fasta_target = lambda wc: f"{D_REF}/bwa_meth/{wc.emseq_ref_name}/{wc.emseq_ref_name}.fa"
    threads: 1
    output:
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.0123",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.amb",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.ann",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.bwt.2bit.64",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t.pac",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[bwa-meth-index] $(date) ref={wildcards.emseq_ref_name} threads={threads}"
        # awk 1 ends each file with a newline, so the next header starts a line
        for fasta in {input.fastas}; do
            if file -bL "$fasta" | grep -qi gzip; then zcat "$fasta"; else cat "$fasta"; fi | awk 1
        done > "{params.fasta_target}"
        samtools faidx "{params.fasta_target}"
        if [ "$(grep -c '^>' "{params.fasta_target}")" -ne "$(wc -l < "{params.fasta_target}.fai")" ]; then
            echo "[bwa-meth-index] a contig name repeats across the FASTAs" >&2
            exit 1
        fi
        bwameth.py index-mem2 "{params.fasta_target}"
        """
# ── Alignment ────────────────────────────────────────────────────────────────
rule emseq_align_bwameth:
    message: "BWA-meth bisulfite alignment to human reference with coordinate-sorted BAM output"
    conda: ENV_EMSEQ
    input:
        r1  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2  = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        ref = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        c2t = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.bwameth.c2t",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}_emseq_align_bwameth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}_emseq_align_bwameth.tsv"
    params:
        temp_prefix = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.emseq_ref_name}",
    threads: workflow.cores
    resources:
        concurrency = 100,
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.bwa_meth.coorsort.bam",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[bwameth] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} threads={threads}"
        mkdir -p "$(dirname "{params.temp_prefix}")"
        bwameth.py --threads {threads} \
          --reference "{input.ref}" \
          "{input.r1}" "{input.r2}" \
        | samtools view -u - \
        | samtools sort -@ 8 -T "{params.temp_prefix}" -o "{output.bam}"
        # -u: uncompressed BAM avoids compress→decompress overhead in pipe
        """
rule emseq_biscuit_index:
    message: "Build BISCUIT bisulfite index from reference FASTA"
    conda: ENV_EMSEQ
    input:
        lambda wc: f"{D_INPUTS}/{EMSEQ_REF_INPUTS[wc.emseq_ref_name]}"
    log:
        cmd = f"{D_LOGS}/{{emseq_ref_name}}_emseq_biscuit_index.log"
    benchmark:
        f"{D_BENCHMARK}/{{emseq_ref_name}}_emseq_biscuit_index.tsv"
    threads: 1
    output:
        fasta = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        fai = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
        index = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.par.sa",
        biscuit_index_done = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.biscuit.index.done"
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[biscuit-index] $(date) ref={wildcards.emseq_ref_name} threads={threads}"
        mkdir -p "$(dirname "{output.fasta}")"
        zcat "{input}" > "{output.fasta}"
        samtools faidx "{output.fasta}"
        biscuit index "{output.fasta}"
        touch "{output.biscuit_index_done}"
        """
rule emseq_align_biscuit:
    message: "BISCUIT bisulfite alignment to human reference with coordinate-sorted BAM output"
    conda: ENV_EMSEQ
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        fasta = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        fai   = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
        index = f"{D_REF}/biscuit/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.par.sa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}_emseq_align_biscuit.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}_emseq_align_biscuit.tsv"
    params:
        tmp_dir = D_DATA,
    threads: workflow.cores
    resources:
        concurrency = 100
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.biscuit.coorsort.bam",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[biscuit-align] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} threads={threads}"
        mkdir -p "{params.tmp_dir}/tmp"
        biscuit align -@ {threads} "{input.fasta}" "{input.r1}" "{input.r2}" \
        | samtools view -u -t "{input.fai}" \
        | samtools sort -@ {threads} -m 2G \
            -T "{params.tmp_dir}/tmp/{wildcards.library_id}_sorttmp" \
            -o "{output.bam}"
        # -u: uncompressed BAM avoids compress→decompress overhead in pipe
        # -t ref.fai: ensures SQ headers are present for downstream tools
        """
# ── Pre-processing ───────────────────────────────────────────────────────────
rule emseq_fastp:
    message: "Adapter trimming and QC metrics with fastp for EM-seq reads"
    conda: ENV_EMSEQ
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_emseq_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_emseq_fastp.tsv"
    params:
        extra = FASTP_EXTRA,
    threads: 8
    output:
        failed = f"{D_EMSEQ}/fastqs/{{library_id}}.failed.fastq.gz",
        html = f"{D_EMSEQ}/qc/{{library_id}}_emseq_fastp.html",
        json = f"{D_EMSEQ}/qc/{{library_id}}_emseq_fastp.json",
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[fastp] $(date) lib={wildcards.library_id} threads={threads}"
        # Quality filtering disabled — base quality recalibration is unreliable
        # for C-to-T converted libraries
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
    message: "FastQC quality assessment on raw or trimmed FASTQ files"
    conda: ENV_EMSEQ
    input:
        fq = f"{D_EMSEQ}/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{processing}}_{{read}}_emseq_fastqc.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{processing}}_{{read}}_emseq_fastqc.tsv"
    threads: 4
    resources:
        concurrency = 10
    output:
        html = f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.html",
        zip  = f"{D_EMSEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[fastqc] $(date) lib={wildcards.library_id} proc={wildcards.processing} read={wildcards.read} threads={threads}"
        fastqc \
          --outdir "$(dirname "{output.html}")" \
          --quiet \
          --threads {threads} \
          "{input.fq}"
        """
# ── Quality Control ──────────────────────────────────────────────────────────
rule emseq_mosdepth:
    message: "Coverage depth profiling with mosdepth on filtered BAM"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
        index = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam.bai",
        regions_bed = KEEP_BED,
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mosdepth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mosdepth.tsv"
    params:
        script       = f"{R_EMSEQ}/scripts/emseq_mosdepth.sh",
        quant_levels = MOSDEPTH_QUANT_LEVELS,
    threads: 4
    resources:
        concurrency = 20
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
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[mosdepth] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        "{params.script}" \
        "{input.bam}" \
        "$(dirname "{output.summary}")" \
        "{wildcards.library_id}.{wildcards.emseq_ref_name}.{wildcards.align_method}" \
        '{params.quant_levels}' \
        "{input.regions_bed}" \
        {threads}
        """
rule emseq_mbias:
    message: "MethylDackel M-bias detection for position-dependent methylation artifacts"
    conda: ENV_EMSEQ
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        fasta = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.tsv"
    threads: 8
    output:
        txt = f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_mbias.txt",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[mbias] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        # Runs on deduped BAM (not filtered) — full read set needed for bias detection
        # --noSVG: text output only, avoids R dependency
        # --keepDiscordant: count pairs without the proper-pair flag, which
        # bwa-meth leaves unset on cfDNA fragments of about 280 bp or longer
        MethylDackel mbias \
          -@ {threads} \
          --noSVG \
          --keepDiscordant \
          "{input.fasta}" "{input.bam}" > "{output.txt}"
        """
# ── Post-alignment Processing ────────────────────────────────────────────────
# Dedup and filter are deliberately split: dedup MARKS duplicates (0x400)
# so the marked BAM is available for M-bias assessment. Filter then REMOVES
# marked duplicates along with secondary/supplementary/failed-QC reads.
rule emseq_dedup:
    message: "Mark duplicates with dupsifter (marks only, does not remove) on bisulfite-aligned BAM"
    conda: ENV_EMSEQ
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.bam",
        fasta = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_dedup.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_dedup.tsv"
    params:
        temp_prefix = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.emseq_ref_name}.{wc.align_method}.coorsort",
        pair_filter = EMSEQ_PAIR_FILTER,
    threads: 8
    resources:
        concurrency = 25
    output:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        index = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
    shell:
        """
        exec >> "{log.cmd}" 2>&1
        echo "[dedup] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        mkdir -p "$(dirname "{params.temp_prefix}")"
        rm -f {params.temp_prefix}.nsort.tmp.* {params.temp_prefix}.csort.tmp.*
        # Both mates mapped and passing EMSEQ_PAIR_FILTER before name-sorting for dupsifter.
        # The two sorts run at once in this pipe, so each has its own temporary-file prefix.
        samtools view -bh -f 1 -F 12 -e '{params.pair_filter}' "{input.bam}" \
        | samtools sort -n -@ {threads} -O BAM -T "{params.temp_prefix}.nsort.tmp" -o - \
        | dupsifter \
        --add-mate-tags \
        --stats-output "{log.cmd}.dupsifter.tsv" \
        "{input.fasta}" - \
        | samtools sort -@ {threads} -O BAM -T "{params.temp_prefix}.csort.tmp" -o "{output.bam}"
        samtools index -@ {threads} "{output.bam}"
        """
rule emseq_filter_bam:
    message: "Filter BAM to inward-facing pairs <= 1000 bp, MAPQ>=30, autosomes, excluding blacklist and duplicates"
    conda: ENV_EMSEQ
    input:
        bam   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam",
        bai   = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.deduped.bam.bai",
        keep_bed    = KEEP_BED,
        exclude_bed = EXCL_BED,
        fai   = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_filter_bam.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_filter_bam.tsv"
    params:
        pair_filter = EMSEQ_PAIR_FILTER,
    threads: 8
    output:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
        bai = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam.bai",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[filter-bam] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        # -f 1 -F 12 -e pair_filter: both mates mapped, EMSEQ_PAIR_FILTER
        # -q 30: MAPQ >= 30
        # -F 3852 (0xF0C): also remove secondary + supplementary + failed QC + duplicates
        # bedtools -sorted stops reading at the last keep_bed contig; samtools then
        # writes trailing non-primary contigs into a closed pipe (SIGPIPE, exit 141).
        # Scope pipefail off to this one pipe so the benign 141 is not fatal.
        set +o pipefail
        samtools view -@ {threads} -u -f 1 -q 30 -F 3852 -e '{params.pair_filter}' "{input.bam}" \
          | bedtools intersect -sorted -g "{input.fai}" -a stdin -b "{input.exclude_bed}" -v -ubam \
          | bedtools intersect -sorted -g "{input.fai}" -a stdin -b "{input.keep_bed}" -ubam \
          > "{output.bam}.tmp"
        set -o pipefail
        # Prune @SQ dictionary to primary contigs. Non-primary (decoy/unplaced/alt)
        # contigs carry zero reads after the keep-bed filter and trail the primary
        # contigs in the dictionary, so dropping them shifts no read's reference tid.
        samtools reheader \
          <(samtools view -H "{output.bam}.tmp" \
              | grep -E '^@HD|^@PG|^@RG|^@CO|^@SQ[[:space:]]SN:chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]') \
          "{output.bam}.tmp" > "{output.bam}"
        rm -f "{output.bam}.tmp"
        samtools index -@ {threads} "{output.bam}" "{output.bai}"
        """
rule emseq_samtools_stats:
    message: "Samtools stats and flagstat on filtered BAM for alignment QC metrics"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_samtools_stats.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_samtools_stats.tsv"
    threads: 8
    output:
        stats    = f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.stats.txt",
        flagstat = f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.flagstat.txt",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[samtools-stats] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        samtools stats -@ {threads} "{input.bam}" > "{output.stats}"
        samtools flagstat -@ {threads} "{input.bam}" > "{output.flagstat}"
        """
# ── Methylation Calling ──────────────────────────────────────────────────────
rule emseq_methyldackel:
    message: "MethylDackel CpG methylation extraction in methylKit format from filtered BAM"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
        fasta = f"{D_REF}/{{align_method}}/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_methyldackel.tsv"
    params:
        out_prefix = lambda wc, input, output: output.bed.rsplit("_CpG.methylKit", 1)[0],
        bounds     = EMSEQ_CALL_BOUNDS,
    threads: 8
    resources:
        concurrency = 25
    output:
        bed = f"{D_EMSEQ}/meth/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methyldackel] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        # --keepDiscordant: count the pairs EMSEQ_PAIR_FILTER keeps without the
        # proper-pair flag (cfDNA fragments of about 280 bp or longer)
        MethylDackel extract \
          -@ {threads} \
          --methylKit \
          --keepDiscordant \
          {params.bounds} \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """
rule emseq_make_single_methylkit_amp_obj:
    message: "Build tabix-backed methylKit object from MethylDackel AMP output"
    conda: ENV_METHYLKIT
    input:
        amp = f"{D_EMSEQ}/meth/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_methyldackel_CpG.methylKit",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_make_single_methylkit_amp_obj.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_make_single_methylkit_amp_obj.tsv"
    params:
        Rscript   = f"{R_EMSEQ}/scripts/make_single_amp_methylkit_obj.R",
        mincov    = EMSEQ_MINCOV,
        build     = lambda wc: wc.emseq_ref_name,
        treatment = 1,
    threads: 1
    output:
        bgz = f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
        tbi = f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz.tbi",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylKit-amp] $(date) lib={wildcards.library_id} ref={wildcards.emseq_ref_name} aln={wildcards.align_method} threads={threads}"
        Rscript "{params.Rscript}" \
          --amp_file "{input.amp}" \
          --library_id "{wildcards.library_id}.{wildcards.emseq_ref_name}.{wildcards.align_method}.methyldackel" \
          --mincov {params.mincov} \
          --out_dir "$(dirname "{output.bgz}")" \
          --treatment {params.treatment} \
          --build {params.build}
        """
rule emseq_multiqc:
    message: "Aggregate QC metrics into MultiQC report for EM-seq workflow"
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
            align_method=emseq_align_methods,
            allow_missing=True,
        ),
        mosdepth_summary = expand(
            f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.summary.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=emseq_align_methods,
            allow_missing=True,
        ),
        mosdepth_dists = expand(
            f"{D_EMSEQ}/qc/mosdepth_{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mosdepth.{{dist}}.dist.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=emseq_align_methods,
            dist=["global","region"],
            allow_missing=True,
        ),
        samstats = expand(
            f"{D_EMSEQ}/qc/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.samtools.{{stat}}.txt",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=emseq_align_methods,
            stat=["stats","flagstat"],
        ),
    log:
        cmd = f"{D_LOGS}/emseq_multiqc.log",
    benchmark:
        f"{D_BENCHMARK}/emseq_multiqc.tsv"
    params:
        extra = "--force",
    threads: 1
    output:
        html = f"{D_EMSEQ}/qc/multiqc.html",
        data = directory(f"{D_EMSEQ}/qc/multiqc_data"),
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[multiqc] $(date) threads={threads}"
        mkdir -p "$(dirname "{output.html}")"
        multiqc \
            {input} \
            {params.extra} \
            --outdir "$(dirname "{output.html}")" \
            --filename "$(basename "{output.html}")" \
        """
