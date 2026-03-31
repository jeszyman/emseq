# ============================================================================
# emseq_analysis.smk — EM-seq downstream analysis module
# Tangled from emseq.org; do not edit directly.
# ============================================================================
# Required wrapper-provided variables (in addition to emseq.smk variables):
#   R_MHAPTOOLS, R_WGBSTOOLS, R_UXM — external tool repository paths
#   CPG_REF          — tabixed CpG reference file
#   MHB_BED          — methylation haplotype block BED
#   HAP_METRICS      — space-separated haplotype metric names
#   DECONV_GENOME    — wgbs_tools genome identifier
#   DECONV_ATLAS     — UXM deconvolution atlas path
#   ENV_HAPLOTYPE    — haplotype conda env YAML path
#   ENV_DECONV       — deconvolution conda env YAML path
rule emseq_analysis_methylkit_unite:
    message: "Unite per-sample methylKit tabix databases into single methylBase for differential analysis"
    wildcard_constraints:
        experiment = "[^.]+",
    conda: ENV_METHYLKIT
    input:
        mkit_lib_db = lambda wc: expand(
            f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
            library_id=meth_map[wc.experiment]["libs"],
            emseq_ref_name=meth_map[wc.experiment]["emseq_ref_name"],
            align_method=meth_map[wc.experiment]["align_method"],
        ),
    log:
        cmd = f"{D_LOGS}/{{experiment}}_emseq_analysis_methylkit_unite.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_emseq_analysis_methylkit_unite.tsv"
    params:
        library_id = lambda wc: " ".join(meth_map[wc.experiment]["libs"]),
        treatment_list = lambda wc: " ".join(map(str, meth_map[wc.experiment]["tx"])),
        script = f"{R_EMSEQ}/scripts/make_methylkit_unite_db.R",
        mincov = lambda wc: meth_map[wc.experiment]["mincov"],
        min_per_group = lambda wc: meth_map[wc.experiment]["mingroup"],
        chunk_size = lambda wc: meth_map[wc.experiment]["chunksize"],
    threads: 32
    output:
        mbase = f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-unite] $(date) experiment={wildcards.experiment} threads={threads}"
        rm -f "{output.mbase}"*
        Rscript "{params.script}" \
          --lib_db_list "{input.mkit_lib_db}" \
          --lib_id_list "{params.library_id}" \
          --treatment_list "{params.treatment_list}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.mbase}")" \
          --suffix {wildcards.experiment} \
          --assembly "hg38" \
          --mincov {params.mincov} \
          --min_per_group {params.min_per_group} \
          --chunk_size {params.chunk_size}
        """
rule emseq_analysis_methylkit_diff:
    message: "Calculate differential methylation from united methylBase using methylKit"
    wildcard_constraints:
        experiment = "[^.]+",
    conda: ENV_METHYLKIT
    input:
        mbase = f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
    log:
        cmd = f"{D_LOGS}/{{experiment}}_emseq_analysis_methylkit_diff.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_emseq_analysis_methylkit_diff.tsv"
    params:
        script     = f"{R_EMSEQ}/scripts/make_methylkit_diff_db.R",
        chunk_size = lambda wc: meth_map[wc.experiment]["chunksize"],
    threads: 32
    output:
        mdiff = f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-diff] $(date) experiment={wildcards.experiment} threads={threads}"
        rm -f "{output.mdiff}" "{output.mdiff}.tbi"
        Rscript "{params.script}" \
          --mbase "{input.mbase}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.mdiff}")" \
          --suffix {wildcards.experiment} \
          --chunk_size {params.chunk_size}
        """
rule emseq_analysis_methylkit_diff_tiled:
    message: "Tile methylation data into windows and calculate differential methylation with methylKit"
    conda: ENV_METHYLKIT
    input:
        mkit_lib_db = lambda wc: expand(
            f"{D_EMSEQ}/dmr/tabix/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.methyldackel.txt.bgz",
            library_id=meth_map[wc.experiment]['libs'],
            emseq_ref_name=meth_map[wc.experiment]['emseq_ref_name'],
            align_method=meth_map[wc.experiment]['align_method'],
        ),
    log:
        cmd = f"{D_LOGS}/{{experiment}}_emseq_analysis_methylkit_diff_tiled.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_emseq_analysis_methylkit_diff_tiled.tsv"
    params:
        library_id = lambda wc: " ".join(meth_map[wc.experiment]['libs']),
        treatment_list = lambda wc: " ".join(map(str, meth_map[wc.experiment]['tx'])),
        script = f"{R_EMSEQ}/scripts/make_methylkit_diff_tiled_db.R",
        mincov = lambda wc: meth_map[wc.experiment]["mincov"],
        chunk_size = lambda wc: meth_map[wc.experiment]["chunksize"],
        min_per_group = lambda wc: meth_map[wc.experiment]["mingroup"],
        win_size = lambda wc: meth_map[wc.experiment]["win_size"],
    threads: 32
    output:
        unite = f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz",
        diff = f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-diff-tiled] $(date) experiment={wildcards.experiment} threads={threads}"
        rm -f "{output.unite}"* "{output.diff}"*
        Rscript "{params.script}" \
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
          --chunk_size {params.chunk_size}
        """
rule emseq_analysis_methylkit_meth_extract:
    message: "Extract percent methylation matrix from united methylBase"
    conda: ENV_METHYLKIT
    input:
        mbase = f"{D_EMSEQ}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
    log:
        cmd = f"{D_LOGS}/{{experiment}}_emseq_analysis_methylkit_meth_extract.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_emseq_analysis_methylkit_meth_extract.tsv"
    params:
        script = f"{R_EMSEQ}/scripts/all_experiment_methylation.R",
    threads: 1
    output:
        tsv = f"{D_EMSEQ}/dmr/diff/{{experiment}}_pos_meth.tsv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-meth-extract] $(date) experiment={wildcards.experiment}"
        Rscript "{params.script}" \
          --db_file "{input.mbase}" \
          --out_file "{output.tsv}"
        """
rule emseq_analysis_annotate_cpg:
    message: "Annotate methylKit diff/base DB with CpG context, gene parts, and nearest TSS"
    conda: ENV_METHYLKIT
    input:
        db = f"{D_EMSEQ}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
    log:
        cmd = f"{D_LOGS}/{{experiment}}_emseq_analysis_annotate_cpg.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_emseq_analysis_annotate_cpg.tsv"
    params:
        script = f"{R_EMSEQ}/scripts/emseq_annotate_methylkit.R",
    threads: 1
    output:
        tsv = f"{D_EMSEQ}/dmr/annotation/{{experiment}}_annotated.tsv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[annotate-cpg] $(date) experiment={wildcards.experiment}"
        mkdir -p "$(dirname "{output.tsv}")"
        Rscript "{params.script}" \
          --db "{input.db}" \
          --out "{output.tsv}"
        """
rule emseq_analysis_mhap_convert:
    message: "Convert filtered BAM to mhap format using mHapTools for haplotype analysis"
    conda: ENV_HAPLOTYPE
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_mhap_convert.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_mhap_convert.tsv"
    params:
        mhaptools = f"{R_MHAPTOOLS}/mhaptools",
        cpg_ref   = CPG_REF,
    threads: 1
    output:
        mhap = f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mhap.gz",
        tbi  = f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mhap.gz.tbi",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[mhap-convert] $(date) lib={wildcards.library_id}"
        mkdir -p "$(dirname "{output.mhap}")"
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        "{params.mhaptools}" convert \
          -i "{input.bam}" \
          -c "{params.cpg_ref}" \
          -o "{output.mhap}"
        tabix -b 2 -e 3 -p bed "{output.mhap}"
        """
rule emseq_analysis_mhaptk_stat:
    message: "Compute MHL, PDR, and Entropy per MHB region using mhaptk stat"
    conda: ENV_HAPLOTYPE
    input:
        mhap = f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mhap.gz",
        tbi  = f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.mhap.gz.tbi",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_mhaptk_stat.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_mhaptk_stat.tsv"
    params:
        cpg_ref = CPG_REF,
        mhb_bed = MHB_BED,
        metrics = HAP_METRICS,
    threads: 1
    output:
        txt = f"{D_EMSEQ}/haplotype/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_mhl.txt",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[mhaptk-stat] $(date) lib={wildcards.library_id}"
        mhaptk stat \
          --mhapPath "{input.mhap}" \
          --cpgPath "{params.cpg_ref}" \
          --bedPath "{params.mhb_bed}" \
          --metrics {params.metrics} \
          --outputFile "{output.txt}"
        """
rule emseq_analysis_chr_reheader:
    message: "Add chr prefix to BAM contig names for wgbs_tools/UXM compatibility"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/bams/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.coorsort.filt.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_chr_reheader.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_chr_reheader.tsv"
    threads: 4
    output:
        bam = f"{D_EMSEQ}/deconv/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.chr.bam",
        bai = f"{D_EMSEQ}/deconv/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.chr.bam.bai",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[chr-reheader] $(date) lib={wildcards.library_id}"
        mkdir -p "$(dirname "{output.bam}")"
        samtools view -H "{input.bam}" \
          | sed 's/SN:\\([0-9XYM]\\)/SN:chr\\1/g' \
          | samtools reheader - "{input.bam}" > "{output.bam}"
        samtools index -@ {threads} "{output.bam}"
        """
rule emseq_analysis_wgbstools_init_genome:
    message: "Initialize wgbs_tools genome reference from indexed FASTA"
    conda: ENV_DECONV
    input:
        fasta = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa",
        fai   = f"{D_REF}/bwa_meth/{{emseq_ref_name}}/{{emseq_ref_name}}.fa.fai",
    log:
        cmd = f"{D_LOGS}/{{emseq_ref_name}}_emseq_analysis_wgbstools_init_genome.log",
    benchmark:
        f"{D_BENCHMARK}/{{emseq_ref_name}}_emseq_analysis_wgbstools_init_genome.tsv"
    params:
        wgbstools = f"{R_WGBSTOOLS}/wgbstools",
        genome    = DECONV_GENOME,
    threads: 1
    output:
        done = f"{D_REF}/wgbstools/{{emseq_ref_name}}.init.done",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[wgbstools-init] $(date) ref={wildcards.emseq_ref_name} genome={params.genome}"
        mkdir -p "$(dirname "{output.done}")"
        "{params.wgbstools}" init_genome \
          --fasta_path "{input.fasta}" \
          -f \
          "{params.genome}"
        touch "{output.done}"
        """
rule emseq_analysis_bam2pat:
    message: "Convert chr-prefixed BAM to methylation pat file using wgbs_tools"
    conda: ENV_DECONV
    input:
        bam = f"{D_EMSEQ}/deconv/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.chr.bam",
        bai = f"{D_EMSEQ}/deconv/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.chr.bam.bai",
        genome_init = f"{D_REF}/wgbstools/{{emseq_ref_name}}.init.done",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_bam2pat.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{emseq_ref_name}}.{{align_method}}_emseq_analysis_bam2pat.tsv"
    params:
        wgbstools = f"{R_WGBSTOOLS}/wgbstools",
        genome    = DECONV_GENOME,
        out_dir   = f"{D_EMSEQ}/deconv",
    threads: 1
    output:
        pat = f"{D_EMSEQ}/deconv/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.chr.pat.gz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[bam2pat] $(date) lib={wildcards.library_id}"
        "{params.wgbstools}" bam2pat \
          --genome "{params.genome}" \
          --out_dir "{params.out_dir}" \
          "{input.bam}"
        """
rule emseq_analysis_uxm_deconv:
    message: "Run UXM tissue deconvolution on all sample pat files"
    conda: ENV_DECONV
    input:
        pats = expand(
            f"{D_EMSEQ}/deconv/{{library_id}}.{{emseq_ref_name}}.{{align_method}}.chr.pat.gz",
            library_id=emseq_library_ids,
            emseq_ref_name=emseq_ref_names,
            align_method=["bwa_meth"],
        ),
    log:
        cmd = f"{D_LOGS}/emseq_analysis_uxm_deconv.log",
    benchmark:
        f"{D_BENCHMARK}/emseq_analysis_uxm_deconv.tsv"
    params:
        uxm   = f"{R_UXM}/uxm",
        atlas = DECONV_ATLAS,
    threads: 1
    output:
        csv = f"{D_EMSEQ}/deconv/uxm_results.csv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[uxm-deconv] $(date) n_samples=$(echo {input.pats} | wc -w)"
        export PATH="{R_WGBSTOOLS}:$PATH"
        "{params.uxm}" deconv \
          --atlas "{params.atlas}" \
          --output "{output.csv}" \
          {input.pats}
        """
