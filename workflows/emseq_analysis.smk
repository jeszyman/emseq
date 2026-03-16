# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/emseq/emseq.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-16 11:33:40
# ============================================================

# emseq_analysis.smk — EM-seq downstream analysis module
# Tangled from emseq.org; do not edit directly.
# ============================================================================
# Required wrapper-provided variables (must be in scope before include):
#   D_EMSEQ, D_LOGS, D_BENCHMARK, D_DATA, R_EMSEQ
#   ENV_METHYLKIT, ENV_EMSEQ
#   meth_map, emseq_library_ids, emseq_ref_names, KEEP_BED, EXCL_BED

# --- Analysis-specific config ---
R_MHAPTOOLS   = config["repos"]["mhaptools"]
R_WGBSTOOLS   = config["repos"]["wgbs_tools"]
R_UXM         = config["repos"]["uxm_deconv"]
CPG_REF       = config["haplotype"]["cpg-ref"]
MHB_BED       = config["haplotype"]["mhb-bed"]
HAP_METRICS   = " ".join(config["haplotype"]["metrics"])
DECONV_GENOME = config["deconv"]["genome-name"]
DECONV_ATLAS  = config["deconv"]["atlas"]

ENV_HAPLOTYPE = config["envs"]["haplotype"]
ENV_DECONV    = config["envs"]["deconv"]
rule emseq_analysis_methylkit_unite:
    message: "Unite per-sample methylKit tabix databases into single methylBase for differential analysis"
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
