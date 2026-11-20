#################################
###   EM-seq validation rules  ###
#################################
#
# Rules for validations that run on data held outside the repository.
#
# Required wrapper-provided variables (must be in scope before include):
#   ENV_EMSEQ                      conda environment YAML path
#   R_EMSEQ                        repository root path
#   D_EMSEQ, D_REF, D_LOGS, D_BENCHMARK   directory layout
#   emseq_library_ids              list of library IDs
#   VALIDATION_DATASET             label of the data set in result paths
#   VALIDATION_INDEX               dict from library ID to its index in result files
#   VALIDATION_CONTROL_CONTEXTS    dict from control reference name to its cytosine contexts

# MethylDackel mbias options for one cytosine context
VALIDATION_CONTEXT_FLAGS = {"CpG": "", "CHG": "--noCpG --CHG", "CHH": "--noCpG --CHH"}

VALIDATION_MBIAS_TABLES = [
    f"{D_EMSEQ}/validation/mbias/{library_id}.{source}.{context}.txt"
    for library_id in emseq_library_ids
    for source, contexts in VALIDATION_CONTROL_CONTEXTS.items()
    for context in contexts
]

rule emseq_validation_mbias:
    message: "Methylation by read position, with MethylDackel's suggested bounds, on a control alignment"
    conda: ENV_EMSEQ
    input:
        bam   = f"{D_EMSEQ}/spike/{{library_id}}.{{source}}.bwa_meth.coorsort.bam",
        fasta = f"{D_REF}/bwa_meth/{{source}}/{{source}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{source}}.{{context}}_emseq_validation_mbias.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{source}}.{{context}}_emseq_validation_mbias.tsv"
    params:
        flags = lambda wc: VALIDATION_CONTEXT_FLAGS[wc.context],
        svg_prefix = lambda wc, output: output.txt.rsplit(".txt", 1)[0],
    wildcard_constraints:
        source = "puc19|unmeth_lambda",
        context = "CpG|CHG|CHH",
    threads: 2
    output:
        txt    = f"{D_EMSEQ}/validation/mbias/{{library_id}}.{{source}}.{{context}}.txt",
        bounds = f"{D_EMSEQ}/validation/mbias/{{library_id}}.{{source}}.{{context}}.bounds.txt",
    shell:
        """
        exec 2>> "{log.cmd}"
        echo "[validation-mbias] $(date) lib={wildcards.library_id} source={wildcards.source} context={wildcards.context}" >&2
        # An output prefix (no --noSVG) makes MethylDackel print its suggested
        # bounds, on stderr; --txt sends the table to stdout.
        # --keepDiscordant: bwa-meth leaves the proper-pair flag unset on
        # nearly all control pairs.
        # -D: the control alignments are deeper than the default cap of 2000.
        MethylDackel mbias \
          -@ {threads} \
          --txt \
          --keepDiscordant \
          -D 1000000 \
          {params.flags} \
          "{input.fasta}" "{input.bam}" "{params.svg_prefix}" \
          > "{output.txt}" 2> "{output.bounds}.stderr"
        cat "{output.bounds}.stderr" >&2
        # grep exits 1 when MethylDackel printed no suggestion; the values
        # script reports that with the file named
        grep "Suggested inclusion options" "{output.bounds}.stderr" > "{output.bounds}" || true
        rm -f "{output.bounds}.stderr" "{params.svg_prefix}"_O?.svg "{params.svg_prefix}"_CTO?.svg
        """

rule emseq_validation_values:
    message: "Numbers and position table for the read-end validation report"
    conda: ENV_EMSEQ
    input:
        tables = VALIDATION_MBIAS_TABLES,
    log:
        cmd = f"{D_LOGS}/emseq_validation_values.log",
    params:
        script = f"{R_EMSEQ}/scripts/emseq_validation_end_handling.py",
        mbias_dir = f"{D_EMSEQ}/validation/mbias",
        out_dir = f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}",
        index = " ".join(f"{lib}={idx}" for lib, idx in VALIDATION_INDEX.items()),
        repo = R_EMSEQ,
    output:
        values   = f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}/values.json",
        position = f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}/position.tsv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[validation-values] $(date) dataset={VALIDATION_DATASET}"
        python "{params.script}" \
          --mbias-dir "{params.mbias_dir}" \
          --out-dir "{params.out_dir}" \
          --library-index {params.index} \
          --methyldackel-version "$(MethylDackel --version 2>&1 | head -n 1)" \
          --pipeline-commit "$(git -C "{params.repo}" rev-parse --short HEAD)"
        """

rule emseq_validation_figure:
    message: "Figure of control methylation by read position"
    conda: ENV_EMSEQ
    input:
        position = f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}/position.tsv",
        values   = f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}/values.json",
    log:
        cmd = f"{D_LOGS}/emseq_validation_figure.log",
    params:
        script = f"{R_EMSEQ}/scripts/emseq_validation_end_handling_plot.R",
    output:
        png = f"{R_EMSEQ}/resources/figures/validation_end_handling_{VALIDATION_DATASET}.png",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[validation-figure] $(date) dataset={VALIDATION_DATASET}"
        Rscript "{params.script}" \
          --position "{input.position}" \
          --values "{input.values}" \
          --out "{output.png}"
        """
