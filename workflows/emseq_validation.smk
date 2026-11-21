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
#   VALIDATION_HOST_REF            name of the sample reference
#   VALIDATION_FRAGMENT_SOURCES    sources tallied by distance from the strand ends (puc19, unmeth_lambda, host)
#   VALIDATION_CHECK_INDEX         index of the library that checks the measures, or None
#   VALIDATION_END_SETTINGS        dict from calling setting to its MethylDackel flags, empty to skip the comparison
#   VALIDATION_FRAGMENT_PANELS     list of (data set, source, panel title, excluded library index) for the strand-end figure

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
        # --keepDiscordant: count pairs without the proper-pair flag, as the
        # pipeline's calling rules do.
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

rule emseq_validation_control_reads:
    message: "Usable read count and fragment lengths of a control alignment"
    conda: ENV_EMSEQ
    input:
        bam = f"{D_EMSEQ}/spike/{{library_id}}.{{source}}.bwa_meth.coorsort.bam",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{source}}_emseq_validation_control_reads.log",
    wildcard_constraints:
        source = "puc19|unmeth_lambda",
    threads: 2
    output:
        reads = f"{D_EMSEQ}/validation/control/{{library_id}}.{{source}}.reads.txt",
        tlen  = f"{D_EMSEQ}/validation/control/{{library_id}}.{{source}}.tlen.tsv",
    shell:
        """
        exec 2>> "{log.cmd}"
        echo "[validation-control-reads] $(date) lib={wildcards.library_id} source={wildcards.source}" >&2
        # Usable reads are the ones MethylDackel counts by default: MAPQ 10 or
        # more, primary, not QC-fail, not duplicate (-F 0xF00). bwa-meth marks
        # short partial alignments of host reads to a control as QC-fail.
        samtools view -@ {threads} -c -q 10 -F 0xF00 "{input.bam}" > "{output.reads}"
        # Fragment length (TLEN) of usable pairs, taken from the mate with a
        # positive TLEN, as length and count
        samtools view -@ {threads} -q 10 -F 0xF00 -f 1 "{input.bam}" \
          | awk '$9 > 0 && $9 <= 1000 {{n[$9]++}} END {{for (l in n) print l "\t" n[l]}}' \
          | sort -k1,1n > "{output.tlen}"
        """

def _validation_fragment_bam(wildcards):
    """Filtered host alignment for the sample reads, spike alignment for a control."""
    if wildcards.source == "host":
        return f"{D_EMSEQ}/bams/{wildcards.library_id}.{VALIDATION_HOST_REF}.bwa_meth.coorsort.filt.bam"
    return f"{D_EMSEQ}/spike/{wildcards.library_id}.{wildcards.source}.bwa_meth.coorsort.bam"

def _validation_fragment_fasta(wildcards):
    name = VALIDATION_HOST_REF if wildcards.source == "host" else wildcards.source
    return f"{D_REF}/bwa_meth/{name}/{name}.fa"

rule emseq_validation_fragment_end:
    message: "CpG calls by distance from the 5' and 3' ends of the original strand"
    conda: ENV_EMSEQ
    input:
        bam   = _validation_fragment_bam,
        fasta = _validation_fragment_fasta,
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{source}}_emseq_validation_fragment_end.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{source}}_emseq_validation_fragment_end.tsv"
    params:
        script = f"{R_EMSEQ}/scripts/emseq_validation_fragment_end.py",
    wildcard_constraints:
        source = "puc19|unmeth_lambda|host",
    threads: 1
    output:
        ends      = f"{D_EMSEQ}/validation/fragment_end/{{library_id}}.{{source}}.ends.tsv",
        agreement = f"{D_EMSEQ}/validation/fragment_end/{{library_id}}.{{source}}.agreement.tsv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[validation-fragment-end] $(date) lib={wildcards.library_id} source={wildcards.source}"
        python "{params.script}" \
          --bam "{input.bam}" \
          --fasta "{input.fasta}" \
          --out-ends "{output.ends}" \
          --out-agreement "{output.agreement}"
        """

VALIDATION_FRAGMENT_TABLES = [
    f"{D_EMSEQ}/validation/fragment_end/{library_id}.{source}.ends.tsv"
    for library_id in emseq_library_ids
    for source in VALIDATION_FRAGMENT_SOURCES
]
VALIDATION_EXTRACT_TABLES = [
    f"{D_EMSEQ}/validation/extract/{library_id}.{end_setting}_CpG.bedGraph"
    for library_id in emseq_library_ids
    for end_setting in VALIDATION_END_SETTINGS
]

def _validation_values_extra_args():
    """Arguments of the values script for the plasma measures that the config switches on."""
    args = []
    if VALIDATION_FRAGMENT_SOURCES:
        args += ["--fragment-dir", f"{D_EMSEQ}/validation/fragment_end",
                 "--fragment-sources", *VALIDATION_FRAGMENT_SOURCES]
    if VALIDATION_CHECK_INDEX is not None:
        args += ["--check-index", str(VALIDATION_CHECK_INDEX)]
    if VALIDATION_END_SETTINGS:
        args += ["--extract-dir", f"{D_EMSEQ}/validation/extract",
                 "--extract-settings", *VALIDATION_END_SETTINGS]
    return " ".join(args)

rule emseq_validation_extract:
    message: "CpG calls on the sample alignment under one read-end setting"
    conda: ENV_EMSEQ
    input:
        bam   = lambda wc: f"{D_EMSEQ}/bams/{wc.library_id}.{VALIDATION_HOST_REF}.bwa_meth.coorsort.filt.bam",
        fasta = f"{D_REF}/bwa_meth/{VALIDATION_HOST_REF}/{VALIDATION_HOST_REF}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{end_setting}}_emseq_validation_extract.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{end_setting}}_emseq_validation_extract.tsv"
    params:
        flags = lambda wc: VALIDATION_END_SETTINGS[wc.end_setting],
        out_prefix = lambda wc, output: output.bedgraph.rsplit("_CpG.bedGraph", 1)[0],
    threads: 4
    output:
        bedgraph = f"{D_EMSEQ}/validation/extract/{{library_id}}.{{end_setting}}_CpG.bedGraph",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[validation-extract] $(date) lib={wildcards.library_id} setting={wildcards.end_setting}"
        MethylDackel extract \
          -@ {threads} \
          --keepDiscordant \
          {params.flags} \
          "{input.fasta}" \
          "{input.bam}" \
          -o "{params.out_prefix}"
        """

rule emseq_validation_values:
    message: "Numbers and position table for the read-end validation report"
    conda: ENV_EMSEQ
    input:
        tables = VALIDATION_MBIAS_TABLES,
        fragment = VALIDATION_FRAGMENT_TABLES,
        extract = VALIDATION_EXTRACT_TABLES,
        control = expand(
            f"{D_EMSEQ}/validation/control/{{library_id}}.{{source}}.{{kind}}",
            library_id=emseq_library_ids,
            source=list(VALIDATION_CONTROL_CONTEXTS),
            kind=["reads.txt", "tlen.tsv"],
        ),
    log:
        cmd = f"{D_LOGS}/emseq_validation_values.log",
    params:
        script = f"{R_EMSEQ}/scripts/emseq_validation_end_handling.py",
        mbias_dir = f"{D_EMSEQ}/validation/mbias",
        control_dir = f"{D_EMSEQ}/validation/control",
        out_dir = f"{R_EMSEQ}/results/end-handling-{VALIDATION_DATASET}",
        index = " ".join(f"{lib}={idx}" for lib, idx in VALIDATION_INDEX.items()),
        extra = _validation_values_extra_args(),
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
          --control-dir "{params.control_dir}" \
          --out-dir "{params.out_dir}" \
          --library-index {params.index} \
          --methyldackel-version "$(MethylDackel --version 2>&1 | head -n 1)" \
          --pipeline-commit "$(git -C "{params.repo}" rev-parse --short HEAD)" \
          {params.extra}
        """

rule emseq_validation_fragment_figure:
    message: "Figure of methylation by distance from the strand ends"
    conda: ENV_EMSEQ
    input:
        tables = [f"{R_EMSEQ}/results/end-handling-{dataset}/fragment_end.tsv"
                  for dataset, _, _, _ in VALIDATION_FRAGMENT_PANELS],
    log:
        cmd = f"{D_LOGS}/emseq_validation_fragment_figure.log",
    params:
        script = f"{R_EMSEQ}/scripts/emseq_validation_fragment_end_plot.R",
        sources = " ".join(source for _, source, _, _ in VALIDATION_FRAGMENT_PANELS),
        labels = " ".join(f"'{label}'" for _, _, label, _ in VALIDATION_FRAGMENT_PANELS),
        excluded = " ".join(str(index) for _, _, _, index in VALIDATION_FRAGMENT_PANELS),
    output:
        png = f"{R_EMSEQ}/resources/figures/validation_fragment_end_{VALIDATION_DATASET}.png",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[validation-fragment-figure] $(date) dataset={VALIDATION_DATASET}"
        Rscript "{params.script}" \
          --table {input.tables} \
          --source {params.sources} \
          --label {params.labels} \
          --exclude-index {params.excluded} \
          --out "{output.png}"
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
