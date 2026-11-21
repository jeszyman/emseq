#!/usr/bin/env bash
# Run the four steps of JEEPERS (Wang 2026, sibling mode), unmodified, on one
# prepared BAM, and check that it returned every fragment it was given.
#
# Usage: emseq_jeepers_run.sh <in.bam> <regions.bed> <outdir> <0.9|default> [jeepers_dir]
#
#   in.bam       BAM from emseq_jeepers_input.py
#   regions.bed  CpG region file from emseq_jeepers_regions.py
#   0.9|default  0.9 passes -m 0.9 to the correction step (the value the
#                JEEPERS README states); default passes no -m, which leaves
#                the code's own default of 0.1
#   jeepers_dir  folder holding the JEEPERS scripts
#                (default /mnt/data/tools/jeepers/JEEPERS_code)
#
# Environment:
#   JEEPERS_CORES       processes for each JEEPERS step (default 16)
#   JEEPERS_PYTHON      interpreter for the JEEPERS scripts (default python)
#   SAMTOOLS            samtools executable (default samtools)
#   JEEPERS_STOP_AFTER  "merge" stops after the merge step, for a run that
#                       needs only the table of fragments and their calls
#
# Leaves in <outdir>: t_Merge.txt (fragments given to the correction step),
# t_JEEPERS_Sibling.txt (fragments after correction), one log per step,
# run.tsv (step, seconds, rows), and blocks.tsv (rows per block in each
# table).
#
# Exit 1 when a table is empty, when a log holds "Traceback" or "ERROR", or
# when a block has a different number of rows after correction than before.
# The JEEPERS scripts exit 0 with empty output when they cannot read their
# input, so their exit codes are not the check.
set -uo pipefail

if [[ $# -lt 4 ]]; then
    sed -n '2,6p' "$0" >&2
    exit 2
fi
in_bam=$(realpath "$1")
regions=$(realpath "$2")
outdir=$3
setting=$4
jeepers=${5:-/mnt/data/tools/jeepers/JEEPERS_code}
cores=${JEEPERS_CORES:-16}
py=${JEEPERS_PYTHON:-python}
samtools=${SAMTOOLS:-samtools}

case "$setting" in
    0.9) meth_args=(-m 0.9) ;;
    default) meth_args=() ;;
    *) echo "setting must be 0.9 or default, got: $setting" >&2; exit 2 ;;
esac

mkdir -p "$outdir"
outdir=$(realpath "$outdir")
jeepers=$(realpath "$jeepers")
cd "$outdir" || exit 2
printf 'step\tseconds\trows\n' > run.tsv

# step <name> <file whose rows are counted> <command...>
step() {
    local name=$1 counted=$2
    shift 2
    local t0
    t0=$(date +%s)
    "$@" > "$name.log" 2>&1
    local rows=0
    [[ -f $counted ]] && rows=$(wc -l < "$counted")
    printf '%s\t%s\t%s\n' "$name" "$(( $(date +%s) - t0 ))" "$rows" >> run.tsv
    if [[ $rows -eq 0 ]]; then
        echo "FAIL: step $name left an empty table ($counted)" >&2
        exit 1
    fi
}

# rows per block (column 4), skipping a header line if one is present
block_rows() { awk -F'\t' '$4 != "Block" {n[$4]++} END {for (b in n) print b "\t" n[b]}' "$1" | sort -k1,1g; }

step format tagged.bam "$py" "$jeepers/FormatChange/bamFormatChange.py" --ibam "$in_bam" --obam "$outdir/tagged.bam" --cpu_num 1
"$samtools" index tagged.bam
step bam2txt t.B2T.txt "$py" "$jeepers/Bam2Txt.py" --ibam "$outdir/tagged.bam" --ibed "$regions" --ofile "$outdir/t.B2T.txt" -c "$cores" -d T
step merge t_Merge.txt "$py" "$jeepers/MergeR1R2.py" -i "$outdir/t.B2T.txt" -o "$outdir/t_Merge.txt" -c "$cores"

status=0
logs=(format.log bam2txt.log merge.log)
if [[ ${JEEPERS_STOP_AFTER:-} != merge ]]; then
    step correct t_JEEPERS_Sibling.txt "$py" "$jeepers/JEEPERS_correction.py" -i "$outdir/t_Merge.txt" -o "$outdir/t_JEEPERS.txt" -c "$cores" -n 4 "${meth_args[@]}"
    logs+=(correct.log)
    join -t $'\t' -a 1 -a 2 -e 0 -o 0,1.2,2.2 <(block_rows t_Merge.txt | sort -k1,1) <(block_rows t_JEEPERS_Sibling.txt | sort -k1,1) \
        | sort -k1,1g | { printf 'block\tmerge_rows\tsibling_rows\n'; cat; } > blocks.tsv
    if awk -F'\t' 'NR > 1 && $2 != $3 {bad = 1} END {exit !bad}' blocks.tsv; then
        echo "FAIL: rows after correction differ from rows before in at least one block (blocks.tsv)" >&2
        status=1
    fi
fi
if grep -l -e Traceback -e ERROR "${logs[@]}" >&2; then
    echo "FAIL: Traceback or ERROR in the logs named above" >&2
    status=1
fi
exit $status
