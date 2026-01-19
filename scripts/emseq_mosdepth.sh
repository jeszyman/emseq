#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# mosdepth-wrapper.sh
#
# This script wraps the `mosdepth` tool to compute read depth over a BAM file,
# optimized for EM-seq cfDNA data. It configures the run to:
#   - use median depth (`--use-median`)
#   - run in fast mode (no per-base depth)
#   - report thresholds and quantized bins
#   - restrict to user-provided target regions
#
# Output files are written using a prefix of "mosdepth_<OUT_PREFIX>" in <OUT_DIR>.
# Designed for use in explicit I/O workflows like Snakemake or manual batch calls.
# -----------------------------------------------------------------------------
print_usage() {
    cat <<EOF
USAGE: mosdepth-wrapper.sh <BAM> <OUT_DIR> <OUT_PREFIX> <QUANT_LEVELS> <REGIONS_BED> [THREADS]
DESCRIPTION:
  Run mosdepth on a BAM file with EM-seq-appropriate settings.
  QUANT_LEVELS is a comma-separated string of coverage cutoffs (e.g. 1,5,10,20).
  REGIONS_BED is a BED file defining target regions (e.g. autosomes minus blacklist).
  The OUT_PREFIX will be prepended with 'mosdepth_' before being passed to mosdepth.
  Output files (e.g. mosdepth_<OUT_PREFIX>.summary.txt) will be written to OUT_DIR.
  THREADS is optional (default: 8).
EOF
}
main() {
    parse_args "$@"
    run_mosdepth
}
parse_args() {
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        print_usage
        exit 0
    fi
    if [[ $# -lt 5 ]]; then
        echo "Error: Missing required arguments." >&2
        print_usage
        exit 1
    fi
    declare -g bam_file="$1"
    declare -g out_dir="$2"
    declare -g user_prefix="$3"
    declare -g quant_levels="$4"
    declare -g regions_bed="$5"
    declare -g threads="${6:-8}"
    [[ -f "$bam_file" ]] || { echo "Error: BAM file not found: $bam_file" >&2; exit 1; }
    [[ -f "$regions_bed" ]] || { echo "Error: Regions BED not found: $regions_bed" >&2; exit 1; }
    mkdir -p "$out_dir"
    declare -g out_prefix="mosdepth_${user_prefix}"
    declare -g out_path="${out_dir%/}/${out_prefix}"
    declare -g quant_str="0:${quant_levels/,/:}"
}
run_mosdepth() {
    echo "[INFO] PID $$ running mosdepth on $bam_file" >&2
    echo "[INFO] Output prefix: $out_path" >&2
    echo "[INFO] Regions BED: $regions_bed" >&2
    echo "[INFO] Quantize string: $quant_str" >&2
    echo "[INFO] Threads: $threads" >&2
    mosdepth \
        --threads "$threads" \
        --no-per-base \
        --fast-mode \
        --use-median \
        --quantize "$quant_str" \
        --by "$regions_bed" \
        --thresholds "$quant_levels" \
        "$out_path" "$bam_file"
    echo "[INFO] mosdepth complete for PID $$" >&2
}
main "$@"
