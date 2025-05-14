#!/usr/bin/env bash
set -euo pipefail

print_usage() {
    cat <<EOF
USAGE: fastp-emseq-wrapper.sh <INPUT R1 FASTQ.GZ> <OUTPUT R1> <OUTPUT R2> <FAILED OUT> <LOG TXT> <LOG JSON> <LOG HTML> [THREADS]

DESCRIPTION:
  fastp wrapper for EM-seq data using Snakemake-style explicit I/O.
  Default threads is 16 (max for fastp).
EOF
}

main() {
    parse_args "$@"

    echo "Running fastp on: $in_r1 and $in_r2" | tee "$log_txt"
    echo "Output files: $out_r1, $out_r2, $failed_out" | tee -a "$log_txt"
    echo "QC logs: $log_json, $log_html" | tee -a "$log_txt"
    echo "Threads: $threads" | tee -a "$log_txt"

    fastp_wrap &>> "$log_txt"

    echo "fastp completed successfully." | tee -a "$log_txt"
}

parse_args() {
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        print_usage
        exit 0
    fi

    if [[ $# -lt 7 ]]; then
        echo "Error: Missing required arguments." >&2
        print_usage
        exit 1
    fi

    declare -g in_r1="$1"
    declare -g out_r1="$2"
    declare -g out_r2="$3"
    declare -g failed_out="$4"
    declare -g log_txt="$5"
    declare -g log_json="$6"
    declare -g log_html="$7"
    declare -g threads="${8:-16}"

    declare -g in_r2="${in_r1/_R1/_R2}"
    [[ -f "$in_r2" ]] || { echo "Error: R2 file '$in_r2' does not exist." >&2; exit 1; }
}

fastp_wrap() {
    fastp \
        --detect_adapter_for_pe \
        --disable_quality_filtering \
        --failed_out "$failed_out" \
        --in1 "$in_r1" \
        --in2 "$in_r2" \
        --json "$log_json" \
        --html "$log_html" \
        --out1 "$out_r1" \
        --out2 "$out_r2" \
        --thread "$threads"
}

main "$@"
  #!/usr/bin/env bash
  set -euo pipefail

  print_usage() {
      cat <<EOF
  USAGE: fastp-emseq-wrapper.sh <INPUT R1 FASTQ.GZ> <OUTPUT R1> <OUTPUT R2> <FAILED OUT> <LOG TXT> <LOG JSON> <LOG HTML> [THREADS]

  DESCRIPTION:
    fastp wrapper for EM-seq data using Snakemake-style explicit I/O.
    Default threads is 16 (max for fastp).
  EOF
  }

  main() {
      parse_args "$@"

      echo "Running fastp on: $in_r1 and $in_r2" | tee "$log_txt"
      echo "Output files: $out_r1, $out_r2, $failed_out" | tee -a "$log_txt"
      echo "QC logs: $log_json, $log_html" | tee -a "$log_txt"
      echo "Threads: $threads" | tee -a "$log_txt"

      fastp_wrap &>> "$log_txt"

      echo "fastp completed successfully." | tee -a "$log_txt"
  }

  parse_args() {
      if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
          print_usage
          exit 0
      fi

      if [[ $# -lt 7 ]]; then
          echo "Error: Missing required arguments." >&2
          print_usage
          exit 1
      fi

      declare -g in_r1="$1"
      declare -g out_r1="$2"
      declare -g out_r2="$3"
      declare -g failed_out="$4"
      declare -g log_txt="$5"
      declare -g log_json="$6"
      declare -g log_html="$7"
      declare -g threads="${8:-16}"

      declare -g in_r2="${in_r1/_R1/_R2}"
      [[ -f "$in_r2" ]] || { echo "Error: R2 file '$in_r2' does not exist." >&2; exit 1; }
  }

  fastp_wrap() {
      fastp \
          --detect_adapter_for_pe \
          --disable_quality_filtering \
          --failed_out "$failed_out" \
          --in1 "$in_r1" \
          --in2 "$in_r2" \
          --json "$log_json" \
          --html "$log_html" \
          --out1 "$out_r1" \
          --out2 "$out_r2" \
          --thread "$threads"
  }

  main "$@"
