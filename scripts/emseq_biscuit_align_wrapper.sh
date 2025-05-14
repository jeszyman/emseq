#!/usr/bin/env bash
set -euo pipefail

print_usage() {
    cat <<EOF
USAGE: biscuit_align_wrapper.sh <R1 FASTQ.GZ> <BISCUIT REF FASTA> <OUTPUT BAM> <LOG DIR> [THREADS]

DESCRIPTION:
  Wrapper for Biscuit alignment of paired-end EM-seq data.
  Produces a sorted BAM file.
EOF
}

main() {
    parse_args "$@"

    echo "Running biscuit align on: $in_r1 and $in_r2" | tee "$log"
    echo "Reference genome: $biscuit_fa" | tee -a "$log"
    echo "Output BAM: $out_bam" | tee -a "$log"
    echo "Threads: $threads" | tee -a "$log"

    biscuit_align

    echo "Biscuit alignment completed successfully." | tee -a "$log"
}

parse_args() {
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        print_usage
        exit 0
    fi

    if [[ $# -lt 4 ]]; then
        echo "Error: Missing required arguments." >&2
        print_usage
        exit 1
    fi

    declare -g in_r1="$1"
    declare -g biscuit_fa="$2"
    declare -g out_bam="$3"
    declare -g log_dir="$4"
    declare -g threads="${5:-20}"

    [[ -f "$in_r1" ]] || { echo "Error: R1 file '$in_r1' does not exist." >&2; exit 1; }
    [[ -f "$biscuit_fa" ]] || { echo "Error: Reference genome '$biscuit_fa' not found." >&2; exit 1; }

    in_r2="${in_r1/_R1/_R2}"
    declare -g in_r2
    [[ -f "$in_r2" ]] || { echo "Error: R2 file '$in_r2' does not exist." >&2; exit 1; }

    base=$(basename "${in_r1%%_R1*}")
    declare -g base
    declare -g log="${log_dir}/${base}-biscuit-align.log"

    mkdir -p "$log_dir"
}

biscuit_align() {
    biscuit align \
        -@ "$threads" \
        -biscuit-ref "$biscuit_fa" \
        "$in_r1" "$in_r2" \
        | samtools sort -@ "$threads" -o "$out_bam" &>> "$log"
}

main "$@"
  #!/usr/bin/env bash
  set -euo pipefail

  print_usage() {
      cat <<EOF
  USAGE: biscuit_align_wrapper.sh <R1 FASTQ.GZ> <BISCUIT REF FASTA> <OUTPUT BAM> <LOG DIR> [THREADS]

  DESCRIPTION:
    Wrapper for Biscuit alignment of paired-end EM-seq data.
    Produces a sorted BAM file.
  EOF
  }

  main() {
      parse_args "$@"

      echo "Running biscuit align on: $in_r1 and $in_r2" | tee "$log"
      echo "Reference genome: $biscuit_fa" | tee -a "$log"
      echo "Output BAM: $out_bam" | tee -a "$log"
      echo "Threads: $threads" | tee -a "$log"

      biscuit_align

      echo "Biscuit alignment completed successfully." | tee -a "$log"
  }

  parse_args() {
      if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
          print_usage
          exit 0
      fi

      if [[ $# -lt 4 ]]; then
          echo "Error: Missing required arguments." >&2
          print_usage
          exit 1
      fi

      declare -g in_r1="$1"
      declare -g biscuit_fa="$2"
      declare -g out_bam="$3"
      declare -g log_dir="$4"
      declare -g threads="${5:-20}"

      [[ -f "$in_r1" ]] || { echo "Error: R1 file '$in_r1' does not exist." >&2; exit 1; }
      [[ -f "$biscuit_fa" ]] || { echo "Error: Reference genome '$biscuit_fa' not found." >&2; exit 1; }

      in_r2="${in_r1/_R1/_R2}"
      declare -g in_r2
      [[ -f "$in_r2" ]] || { echo "Error: R2 file '$in_r2' does not exist." >&2; exit 1; }

      base=$(basename "${in_r1%%_R1*}")
      declare -g base
      declare -g log="${log_dir}/${base}-biscuit-align.log"

      mkdir -p "$log_dir"
  }

  biscuit_align() {
      biscuit align \
          -@ "$threads" \
          -biscuit-ref "$biscuit_fa" \
          "$in_r1" "$in_r2" \
          | samtools sort -@ "$threads" -o "$out_bam" &>> "$log"
  }

  main "$@"
