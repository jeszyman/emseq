#!/usr/bin/env bash
# Read-end settings on a two-group cohort held outside the repository.
#
# extract: copies each library's filtered alignment (and index) from a bucket
#   folder, calls CpG methylation with MethylDackel under each read-end setting
#   (the flags of VALIDATION_END_SETTINGS in workflows/validation.smk), gzips
#   the calls, and deletes the alignment. Libraries at a time is the first
#   option, so the alignments on disk at once stay bounded.
# dmr: converts each call table to MethylDackel's methylKit table, makes the
#   per-library methylKit object, then calls differentially methylated tiles
#   between the groups once per setting with the pipeline's methylKit script.
#   Healthy libraries are named NH_*, patients any other name.
# compare: writes the comparison tables (emseq_validation_settings_compare.R),
#   naming libraries by the index map kept in the project repository.
#
# Usage:
#   emseq_validation_settings_cohort.sh extract <work dir> <bucket folder> <libraries at a time>
#   emseq_validation_settings_cohort.sh dmr <work dir> <conversions at a time> <cores per setting>
#   emseq_validation_settings_cohort.sh compare <work dir> <index map>
# The work dir holds ref/hg38.fa (+ .fai), libraries.txt (one library name per
# line, written from the bucket listing when absent), and receives extract/,
# tabix/, dmr/, global/, compare/, logs/, timing.tsv, and failed.txt.
# MethylDackel, pigz, gsutil, Rscript (with methylKit), and python3 are taken from PATH.
# Bucket alignments are named <library>.hg38.bwa_meth.coorsort.filt.bam.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
declare -A FLAGS=([none]="" [neb]="--nOT 0,0,0,5 --nOB 0,0,5,0" [all10]="--nOT 10,11,10,11 --nOB 10,11,10,11")
SETTINGS=(none neb all10)
mode=${1:?usage: see header}; W=${2:?work dir}
cd "$W" || exit 1
mkdir -p extract bams logs global tabix dmr tmp
export TMPDIR=$W/tmp
[[ -f timing.tsv ]] || printf 'unit\tstep\tseconds\n' > timing.tsv
stamp() { printf '%s\t%s\t%s\n' "$1" "$2" $(( $(date +%s) - $3 )) >> timing.tsv; }

extract_one() {
  local lib=$1 bucket=$2 bam=bams/$1.bam t0 s
  ls extract/"$lib".{none,neb,all10}_CpG.bedGraph.gz > /dev/null 2>&1 && return 0
  t0=$(date +%s)
  gsutil -q cp "$bucket/$lib.hg38.bwa_meth.coorsort.filt.bam" "$bam" \
    && gsutil -q cp "$bucket/$lib.hg38.bwa_meth.coorsort.filt.bam.bai" "$bam.bai" \
    || { echo "FAIL copy $lib" >> failed.txt; return 1; }
  stamp "$lib" copy "$t0"
  for s in "${SETTINGS[@]}"; do
    (
      t1=$(date +%s)
      # shellcheck disable=SC2086
      MethylDackel extract -@ 4 --keepDiscordant ${FLAGS[$s]} ref/hg38.fa "$bam" -o "extract/$lib.$s" > "logs/$lib.$s.log" 2>&1 \
        && pigz -p 2 -f "extract/$lib.${s}_CpG.bedGraph" \
        || echo "FAIL extract $lib $s" >> failed.txt
      stamp "$lib" "$s" "$t1"
    ) &
  done
  wait
  rm -f "$bam" "$bam.bai"
}

single_one() {
  local lib=$1 s=$2 tx t0; t0=$(date +%s)
  [[ $lib == NH_* ]] && tx=0 || tx=1
  [[ -s tabix/$s/$lib.txt.bgz.tbi ]] && return 0
  mkdir -p "meth/$s" "tabix/$s"
  python3 "$HERE/emseq_bedgraph_to_methylkit.py" "extract/$lib.${s}_CpG.bedGraph.gz" ref/hg38.fa "meth/$s/$lib.methylKit" 2> "global/$lib.$s.txt" \
    && Rscript "$HERE/make_single_amp_methylkit_obj.R" --amp_file "meth/$s/$lib.methylKit" --library_id "$lib" --treatment "$tx" \
         --mincov 1 --build hg38 --out_dir "tabix/$s" > "logs/$lib.$s.single.log" 2>&1 \
    || echo "FAIL single $lib $s" >> failed.txt
  # the per-library step leaves the uncompressed table beside the tabix file
  rm -f "meth/$s/$lib.methylKit" "tabix/$s/$lib.txt"
  stamp "$lib.$s" single "$t0"
}

dmr_one() {
  local s=$1 cores=$2 t0 dbs ids txs; t0=$(date +%s)
  dbs=$(sed "s|.*|tabix/$s/&.txt.bgz|" libraries.txt | tr '\n' ' ')
  ids=$(tr '\n' ' ' < libraries.txt)
  txs=$(sed 's/^NH_.*/0/; t; s/.*/1/' libraries.txt | tr '\n' ' ')
  Rscript "$HERE/make_methylkit_diff_tiled_db.R" --lib_db_list "$dbs" --lib_id_list "$ids" --treatment_list "$txs" \
    --cores "$cores" --out_dir "dmr/$s" --suffix 2l_vs_nh --assembly hg38 --mincov 10 --win_size 10000 \
    --chunk_size 1e9 --min_per_group 3 > "logs/dmr.$s.log" 2>&1 || echo "FAIL dmr $s" >> failed.txt
  stamp "$s" dmr "$t0"
}

case $mode in
  extract)
    bucket=${3:?bucket folder}; at_once=${4:?libraries at a time}
    [[ -s libraries.txt ]] || gsutil ls "$bucket/" | grep '\.bam$' | sed 's|.*/||; s|\.hg38\..*||' > libraries.txt
    n=0
    while read -r lib; do
      extract_one "$lib" "$bucket" &
      n=$((n + 1)); (( n % at_once == 0 )) && wait
    done < libraries.txt
    wait ;;
  dmr)
    at_once=${3:?conversions at a time}; cores=${4:?cores per setting}
    n=0
    for s in "${SETTINGS[@]}"; do
      while read -r lib; do
        single_one "$lib" "$s" &
        n=$((n + 1)); (( n % at_once == 0 )) && wait
      done < libraries.txt
    done
    wait
    # methylKit normalizes one library at a time, so the settings run side by side
    for s in "${SETTINGS[@]}"; do dmr_one "$s" "$cores" & done
    wait ;;
  compare)
    Rscript "$HERE/emseq_validation_settings_compare.R" . compare "${3:?index map}" ;;
  *) echo "unknown mode $mode" >&2; exit 2 ;;
esac
[[ -s failed.txt ]] && exit 1
echo "done $mode" > "finished_$mode"
