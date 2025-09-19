#!/usr/bin/env bash
# setup_test_data.sh — chr22 subset + lambda + pUC19 FASTA, plus 4 tiny paired WGBS FASTQs
# Paired-only; hard-fails if a run lacks _1/_2. Cleans tests/full/ before writing.
set -euo pipefail

# --- self-locate & cd to repo root (script is assumed to live one level below root, e.g., tools/) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
cd "$REPO_DIR" || { echo "ERR: failed to cd to repo root: $REPO_DIR" >&2; exit 1; }

# --- config ---
R_EMSEQ="${R_EMSEQ:-$PWD}"
TEST_DIR="${R_EMSEQ}/tests/full"
OUT_DIR="${TEST_DIR}/inputs"
OUT_FA="${OUT_DIR}/chr22.test.fa.gz"
OUT_LAMBDA="${OUT_DIR}/lambda.fa.gz"
OUT_PUC19="${OUT_DIR}/pUC19.fa.gz"
OUT_BLK="${OUT_DIR}/hg38-blacklist.v2.bed.gz"

# NEW: bed outputs
KEEP_BED="${OUT_DIR}/chr22.keep.bed"
EXCL_BED="${OUT_DIR}/chr22.exclude.blacklist.bed.gz"

REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
FA_HEAD_LINES=4000000

# full (small) references
LAMBDA_URL="https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=9626243&db=nuccore&report=fasta&retmode=text"

# pUC19 via NCBI efetch (hard-coded tool/email)
PUC19_EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137.2&rettype=fasta&retmode=text&tool=emseq_setup&email=anon@example.com"

BLK_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"

NREADS="${NREADS:-60000}"   # reads per mate for FASTQ clipping
FQ_HEAD_LINES=$((NREADS * 4))

# ENA WGBS runs (paired) — adjust as needed
ACCESSIONS=( "ERR022484" "ERR022487" "ERR022003" "ERR022483" )

# --- deps ---
need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 1; }; }
need git; need wget; need curl; need zcat; need gzip; need wc; need head; need grep; need sed
# NEW: for beds
need samtools; need bedtools; need bgzip; need tabix; need awk; need sort; need cut

# --- ensure .gitignore rules: ignore tests/full/* but NOT tests/full/inputs/** ---
ensure_gitignore() {
  local gi="${REPO_DIR}/.gitignore"
  local req1='tests/full/*'
  local req2='!tests/full/inputs/'
  local req3='!tests/full/inputs/**'
  touch "$gi"
  local have1=0 have2=0 have3=0
  grep -Fxq "$req1" "$gi" && have1=1 || true
  grep -Fxq "$req2" "$gi" && have2=1 || true
  grep -Fxq "$req3" "$gi" && have3=1 || true
  if (( ! have1 || ! have2 || ! have3 )); then
    {
      echo ''
      echo '# --- test data (auto-managed by setup_test_data.sh) ---'
      (( have1 )) || echo "$req1"
      (( have2 )) || echo "$req2"
      (( have3 )) || echo "$req3"
    } >> "$gi"
    echo "[gitignore] ensured patterns for tests/full with inputs preserved"
  fi
}

# ENA dir that contains both _1/_2
ena_dir_for() {
  local acc="$1" first6="${acc:0:6}" last3="${acc: -3}"
  local candidates=(
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${acc}/"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${last3}/${acc}/"
  )
  for d in "${candidates[@]}"; do
    if wget -q --spider "${d}${acc}_1.fastq.gz" && wget -q --spider "${d}${acc}_2.fastq.gz"; then
      echo "$d"; return 0
    fi
  done
  return 1
}

# wget stream → zcat → head → gzip (silence SIGPIPE noise)
fetch_head() {
  ( set +o pipefail
    wget -qO- "$1" | zcat 2>/dev/null | head -n "$2" | gzip > "$3"
  ) || true
}

# download full compressed file (blacklist)
fetch_all_gz() { wget -qO "$2" "$1"; }

# --- add near fetch_puc19() ---
fetch_plain_fasta_gzip() {
  local url out tmp
  url="$1"
  out="$2"
  tmp="${out%.gz}.tmp"
  rm -f "$tmp" "$out"
  curl -fsSL "$url" -o "$tmp"
  [[ -s "$tmp" ]] || { echo "ERR: lambda fetch produced empty file"; exit 1; }
  head -n1 "$tmp" | grep -q '^>' || { echo "ERR: lambda FASTA header missing"; exit 1; }
  gzip -f "$tmp"
  mv "${tmp}.gz" "$out"
}

# robust pUC19: fetch plain FASTA to tmp with curl, validate header, then gzip
fetch_puc19() {
  local url out tmp
  url="$1"; out="$2"; tmp="${out%.gz}.tmp"
  rm -f "$tmp" "$out"
  curl -fsSL "$url" -o "$tmp"
  [[ -s "$tmp" ]] || { echo "ERR: pUC19 fetch produced empty file"; exit 1; }
  head -n1 "$tmp" | grep -q '^>' || { echo "ERR: pUC19 FASTA header missing"; exit 1; }
  gzip -f "$tmp"
  mv "${tmp}.gz" "$out"
}

check_fastq() {
  local f="$1"
  [[ -s "$f" ]] || { echo "ERR: empty $f" >&2; return 1; }
  local n; n=$(zcat "$f" 2>/dev/null | wc -l)
  (( n % 4 == 0 )) || echo "WARN: $f has $n lines (not multiple of 4)"
  echo "[ok] $(basename "$f"): $(( n / 4 )) reads"
}

# --- run ---
ensure_gitignore

echo "[clean] removing ${TEST_DIR}"
rm -rf "${TEST_DIR}"
mkdir -p "${OUT_DIR}"

# keep inputs dir tracked even if empty (optional helper file)
[[ -e "${OUT_DIR}/.gitkeep" ]] || : > "${OUT_DIR}/.gitkeep"

echo "[ref] chr22 subset → ${OUT_FA}"
fetch_head "${REF_URL}" "${FA_HEAD_LINES}" "${OUT_FA}"
[[ -s "${OUT_FA}" ]] || { echo "ERR: failed to write ${OUT_FA}"; exit 1; }

echo "[ref] lambda (full) → ${OUT_LAMBDA}"
fetch_plain_fasta_gzip "${LAMBDA_URL}" "${OUT_LAMBDA}"
[[ -s "${OUT_LAMBDA}" ]] || { echo "ERR: failed to write ${OUT_LAMBDA}"; exit 1; }

echo "[ref] pUC19 (NCBI efetch) → ${OUT_PUC19}"
fetch_puc19 "${PUC19_EFETCH}" "${OUT_PUC19}"
[[ -s "${OUT_PUC19}" ]] || { echo "ERR: failed to write ${OUT_PUC19}"; exit 1; }

echo "[ref] hg38 blacklist → ${OUT_BLK}"
fetch_all_gz "${BLK_URL}" "${OUT_BLK}"
[[ -s "${OUT_BLK}" ]] || { echo "ERR: failed to write ${OUT_BLK}"; exit 1; }


# --- tiny paired FASTQs ---
i=1
for acc in "${ACCESSIONS[@]}"; do
  dir="$(ena_dir_for "$acc")" || { echo "ERR: ${acc} is not paired on ENA" >&2; exit 1; }
  r1="${dir}${acc}_1.fastq.gz"
  r2="${dir}${acc}_2.fastq.gz"

  id=$(printf "lib%03d" "$i")
  out1="${OUT_DIR}/${id}.raw_R1.fastq.gz"
  out2="${OUT_DIR}/${id}.raw_R2.fastq.gz"
  echo "[fq] ${id} (${acc}) → ${out1}, ${out2}"

  fetch_head "$r1" "${FQ_HEAD_LINES}" "$out1"
  fetch_head "$r2" "${FQ_HEAD_LINES}" "$out2"
  check_fastq "$out1"
  check_fastq "$out2"
  i=$((i+1))
done

# --- NEW: build keep/exclude BEDs in inputs/ from chr22.test.fa.gz + blacklist ---
echo "[beds] building keep/exclude from ${OUT_FA}"
FA_UNGZ="${OUT_DIR}/chr22.test.fa"
zcat "${OUT_FA}" > "${FA_UNGZ}"
samtools faidx "${FA_UNGZ}"

# KEEP: spans for all contigs present in the .fai (chr22-only here)
awk 'BEGIN{OFS="\t"} {print $1,0,$2}' "${FA_UNGZ}.fai" \
  | sort -k1,1 -k2,2n > "${KEEP_BED}" |

# EXCLUDE: clip blacklist to contigs present in .fai and bgzip
awk 'NR==FNR{ok[$1]=1; next} ok[$1]' <(cut -f1 "${FA_UNGZ}.fai") <(zcat "${OUT_BLK}") \
  | sort -k1,1 -k2,2n | bgzip > "${EXCL_BED}" |
tabix -p bed "${EXCL_BED}" || true
echo "[beds] keep=${KEEP_BED}  exclude=${EXCL_BED}"

echo "Done. Outputs in: ${OUT_DIR}"
echo "  - FASTA subset: ${OUT_FA} (+ ${FA_UNGZ} + .fai)"
echo "  - Keep BED:     ${KEEP_BED}"
echo "  - Exclude BED:  ${EXCL_BED}"
