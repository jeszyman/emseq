## Bioinformatics
- methylKit: `filterByCoverage` fails on `methylBaseDB` — always use `save.db=FALSE` in `tileMethylCounts`. Always filter NA chr/start/end rows from united methylBase before tiling.
- methylKit DMR calling is restricted to **two-group** treatment: `make_methylkit_unite_db.R` and `make_methylkit_diff_tiled_db.R` `stop()` if `tx` has ≠2 unique values. For >2 groups run pairwise comparisons; for continuous/survival outcomes binarize (or use Cox / DSS general design). methylKit's `meth.diff` (max−min of per-group means) and numeric-trend test are unreliable for >2 or continuous treatment.
- Both methylKit scripts run `normalizeCoverage()` before unite/tile, and stage inputs as symlinks into a per-experiment scratch dir (`.scratch_<suffix>/`) so destrand/normalize intermediates (written next to each input, keyed by library only) don't collide across concurrent experiments.
- Data-integrity helpers: `scripts/emseq_checks.R` — `log_n`, `check_n`, `check_nonempty`, `assert_cols`; source in analysis scripts to guard joins/filters against silent data loss.

## Environment
- Always use `mamba` (not `conda`) for env installs/updates/creates — conda solver is extremely slow for large envs.
- `bwa-mem2 index` on full hg38 requires ~92GB RAM for suffix array. Use `bwa index` (via `bwameth.py index`) on machines with <64GB RAM.
- Full EM-seq pipeline for 5 samples at production depth needs ~300GB+ working disk (inputs + trimmed + BAMs + refs).

## Workflow architecture
- Snakemake workflows follow modular pattern: reusable modules (`emseq.smk`, `emseq_analysis.smk`) included by project-specific wrappers (`aero.smk`, `test.smk`).
- Wrappers handle: config, ref downloads, input symlinks, rule all. Modules handle: indexing, alignment, QC, methylation calling.
- `emseq_analysis.smk` and other tangled files are auto-generated from `emseq.org` — edits must be reflected in org source.
- `tools/setup_repos.sh` clones external analysis repos (mHapTools, wgbs_tools, UXM_deconv).
- Snakemake runs every `shell:` in bash strict mode (`set -euo pipefail`) by default — rely on it; do NOT add a workflow `shell.prefix` (it *replaces* the default and silently drops pipefail/nounset). For a one-off samtools SIGPIPE (exit 141), scope `set +o pipefail` to that command.
- `emseq_filter_bam` prunes the BAM `@SQ` header to primary contigs (samtools reheader) so downstream methyldackel/mosdepth/stats see a primary-only BAM even with a decoy reference.
- conda env YAMLs keep their `name:` field (it's ignored under `--use-conda`, which builds at a content-hash prefix; prefix the name with the repo stem only if reusing as a module elsewhere).

## Data
- AERO FASTQ data lives in GCS bucket `chaudhuri-lab-bucket1`. Use `gsutil` to pull (not copy via gcsfuse mount).
- AERO_24, AERO_25, AERO_101 in flowcell dirs under `20250224_LH00386_0211_A22MWJWLT4_06_19_25/` and `20250224_LH00386_0212_B22MWJ5LT4_fastq/`.
- NH12, NH13 in `20250224_LH00386_0212_B22MWJ5LT4_fastq/`.
