## Bioinformatics
- methylKit: `filterByCoverage` fails on `methylBaseDB` — always use `save.db=FALSE` in `tileMethylCounts`. Always filter NA chr/start/end rows from united methylBase before tiling.

## Environment
- Always use `mamba` (not `conda`) for env installs/updates/creates — conda solver is extremely slow for large envs.
- `bwa-mem2 index` on full hg38 requires ~92GB RAM for suffix array. Use `bwa index` (via `bwameth.py index`) on machines with <64GB RAM.
- Full EM-seq pipeline for 5 samples at production depth needs ~300GB+ working disk (inputs + trimmed + BAMs + refs).

## Workflow architecture
- Snakemake workflows follow modular pattern: reusable modules (`emseq.smk`, `emseq_analysis.smk`) included by project-specific wrappers (`aero.smk`, `test.smk`).
- Wrappers handle: config, ref downloads, input symlinks, rule all. Modules handle: indexing, alignment, QC, methylation calling.
- `emseq_analysis.smk` and other tangled files are auto-generated from `emseq.org` — edits must be reflected in org source.
- `tools/setup_repos.sh` clones external analysis repos (mHapTools, wgbs_tools, UXM_deconv).

## Data
- AERO FASTQ data lives in GCS bucket `chaudhuri-lab-bucket1`. Use `gsutil` to pull (not copy via gcsfuse mount).
- AERO_24, AERO_25, AERO_101 in flowcell dirs under `20250224_LH00386_0211_A22MWJWLT4_06_19_25/` and `20250224_LH00386_0212_B22MWJ5LT4_fastq/`.
- NH12, NH13 in `20250224_LH00386_0212_B22MWJ5LT4_fastq/`.
