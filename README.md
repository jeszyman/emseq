The EM-seq repository contains modular workflows intended to be run from within a over-wrapping snakemake workflow.  

Current stable version tested with minimal in-repo example is tagged emseq.v3.1.0.  

![img](resources/test_smk.png)  


# Continuous Integration

[\![test-data](![img](https://github.com/jeszyman/emseq/actions/workflows/test-data.yml/badge.svg?branch=main))](<https://github.com/jeszyman/emseq/actions/workflows/test-data.yml>)  
[\![ci-snakemake](![img](https://github.com/jeszyman/emseq/actions/workflows/ci-snakemake.yml/badge.svg?branch=main))](<https://github.com/jeszyman/emseq/actions/workflows/ci-snakemake.yml>)  


# Change Log

-   Development since last tag
-   wf/emseq/v3.1.0  
    -   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-19 Fri] </span></span> Added a first github workflow test
    -   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-19 Fri] </span></span> Robust annotation of methylkit outputs validated as rscript
    -   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-19 Fri] </span></span> Added mature github tests for building in-repo test data and smk\_dry
-   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-19 Fri] </span></span> Updated EM-seq main pipeline to wf/emseq/v3.0.0.  
    -   Includes in-repo small test data for a complete run of emseq.smk
    -   Includes test.smk wrapper and corresponding test.yaml for in-repo small test run
    -   emseq.smk expanded to include differential methylation from nested list map
    -   Many small fixes for consistent naming and run condition optimization
-   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-18 Thu] </span></span> Updated EM-seq main pipeline to wf/emseq/v2.0.0. Mainly improved and simplified variable naming.

