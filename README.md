
# Table of Contents

1.  [Prerequisites](#orged28da9)
2.  [Change Log](#orga4121ee)

<README.md>  

    python3 ~/repos/basecamp/scripts/emacs_export_header_to_markdown.py --org_file ~/repos/emseq/emseq.org --node_id bda70cff-0713-4e32-8da1-ee83924b8f00

The EM-seq repository contains modular workflows intended to be run from within a over-wrapping snakemake workflow.  

Current stable version tested with minimal in-repo example is tagged emseq.v3.0.0.  

![img](resources/test_smk.png)  


<a id="orged28da9"></a>

# Prerequisites


<a id="orga4121ee"></a>

# Change Log

-   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-19 Fri] </span></span> Updated EM-seq main pipeline to wf/emseq/v3.0.0.  
    -   Includes in-repo small test data for a complete run of emseq.smk
    -   Includes test.smk wrapper and corresponding test.yaml for in-repo small test run
    -   emseq.smk expanded to include differential methylation from nested list map
    -   Many small fixes for consistent naming and run condition optimization
-   <span class="timestamp-wrapper"><span class="timestamp">[2025-09-18 Thu] </span></span> Updated EM-seq main pipeline to wf/emseq/v2.0.0. Mainly improved and simplified variable naming.

