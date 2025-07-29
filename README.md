The EM-seq repository contains modular workflows intended to be run from within a over-wrapping snakemake workflow.  

A master data dir, emseq\_dir directs outputs  
Outside of emseq\_dir, BISCUIT index is genrated in data\_dir/ref/biscuit  

Adapter and quality trimming was performed using fastp without any fixed-length end trimming. Alignment used biscuit with Ensemble hg38 primary assembly (<https:/ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz>). Aligned BAMs were deduplicated using dupsifter with mate-tag annotation. Methylation pileups were generated from deduplicated BAMs using biscuit pileup, and the resulting VCFs were converted to bismark-formatted bed files. Bismark-style BEDs were converted to methylKit objects using a custom R script, and downstream differential methylation analysis was performed in methylKit.  

Copy number analysis was performed with ichorCNA. Deduplicated BAMs were used as input. Read counts were generated in 1 Mb windows using readCounter, filtering for base quality and including autosomes and sex chromosomes.  

Copy number estimation was performed using ichorCNA, specifying hg38 reference annotations, GC and mappability corrections, and the default panel of normals. The model estimated tumor fraction, ploidy, and subclonal prevalence using autosomal chromosomes for training. ichorCNA was run with package-specified low input settings including reduced copy number states, no subclonal events and initial ploidy set at diploid.  

Directory Structure  

    -- Main
       |-- config
       |-- resources
       |-- results
       |-- scripts
       |-- test
       |-- tools
       |-- workflows


# Prerequisites


# Change Log

