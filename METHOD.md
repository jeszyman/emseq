
# Introduction

Enzymatic methyl sequencing (EM-seq) measures cytosine methylation at single-base resolution without bisulfite. The TET2 enzyme oxidizes 5-methylcytosine and 5-hydroxymethylcytosine, and the APOBEC3A deaminase then converts the unprotected cytosines to uracil, which is copied as thymine during amplification (<a href="#citeproc_bib_item_12">Vaisvila et al. 2021</a>). A methylated cytosine therefore reads as C and an unmethylated cytosine as T, as in a bisulfite library, and EM-seq reads are analyzed with the aligners and callers written for bisulfite sequencing (<a href="#citeproc_bib_item_10">New England Biolabs 2023</a>). The enzymatic conversion works from picograms of DNA (<a href="#citeproc_bib_item_12">Vaisvila et al. 2021</a>), which suits plasma cell-free DNA (cfDNA), where a blood draw yields nanograms.  

This pipeline takes paired-end EM-seq reads to per-CpG methylation calls with their quality measures, and from the calls to differential methylation, methylation haplotypes, and tissue-of-origin fractions. It was written for plasma cfDNA. Plasma fragments are short, with a modal length near 167 bp, and they are cut at preferred positions and carry single-stranded overhangs. Several steps that are routine for sheared genomic DNA behave differently on such fragments, and the pipeline's departures from common practice are stated with the measurement behind them under [Design decisions](#design-decisions).  

This report first follows the data through the pipeline ([Figure 1](#fig-core-dataflow)), from trimming and alignment to methylation calling, the spike-in controls, and quality control, and then describes the downstream analyses and the output files. The closing sections give the design decisions, the limitations, and the prior work that the pipeline draws on.  


# Overview and data flow

The input is paired-end reads from an EM-seq library, as one pair of FASTQ files per library. The outputs are a per-CpG table of methylated and unmethylated calls for each library, quality-control reports, and cohort tables of differential methylation, methylation haplotypes, and tissue-of-origin fractions.  

<a id="fig-core-dataflow"></a>
![**The core processing rules and their dependencies.** Rule graph of `workflows/emseq.smk` as drawn by Snakemake for the in-repository test libraries, which are aligned with both bwa-meth and BISCUIT.](resources/test_smk.png)
*Figure 1: **The core processing rules and their dependencies.** Rule graph of `workflows/emseq.smk` as drawn by Snakemake for the in-repository test libraries, which are aligned with both bwa-meth and BISCUIT.*  

The stages of [Figure 1](#fig-core-dataflow) are listed in order with the Snakemake rules that run them. Core processing is in [`emseq.smk`](workflows/emseq.smk), downstream analysis is in [`emseq_analysis.smk`](workflows/emseq_analysis.smk), and a project wrapper sets the configuration and includes both (<a href="#citeproc_bib_item_9">Mölder et al. 2021</a>).  

1.  **Reference indexing**, run once per reference (`emseq_bwa_meth_index`, `emseq_biscuit_index`). The reference FASTA is converted in silico and indexed for the chosen aligner.
2.  **Trimming** (`emseq_fastp`). Adapter sequence is removed from both mates.
3.  **Alignment** (`emseq_align_bwameth`, `emseq_align_biscuit`). Trimmed pairs are aligned to the converted reference and sorted by coordinate.
4.  **Duplicate marking** (`emseq_dedup`). Pairs that pass the pair test are kept, and duplicates among them are flagged.
5.  **Filtering** (`emseq_filter_bam`). Flagged duplicates, low mapping quality, and reads outside the kept regions are removed.
6.  **Methylation calling** (`emseq_methyldackel`, `emseq_make_single_methylkit_amp_obj`). Methylated and unmethylated calls are counted at each CpG and stored as an indexed table.
7.  **Spike-in controls** (`emseq_align_bwameth_spike`, `emseq_methyldackel_spike`). The trimmed reads are aligned to the lambda and pUC19 references, and control methylation is counted.
8.  **Quality control** (`emseq_fastqc`, `emseq_mosdepth`, `emseq_mbias`, `emseq_samtools_stats`, `emseq_multiqc`). Read quality, coverage, methylation by read position, and alignment statistics are collected into one report.
9.  **Downstream analysis** (the rules of `emseq_analysis.smk`). Libraries are compared between two groups, read-level methylation patterns are summarized, and tissue-of-origin fractions are estimated.


# Trimming

Adapters are removed with [fastp](https://github.com/OpenGene/fastp) (<a href="#citeproc_bib_item_3">Chen et al. 2018</a>) with adapter detection for paired-end reads switched on (`--detect_adapter_for_pe`). Adapter removal matters more in a converted library than in an ordinary one. The EM-seq adaptor is synthesized with methylated cytosines so that it survives conversion, and an untrimmed adaptor base therefore enters the calls as a methylated cytosine at the 3' end of the read. Plasma fragments are often shorter than the read, so a large share of plasma reads run into the adaptor.  

Read-level quality filtering is switched off (`--disable_quality_filtering`), and no fixed number of bases is clipped from either end of either mate. Low-quality bases are instead excluded when methylation is counted, by the base-quality threshold of the caller, so that they still contribute to alignment. A read shorter than 15 bases after adapter removal is discarded, which is the fastp default, and discarded reads are written to a separate file. The reason no bases are clipped is given under [Read-end handling](#read-end-handling).  


# Alignment

Reads are aligned with [bwa-meth](https://github.com/brentp/bwa-meth), which converts the reference and the reads in silico (C to T in read 1 and G to A in read 2), aligns the converted reads with BWA-MEM, and restores the original read sequence in the output so that the caller can read the methylation state. The reference index is built once per reference with `bwameth.py index-mem2`. [BISCUIT](https://github.com/huishenlab/biscuit) (<a href="#citeproc_bib_item_14">Zhou et al. 2024</a>) is available as a second aligner, and the wrapper variable `emseq_align_methods` sets which aligners run. Alignments are sorted by coordinate with samtools (<a href="#citeproc_bib_item_6">Li 2011</a>). A library can be aligned to more than one reference in one run, and the reference name is part of each output file name.  


# Duplicate marking and filtering

Read pairs are first tested for a plausible fragment. A pair is kept when both mates are mapped to the same chromosome, face each other, and span a fragment of at most 1,000 bp. The aligner's proper-pair flag is left unused here, for the reason given under [Pairs without the proper-pair flag](#pairs-without-the-proper-pair-flag). Duplicates among the kept pairs are flagged with [dupsifter](https://github.com/huishenlab/dupsifter).  

Duplicates are flagged in one rule and removed in the next. The filtering rule keeps pairs that pass the same pair test with a mapping quality of 30 or more, and removes secondary and supplementary alignments, alignments flagged as failing quality control, and flagged duplicates (`samtools view -q 30 -F 3852`). It then removes reads that overlap a blacklist BED and keeps reads within a keep BED, such as the autosomes and sex chromosomes of GRCh38. Last, the sequence dictionary of the BAM header is cut down to the primary chromosomes, so that the caller and the coverage tools see a primary-only alignment when the reference carries decoy and unplaced contigs.  


# Methylation calling

Methylation is counted with [MethylDackel](https://github.com/dpryan79/MethylDackel) (`extract --methylKit --keepDiscordant`, version 0.6.1) on the filtered alignment. A call is the base a read shows at a reference CpG: C for methylated and T for unmethylated on a read from the original top strand, and G and A on a read from the original bottom strand. MethylDackel counts a base when its quality is 5 or more and its alignment has a mapping quality of 10 or more, which the filter has already raised to 30. Where the two mates of a pair overlap, the position is counted once for the pair. The option `--keepDiscordant` makes MethylDackel count pairs whose proper-pair flag is unset, which it otherwise ignores.  

The first five sequenced bases of read 2, which lie at the filled-in 3' end of the original strand, are excluded from the counts with the options of the NEB EM-seq pipeline (`--nOT 0,0,0,5 --nOB 0,0,5,0`), and the same bases are excluded on the spike-in controls. The measurement behind that choice is under [Read-end handling](#read-end-handling). The per-CpG counts are then loaded into [methylKit](https://github.com/al2na/methylKit) (<a href="#citeproc_bib_item_1">Akalin et al. 2012</a>) with a minimum coverage set in the configuration (`emseq-mincov`) and written as a tabix-indexed table, one per library.  


# Spike-in controls

The EM-seq kit supplies two control DNAs that are added to each sample before library preparation: unmethylated lambda DNA, and pUC19 DNA methylated at CpG sites. Lambda measures conversion. A cytosine that reads as methylated in lambda escaped conversion, so lambda methylation is the non-conversion rate of the library, and it is expected to be below 1%. pUC19 measures protection. A CpG that reads as unmethylated in pUC19 was converted despite its methylation, and pUC19 CpG methylation is expected to be above 90% (<a href="#citeproc_bib_item_12">Vaisvila et al. 2021</a>; <a href="#citeproc_bib_item_2">Chauhan et al. 2024</a>; <a href="#citeproc_bib_item_10">New England Biolabs 2023</a>).  

The trimmed reads of each library are aligned separately to the lambda and pUC19 references with bwa-meth, mapped reads are kept without duplicate marking, and methylation is counted with the same MethylDackel command as for the sample. A library that was prepared without a control still gives alignments to that control's reference, because short stretches of sample reads match it. bwa-meth flags most of these partial alignments as failing quality control and MethylDackel leaves them out, but some pass. In 8 plasma libraries prepared without pUC19, 240 reads in total were counted on the pUC19 reference, with 4.8% CpG methylation over 651 calls. A pUC19 table with a few hundred calls and low methylation therefore indicates a control that was never added.  


# Quality control

Read quality is reported before and after trimming with FastQC, and fastp writes its own trimming report. Coverage is profiled on the filtered alignment with [mosdepth](https://github.com/brentp/mosdepth) (<a href="#citeproc_bib_item_11">Pedersen and Quinlan 2018</a>), which gives the depth distribution, the depth within the kept regions, and the share of bases at or above each depth threshold set in the configuration. Alignment statistics come from `samtools stats` and `samtools flagstat` on the filtered alignment.  

Methylation by read position is tabulated with `MethylDackel mbias` on the duplicate-marked alignment, before the mapping-quality and region filters. MethylDackel leaves out flagged duplicates by default, so the table excludes them as the calls do. The table gives, for each strand and mate, the methylated and unmethylated calls at each position along the read, which is the form in which a position-dependent bias is looked for in converted libraries (<a href="#citeproc_bib_item_5">Hansen, Langmead, and Irizarry 2012</a>). [MultiQC](https://github.com/MultiQC/MultiQC) collects the reports of all libraries into one page.  


# Downstream analysis

Differential methylation is called with methylKit between the two groups of an experiment, which the configuration defines as a list of libraries with a group label for each (`meth-map`). Per-library coverage is normalized, the libraries are united over the CpGs covered in at least a set number of libraries per group, and a difference is tested at each CpG and in windows of a set width. An experiment must have exactly two groups, and the scripts stop otherwise, because the difference that methylKit reports is unreliable with more than two groups. Significant CpGs and windows are annotated with [annotatr](https://github.com/rcavalcante/annotatr) for CpG islands, shores, shelves, promoters, exons, and introns.  

Methylation haplotypes describe the methylation of consecutive CpGs on single molecules. The filtered alignment is converted to the mhap format with [mHapTools](https://github.com/butyuhao/mHapTools), and the methylated haplotype load, the proportion of discordant reads, and the methylation entropy are computed for each methylation haplotype block (<a href="#citeproc_bib_item_4">Guo et al. 2017</a>) with mhaptk.  

Tissue-of-origin fractions are estimated with UXM against a reference atlas of cell-type methylation (<a href="#citeproc_bib_item_7">Loyfer et al. 2023</a>). Contig names are first rewritten to the form wgbstools expects, and per-read methylation patterns are written with `wgbstools bam2pat` (<a href="#citeproc_bib_item_8">Loyfer, Rosenski, and Kaplan 2026</a>). The in-repository test covers differential methylation and annotation. The haplotype and deconvolution rules need compiled external tools and a genome-wide reference, and they are run on full data sets outside the test.  


# Pipeline outputs

The files a library or an experiment produces are listed in [Table 1](#tab-outputs). Paths are relative to the data directory set in the configuration, and the placeholders stand for the library, the reference, the aligner, and the experiment.  

<table id="tab-outputs" border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
<colgroup>
<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">File</th>
<th scope="col" class="org-left">Content</th>
</tr>
</thead>
<tbody>
<tr>
<td class="org-left"><code>emseq/fastqs/{lib}.trimmed_{R1,R2}.fastq.gz</code></td>
<td class="org-left">Adapter-trimmed reads</td>
</tr>

<tr>
<td class="org-left"><code>emseq/bams/{lib}.{ref}.{aligner}.coorsort.deduped.bam</code></td>
<td class="org-left">Pairs that pass the pair test, with duplicates flagged</td>
</tr>

<tr>
<td class="org-left"><code>emseq/bams/{lib}.{ref}.{aligner}.coorsort.filt.bam</code></td>
<td class="org-left">Filtered alignment that the caller reads</td>
</tr>

<tr>
<td class="org-left"><code>emseq/meth/{lib}.{ref}.{aligner}_methyldackel_CpG.methylKit</code></td>
<td class="org-left">Per-CpG coverage and percent methylated, as MethylDackel writes it</td>
</tr>

<tr>
<td class="org-left"><code>emseq/dmr/tabix/{lib}.{ref}.{aligner}.methyldackel.txt.bgz</code></td>
<td class="org-left">The same calls as a tabix-indexed methylKit table</td>
</tr>

<tr>
<td class="org-left"><code>emseq/spike/{lib}.{control}.{aligner}_methyldackel_CpG.methylKit</code></td>
<td class="org-left">Per-CpG calls on the lambda or pUC19 reference</td>
</tr>

<tr>
<td class="org-left"><code>emseq/qc/{lib}.{ref}.{aligner}_emseq_mbias.txt</code></td>
<td class="org-left">Methylated and unmethylated calls by strand, mate, and read position</td>
</tr>

<tr>
<td class="org-left"><code>emseq/qc/mosdepth_{lib}.{ref}.{aligner}.mosdepth.summary.txt</code></td>
<td class="org-left">Depth summary, with distribution, threshold, and quantized tracks beside it</td>
</tr>

<tr>
<td class="org-left"><code>emseq/qc/multiqc.html</code></td>
<td class="org-left">Collected quality report for all libraries</td>
</tr>

<tr>
<td class="org-left"><code>emseq/dmr/diff/methylDiff_{experiment}.txt.bgz</code></td>
<td class="org-left">Per-CpG differences between the two groups</td>
</tr>

<tr>
<td class="org-left"><code>emseq/dmr/diff/methylBase_{experiment}.tiled.txt.bgz</code></td>
<td class="org-left">United methylation in windows</td>
</tr>

<tr>
<td class="org-left"><code>emseq/dmr/annotation/{experiment}_annotated.tsv</code></td>
<td class="org-left">Differentially methylated regions with genomic annotation</td>
</tr>

<tr>
<td class="org-left"><code>emseq/haplotype/{lib}.{ref}.{aligner}_mhl.txt</code></td>
<td class="org-left">Haplotype metrics per methylation haplotype block</td>
</tr>

<tr>
<td class="org-left"><code>emseq/deconv/uxm_results.csv</code></td>
<td class="org-left">Tissue-of-origin fractions for all libraries</td>
</tr>
</tbody>
</table>

*Table 1: Output files of the core processing and the downstream analysis.*


<a id="design-decisions"></a>

# Design decisions


<a id="read-end-handling"></a>

## Read-end handling

Methylation is counted with the first five sequenced bases of read 2 excluded. End repair fills recessed 3' ends with unmethylated cytosine before conversion, so the last bases of each original strand read as unmethylated whatever the template carried. Read 2 starts on that end, and read 1 reaches it on fragments no longer than the read. In one plasma cfDNA library carrying both spike-in controls (463,461 usable pUC19 reads), pUC19 CpG methylation was 18.1% at the base that begins read 2 and 57.8% at the last base of read 1, compared to 96.8% in the middle of the read. At the two 5' ends it was 99.2% and 97.5% ([Figure 2](#fig-end-handling-controls)). Lambda stayed at or below 2.1% at all four read ends.  

The loss is therefore at the 3' end of the original strand and absent from the 5' end. The clip of 10 bases from both ends of both mates that [Bismark](https://github.com/FelixKrueger/Bismark/blob/master/docs/src/content/docs/usage/library-types.md) recommends for EM-seq would remove correctly called bases at the 5' ends. The [NEB EM-seq pipeline](https://github.com/nebiolabs/EM-seq) instead excludes bases at the start of read 2 when methylation is counted (`--nOT 0,0,0,5 --nOB 0,0,5,0`). On read pairs constructed by hand, those options excluded 4 or 5 bases of read 2 and left the same bases counted from read 1 where the two mates overlapped ([tests/unit/test\_methyldackel\_bounds.py](tests/unit/test_methyldackel_bounds.py)). On fragments no longer than the read, which are common in plasma, an exclusion that names read 2 alone leaves the filled-in bases in the counts.  

<a id="fig-end-handling-controls"></a>
![**The methylated control loses methylation at the 3' end of the original strand.** CpG methylation of pUC19 and lambda by position from the read end, with panels grouped by the strand end that the read end covers. n = 1 plasma cfDNA library. Bands are 95% Wilson intervals.](resources/figures/validation_end_handling_plasma-spiked.png)
*Figure 2: **The methylated control loses methylation at the 3' end of the original strand.** CpG methylation of pUC19 and lambda by position from the read end, with panels grouped by the strand end that the read end covers. n = 1 plasma cfDNA library. Bands are 95% Wilson intervals.*  

The loss was then measured on plasma reads themselves. Each CpG call was placed at its distance from the two ends of its original strand, on fragments of 160 bp or longer, where a base within 80 of one end is nearer that end than the other. The 5' end of a strand is not filled in, and over both strands the 5' and 3' tallies cover the same genomic positions, so the difference between the two ends at equal distance is the loss ([Figure 3](#fig-end-handling-strand-ends)). In 7 plasma libraries, methylation at the 3' end was 3.1% at the last base, 53.5% at base 10, and 71.6% at base 30, compared to 83.6% and 81.5% at bases 10 and 30 from the 5' end. The difference was 30.1 percentage points at base 10 (range 26.4 to 32.5 across libraries), 9.9 at base 30, and 5.7 at base 50, and it exceeded 5 points over 52 bases from the 3' end. On the pUC19 control of the library that carries it, the same difference exceeded 5 points over 12 bases. The loss in plasma is thus several times longer than on the fragmented control, consistent with the single-stranded overhangs of plasma fragments.  

Excluding read ends at calling removes part of the affected calls. The NEB options left 93.1% of the calls within the 52 bases nearest the 3' end in the counts, and a mask of 10 bases at both ends of both mates left 81.6%. Global CpG methylation over all calls was 74.5% with all positions counted, 75.5% with the NEB options, and 76.0% with the 10-base mask. The rise over all positions counted was 0.96 points (95% interval 0.91 to 1.01) for the NEB options and 1.43 points (1.34 to 1.53) for the 10-base mask, paired over the 7 libraries, with 2.0% and 12.0% of calls removed. Both masks therefore leave most of the loss in the counts.  

At equal distance from the 3' end, calls made by read 1 were 7.4 points higher than calls made by read 2 in plasma, and 1.9 points higher on pUC19. Read 1 reaches the 3' end only on fragments no longer than the read, so in plasma the two mates sample fragments of different lengths, and this comparison cannot tell an effect of the sequencing cycle from an effect of fragment length.  

<a id="fig-end-handling-strand-ends"></a>
![**Plasma reads lose methylation over tens of bases at the 3' end of the original strand.** CpG methylation by distance from the 5' and 3' ends of the original strand, on fragments of 160 bp or longer, counted once per fragment. Thick lines pool the libraries and thin lines are single libraries. Plasma, n = 7 libraries. pUC19, n = 1 library. Bands are 95% Wilson intervals.](resources/figures/validation_fragment_end_plasma-pilot.png)
*Figure 3: **Plasma reads lose methylation over tens of bases at the 3' end of the original strand.** CpG methylation by distance from the 5' and 3' ends of the original strand, on fragments of 160 bp or longer, counted once per fragment. Thick lines pool the libraries and thin lines are single libraries. Plasma, n = 7 libraries. pUC19, n = 1 library. Bands are 95% Wilson intervals.*  

The measures were repeated on public EM-seq data from another laboratory, four NEB libraries of BioProject PRJNA591788 (<a href="#citeproc_bib_item_12">Vaisvila et al. 2021</a>): two of NA12878 genomic DNA sheared before library preparation and two of plasma cfDNA, all carrying lambda and pUC19. The controls were tallied on all read pairs and the human reads on a seeded subsample of 30 million read pairs per library. On pUC19, CpG methylation at the base that begins read 2 was 1.6% in the sheared-DNA libraries and 2.6% in the cfDNA libraries, compared to 85.8% and 97.4% in the middle of the read. On human reads the length of the 3' loss depended on the source of the DNA ([Figure 4](#fig-end-handling-public-strand-ends)). In the sheared DNA, the difference between the 5' and 3' ends was 3.6 points at base 10 and exceeded 5 points over 5 bases, and the NEB options left 26.3% of the calls in those bases in the counts. In the cfDNA, the difference was 31.0 points at base 10, 13.0 at base 30, and 7.2 at base 80, the end of the measured range, and the NEB options left 94.5% of the calls within 80 bases of the 3' end in the counts. Part of the difference far from the end came from the 5' side: in the cfDNA, methylation 10 bases from the 5' end was 80.0%, compared to 72.1% over all calls. The NEB options raised global CpG methylation by 1.13 points in the cfDNA and 0.75 points in the sheared DNA. The long 3' loss of the pilot plasma libraries thus recurred in cfDNA prepared in another laboratory and was limited to a few bases in sheared genomic DNA prepared with the same kit.  

<a id="fig-end-handling-public-strand-ends"></a>
![**In public EM-seq data, the 3' loss extends over tens of bases in cfDNA and over a few bases in sheared genomic DNA.** CpG methylation by distance from the 5' and 3' ends of the original strand, on fragments of 160 bp or longer, from a seeded subsample of 30 million read pairs per library. NEB libraries of BioProject PRJNA591788: cfDNA, n = 2 libraries; sheared NA12878 DNA, n = 2 libraries. Each line pools the two libraries of one kind.](resources/figures/validation_fragment_end_public-cfdna-host.png)
*Figure 4: **In public EM-seq data, the 3' loss extends over tens of bases in cfDNA and over a few bases in sheared genomic DNA.** CpG methylation by distance from the 5' and 3' ends of the original strand, on fragments of 160 bp or longer, from a seeded subsample of 30 million read pairs per library. NEB libraries of BioProject PRJNA591788: cfDNA, n = 2 libraries; sheared NA12878 DNA, n = 2 libraries. Each line pools the two libraries of one kind.*  

The two masks were also compared on a cohort large enough to show whether the change they make to the counts reaches a comparison between groups. We called 55 plasma cfDNA libraries of a breast cancer study with all aligned reads, 40 from patients in a second-line breast cancer cohort and 15 from healthy donors, once with all positions counted, once with the NEB options, and once with the 10-base mask. Global CpG methylation over autosomal calls rose by 0.95 points (95% interval 0.92 to 0.98) with the NEB options and by 1.55 points (1.50 to 1.60) with the 10-base mask (paired t test over libraries). The rise in patients minus the rise in healthy donors was -0.09 points (-0.14 to -0.03) and -0.14 points (-0.25 to -0.04) by Welch's t test. The patient minus healthy difference in global methylation was therefore -3.66 points (-4.88 to -2.45) with all positions counted, -3.75 with the NEB options, and -3.81 with the 10-base mask.  

Differentially methylated 10 kb tiles were then called between the two groups with methylKit (<a href="#citeproc_bib_item_1">Akalin et al. 2012</a>), with the breast study's coverage of at least 10 and at least 3 libraries per group, on 282,078 tiles. The methylation difference of each tile correlated with its value under all positions counted at r = 0.998 for the NEB options and r = 0.988 for the 10-base mask, and no tile called under two settings had differences of opposite sign. The masks called more tiles, however. At a q-value below 0.01 and a difference of at least 25 points, 76 tiles were called with all positions counted, 80 with the NEB options, and 94 with the 10-base mask, and the calls shared with all positions counted were 90% and 53% of the union of the two call sets. At a q-value below 0.05 and a difference of at least 10 points, the counts were 1,401, 1,535, and 1,640, and the shares 90% and 77%. Excluding read ends thus moved the group difference in global methylation by less than a fifth of a point, but it changed which tiles a study of this size reports, by about half of the strict call set for the 10-base mask. The cohort comes from one study, the 3' loss itself was not measured in it, and both masks leave most of the loss in the counts, so the effect of the full loss on a group comparison is not known from it.  

The pipeline excludes the NEB bases, which lie in the filled-in 3' end, and does not use the 10-base mask, which also removes correctly called bases at the 5' ends. The rest of the 3' loss stays in the counts. JEEPERS, which corrects the calls at filled-in 3' ends from other fragments that cover the same CpGs, was tested on chromosome 22 of one deep plasma library. Within 10 bases of the 3' end, the share of calls that agreed with the call on the partner strand of the same fragment rose from 61.5% to 84.4% (difference 22.9 points, 95% interval 18.8 to 26.9), with every call kept ([results/jeepers-test-plasma-spiked/second-run/opposite\_strand.tsv](results/jeepers-test-plasma-spiked/second-run/opposite_strand.tsv)). Its run time and memory on a whole genome have not been measured.  


<a id="pairs-without-the-proper-pair-flag"></a>

## Pairs without the proper-pair flag

BWA-MEM sets the proper-pair flag from a range of fragment lengths that it infers for each batch of reads. In a plasma library that range is narrow, because most fragments are close to 167 bp, and pairs from longer fragments are left unflagged although they are correctly placed. A filter on the flag then removes the fragments of two nucleosomes and longer. The pipeline tests each pair directly, for the same chromosome, inward orientation, and a fragment of at most 1,000 bp, in both the duplicate-marking rule and the filtering rule ([tests/unit/test\_pair\_filter.py](tests/unit/test_pair_filter.py)).  

MethylDackel applies the flag a second time. By default it ignores a pair whose flag is unset, which would discard at calling the pairs that the pair test keeps. The calling commands therefore pass `--keepDiscordant` ([tests/unit/test\_methyldackel\_unflagged\_pairs.py](tests/unit/test_methyldackel_unflagged_pairs.py)).  


## Duplicates flagged before removal

Duplicate marking and duplicate removal are separate rules, so that the alignment with duplicates flagged but present remains on disk beside the filtered one. The flagged alignment holds the pairs of all mapping qualities and regions, and it is the input for methylation by read position and for any analysis that needs the pairs the filter removes.  


# Limitations


## End-repair fill-in is uncorrected

The methylation loss at the 3' end of each original strand enters the per-CpG counts. In plasma it exceeded 5 percentage points over 52 bases from the 3' end ([Read-end handling](#read-end-handling)), which is a large part of a fragment near 167 bp, and neither read-end mask tested removes most of it. The loss has been reported to be larger in plasma from patients with cancer than in plasma from controls (<a href="#citeproc_bib_item_13">Wang et al. 2026</a>), so it can differ between the groups a study compares. The plasma measurement rests on 7 libraries from one sequencing run at about 2 million read pairs each, on fragments of 160 bp or longer, and the control measurement on one library.  


## Spike-in reads in the sample alignment

The sample reference holds no lambda or pUC19 sequence, so a control read can align to a human locus and be counted there. The share of control reads that reach the filtered sample alignment has not been measured for this pipeline.  


## Two groups per experiment

Differential methylation compares two groups. A design with more groups is run as pairwise experiments, and a continuous outcome has to be made binary or analyzed outside the pipeline.  


# Prior work

The pipeline assembles published tools, and its processing follows the NEB pipeline most closely. [Table 2](#tab-prior-work) sets it beside the pipelines in common use for EM-seq.  

<table id="tab-prior-work" border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">Pipeline</th>
<th scope="col" class="org-left">Trimming</th>
<th scope="col" class="org-left">Aligner and caller</th>
<th scope="col" class="org-left">Read ends at calling</th>
</tr>
</thead>
<tbody>
<tr>
<td class="org-left">This pipeline</td>
<td class="org-left">fastp, adapters only</td>
<td class="org-left">bwa-meth or BISCUIT, MethylDackel</td>
<td class="org-left">All positions counted</td>
</tr>

<tr>
<td class="org-left"><a href="https://github.com/nebiolabs/EM-seq">NEB EM-seq</a></td>
<td class="org-left">fastp, adapters only</td>
<td class="org-left">bwa-meth, MethylDackel</td>
<td class="org-left">Bases at the start of read 2 excluded</td>
</tr>

<tr>
<td class="org-left"><a href="https://github.com/nf-core/methylseq">nf-core/methylseq</a> (<code>--em_seq</code>)</td>
<td class="org-left">Trim Galore, 10 bases clipped at both ends of both mates</td>
<td class="org-left">Bismark, or bwa-meth with MethylDackel</td>
<td class="org-left">No exclusion with MethylDackel</td>
</tr>

<tr>
<td class="org-left"><a href="https://github.com/maxplanck-ie/snakepipes">snakePipes</a> WGBS</td>
<td class="org-left">fastp</td>
<td class="org-left">bwa-meth, MethylDackel</td>
<td class="org-left">MethylDackel's suggested bounds applied</td>
</tr>
</tbody>
</table>

*Table 2: Read-end handling, alignment, and calling in pipelines used for EM-seq.*

The pipeline differs from the NEB pipeline at read ends, where it counts all positions. The tools it calls are bwa-meth, BISCUIT (<a href="#citeproc_bib_item_14">Zhou et al. 2024</a>), MethylDackel, dupsifter, methylKit (<a href="#citeproc_bib_item_1">Akalin et al. 2012</a>), mosdepth (<a href="#citeproc_bib_item_11">Pedersen and Quinlan 2018</a>), and wgbstools (<a href="#citeproc_bib_item_8">Loyfer, Rosenski, and Kaplan 2026</a>).  


# Abbreviations

-   **BAM:** binary alignment map, the alignment file format
-   **cfDNA:** cell-free DNA
-   **EM-seq:** enzymatic methyl sequencing
-   **MAPQ:** mapping quality
-   **MHB:** methylation haplotype block
-   **NEB:** New England Biolabs
-   **TET2:** ten-eleven translocation methylcytosine dioxygenase 2
-   **TLEN:** the fragment length field of an alignment record
-   **UXM:** the deconvolution tool that estimates tissue-of-origin fractions from read-level methylation


# Glossary

-   **Conversion:** The enzymatic change of unmethylated cytosine to uracil, read as thymine. Under-conversion leaves an unmethylated cytosine reading as methylated.
-   **Fill-in:** The bases a polymerase adds at a recessed 3' end during end repair. They are made with unmethylated cytosine and read as unmethylated.
-   **M-bias:** Methylation tabulated by position along the read. A level that changes toward a read end marks a position-dependent artifact.
-   **Original top and original bottom strand:** The two strands of the input fragment. After conversion they are no longer complementary, and a pair is assigned to one of them by the orientation of read 1.
-   **Proper-pair flag:** The alignment flag an aligner sets when the two mates of a pair lie at an expected distance and orientation.
-   **Spike-in control:** DNA of known methylation added to a sample before library preparation. Lambda is unmethylated and pUC19 is methylated at CpG sites.
-   **Usable read:** A read that MethylDackel counts by default, with a mapping quality of 10 or more and without the secondary, supplementary, duplicate, or failed-quality-control flags.


# References

<div class="csl-bib-body">
  <div class="csl-entry"><a id="citeproc_bib_item_1"></a>Akalin, Altuna, Matthias Kormaksson, Sheng Li, Francine E Garrett-Bakelman, Maria E Figueroa, Ari Melnick, and Christopher E Mason. 2012. “methylKit: A Comprehensive R Package for the Analysis of Genome-Wide DNA Methylation Profiles.” <i>Genome Biology</i> 13 (10): R87. <a href="https://doi.org/10.1186/gb-2012-13-10-r87">https://doi.org/10.1186/gb-2012-13-10-r87</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_2"></a>Chauhan, Pradeep S., Irfan Alahi, Savar Sinha, Elisa M. Ledet, Ryan Mueller, Jessica Linford, Alexander L. Shiang, et al. 2024. “Genomic and Epigenomic Analysis of Plasma Cell-Free Dna Identifies Stemness Features Associated with Worse Survival in Lethal Prostate Cancer.” <i>Clinical Cancer Research</i>, August. <a href="https://doi.org/10.1158/1078-0432.ccr-24-1658">https://doi.org/10.1158/1078-0432.ccr-24-1658</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_3"></a>Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. 2018. “Fastp: An Ultra-Fast All-in-One FASTQ Preprocessor.” <i>Bioinformatics</i> 34 (17): i884–90. <a href="https://doi.org/10.1093/bioinformatics/bty560">https://doi.org/10.1093/bioinformatics/bty560</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_4"></a>Guo, Shicheng, Dinh Diep, Nongluk Plongthongkum, Ho-Lim Fung, Kang Zhang, and Kun Zhang. 2017. “Identification of Methylation Haplotype Blocks Aids in Deconvolution of Heterogeneous Tissue Samples and Tumor Tissue-of-Origin Mapping from Plasma DNA.” <i>Nature Genetics</i> 49 (4): 635–42. <a href="https://doi.org/10.1038/ng.3805">https://doi.org/10.1038/ng.3805</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_5"></a>Hansen, Kasper D, Benjamin Langmead, and Rafael A Irizarry. 2012. “BSmooth: From Whole Genome Bisulfite Sequencing Reads to Differentially Methylated Regions.” <i>Genome Biology</i> 13 (10): R83. <a href="https://doi.org/10.1186/gb-2012-13-10-r83">https://doi.org/10.1186/gb-2012-13-10-r83</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_6"></a>Li, Heng. 2011. “A Statistical Framework for SNP Calling, Mutation Discovery, Association Mapping and Population Genetical Parameter Estimation from Sequencing Data.” <i>Bioinformatics</i> 27 (21): 2987–93. <a href="https://doi.org/10.1093/bioinformatics/btr509">https://doi.org/10.1093/bioinformatics/btr509</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_7"></a>Loyfer, Netanel, Judith Magenheim, Ayelet Peretz, Gordon Cann, Joerg Bredno, Agnes Klochendler, Ilana Fox-Fisher, et al. 2023. “A Dna Methylation Atlas of Normal Human Cell Types.” <i>Nature</i> 613 (7943): 355–64. <a href="https://doi.org/10.1038/s41586-022-05580-6">https://doi.org/10.1038/s41586-022-05580-6</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_8"></a>Loyfer, Netanel, Jonathan Rosenski, and Tommy Kaplan. 2026. “Wgbstools : a Computational Suite for Dna Methylation Sequencing Data Analysis.” <i>Life Science Alliance</i> 9 (4): e202503514. <a href="https://doi.org/10.26508/lsa.202503514">https://doi.org/10.26508/lsa.202503514</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_9"></a>Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, Michael B. Hall, Christopher H. Tomkins-Tinch, Vanessa Sochat, Jan Forster, et al. 2021. “Sustainable Data Analysis with Snakemake.” <i>F1000research</i> 10 (April): 33. <a href="https://doi.org/10.12688/f1000research.29032.2">https://doi.org/10.12688/f1000research.29032.2</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_10"></a>New England Biolabs. 2023. “Nebnext Enzymatic Methyl-Seq Kit E7120 Manual.” <a href="https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7120.pdf?rev=8572f755e2964742bb5af7532adef458&hash=59A7A19B73DD2E1B07B17A9F34901BC6">https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7120.pdf?rev=8572f755e2964742bb5af7532adef458&#38;hash=59A7A19B73DD2E1B07B17A9F34901BC6</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_11"></a>Pedersen, Brent S., and Aaron R. Quinlan. 2018. “Mosdepth: Quick Coverage Calculation for Genomes and Exomes.” <i>Bioinformatics</i> 34 (5): 867–68. <a href="https://doi.org/10.1093/bioinformatics/btx699">https://doi.org/10.1093/bioinformatics/btx699</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_12"></a>Vaisvila, Romualdas, V.K. Chaithanya Ponnaluri, Zhiyi Sun, Bradley W. Langhorst, Lana Saleh, Shengxi Guan, Nan Dai, et al. 2021. “Enzymatic Methyl Sequencing Detects DNA Methylation at Single-Base Resolution from Picograms of DNA.” <i>Genome Research</i> 31 (7): 1280–89. <a href="https://doi.org/10.1101/gr.266551.120">https://doi.org/10.1101/gr.266551.120</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_13"></a>Wang, Rui, Emily G. Hamilton, Diego Almanza, Rene Bonilla, Takeshi Sugio, Angela Hui, Feng Tian, et al. 2026. “Enhanced Cfdna Methylation Profiling by Detection and Correction of Jagged Ends.” <i>Nature Methods</i>, October. <a href="https://doi.org/10.1038/s41592-026-03278-5">https://doi.org/10.1038/s41592-026-03278-5</a>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_14"></a>Zhou, Wanding, Benjamin K Johnson, Jacob Morrison, Ian Beddows, James Eapen, Efrat Katsman, Ayush Semwal, et al. 2024. “Biscuit: An Efficient, Standards-Compliant Tool Suite for Simultaneous Genetic and Epigenetic Inference in Bulk and Single-Cell Studies.” <i>Nucleic Acids Research</i> 52 (6): e32. <a href="https://doi.org/10.1093/nar/gkae097">https://doi.org/10.1093/nar/gkae097</a>.</div>
</div>

