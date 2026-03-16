#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/emseq/emseq.org
# Author:  Jeff Szymanski
# Tangled: 2026-03-16 13:58:24
# ============================================================

# emseq_annotate_methylkit.R — Annotate methylKit DB with CpG context, gene parts, TSS
# Tangled from emseq.org; do not edit directly.

suppressPackageStartupMessages({
  library(argparse)
  library(methylKit)
  library(GenomicRanges); library(GenomeInfoDb)
  library(dplyr); library(tidyr); library(readr)
  library(annotatr)
  suppressWarnings(library(genomation))
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

parser <- ArgumentParser(description="Annotate methylKit DB with CpG context, gene parts, TSS, and IDs")
parser$add_argument("--db",  required=TRUE, help="Path to methylKit tabix DB (.txt.bgz)")
parser$add_argument("--out", required=TRUE, help="Output TSV")
parser$add_argument("--prom-up", type="integer", default=2000, help="Promoter upstream bp [2000]")
parser$add_argument("--prom-dn", type="integer", default=500,  help="Promoter downstream bp [500]")
args <- parser$parse_args()

muffle_outofbound <- function(expr) {
  withCallingHandlers(expr, warning=function(w){
    if (grepl("out-of-bound ranges located on sequences", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}

md <- methylKit:::readMethylDB(args$db)
gr <- as(md, "GRanges"); mcols(gr)$qid <- seq_along(gr)

ann_cpg <- suppressMessages(build_annotations(genome="hg38", annotations="hg38_cpgs"))
a <- annotate_regions(gr, ann_cpg, ignore.strand=TRUE, quiet=TRUE) |>
  as_tibble() |>
  transmute(q_chr=as.character(seqnames), q_start=start, q_end=end, type=annot.type)
keys <- tibble(q_chr=as.character(seqnames(gr)), q_start=start(gr), q_end=end(gr), qid=mcols(gr)$qid)
cpg_flags <- a |>
  right_join(keys, by=c("q_chr","q_start","q_end")) |>
  group_by(qid) |>
  summarise(
    is_island  = any(type=="hg38_cpg_islands",  na.rm=TRUE),
    is_shore   = any(type=="hg38_cpg_shores",   na.rm=TRUE),
    is_shelf   = any(type=="hg38_cpg_shelves",  na.rm=TRUE),
    is_opensea = any(type=="hg38_cpg_inter",    na.rm=TRUE),
    .groups="drop"
  )

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
ex    <- exons(txdb)
intr  <- unlist(intronsByTranscript(txdb), use.names=FALSE)
utr5  <- unlist(fiveUTRsByTranscript(txdb),  use.names=FALSE)
utr3  <- unlist(threeUTRsByTranscript(txdb), use.names=FALSE)
prom  <- promoters(txdb, upstream=args$prom_up, downstream=args$prom_dn)
tsspt <- promoters(txdb, upstream=0, downstream=1)

for (x in list(ex,intr,utr5,utr3,prom,tsspt)) seqlevelsStyle(x) <- seqlevelsStyle(gr)
muffle_outofbound({
  ex    <- keepStandardChromosomes(ex,   pruning.mode="coarse")
  intr  <- keepStandardChromosomes(intr, pruning.mode="coarse")
  utr5  <- keepStandardChromosomes(utr5, pruning.mode="coarse")
  utr3  <- keepStandardChromosomes(utr3, pruning.mode="coarse")
  prom  <- keepStandardChromosomes(prom, pruning.mode="coarse")
  tsspt <- keepStandardChromosomes(tsspt,pruning.mode="coarse")
  wanted <- intersect(seqlevels(gr), seqlevels(ex))
  ex    <- keepSeqlevels(ex,    wanted, pruning.mode="coarse")
  intr  <- keepSeqlevels(intr,  wanted, pruning.mode="coarse")
  utr5  <- keepSeqlevels(utr5,  wanted, pruning.mode="coarse")
  utr3  <- keepSeqlevels(utr3,  wanted, pruning.mode="coarse")
  prom  <- keepSeqlevels(prom,  wanted, pruning.mode="coarse")
  tsspt <- keepSeqlevels(tsspt, wanted, pruning.mode="coarse")
})

is_prom     <- countOverlaps(gr, prom) > 0L
is_5utr     <- countOverlaps(gr, utr5) > 0L
is_3utr     <- countOverlaps(gr, utr3) > 0L
is_exon_raw <- countOverlaps(gr, ex)   > 0L
is_exon     <- is_exon_raw & !is_5utr & !is_3utr
is_intron   <- countOverlaps(gr, intr) > 0L

gene_part_primary <- ifelse(is_prom, "promoter",
                        ifelse(is_5utr, "5UTR",
                        ifelse(is_3utr, "3UTR",
                        ifelse(is_exon, "exon",
                        ifelse(is_intron, "intron", "intergenic")))))

nn <- distanceToNearest(gr, tsspt, ignore.strand=FALSE)
dist_tbl <- tibble(
  qid         = mcols(gr)$qid[queryHits(nn)],
  dist_to_TSS = mcols(nn)$distance,
  tss_idx     = subjectHits(nn)
)

tx_for_tss <- mcols(tsspt)$tx_id[unique(dist_tbl$tss_idx)]
tx_map <- suppressMessages(AnnotationDbi::select(
  x=txdb, keys=as.character(tx_for_tss),
  keytype="TXID", columns=c("TXID","GENEID")
))
tss2tx <- tibble(
  tss_idx = unique(dist_tbl$tss_idx),
  TXID    = as.character(mcols(tsspt)$tx_id[unique(dist_tbl$tss_idx)])
) |>
  left_join(tx_map, by="TXID") |>
  distinct(tss_idx, .keep_all=TRUE) |>
  rename(ENTREZID = GENEID)

gene_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(na.omit(tss2tx$ENTREZID)),
  keytype = "ENTREZID",
  columns = c("SYMBOL","ENSEMBL")
)) |>
  distinct(ENTREZID, .keep_all=TRUE)

tss2tx <- tss2tx |> left_join(gene_map, by="ENTREZID")

dist_tbl <- dist_tbl |>
  left_join(tss2tx |> dplyr::select(tss_idx, ENTREZID, SYMBOL, ENSEMBL), by="tss_idx") |>
  dplyr::select(-tss_idx)

gene_flags <- tibble(
  qid = mcols(gr)$qid,
  is_promoter = is_prom, is_5utr = is_5utr, is_3utr = is_3utr,
  is_exon = is_exon, is_intron = is_intron,
  gene_part_primary = gene_part_primary
) |>
  left_join(dist_tbl, by="qid")

out <- as_tibble(getData(md)) |>
  mutate(qid = dplyr::row_number()) |>
  left_join(cpg_flags,  by="qid") |>
  left_join(gene_flags, by="qid") |>
  dplyr::select(
    chr, start, end, strand, qid,
    dplyr::starts_with("is_"),
    gene_part_primary, dist_to_TSS,
    ENTREZID, SYMBOL, ENSEMBL,
    dplyr::everything()
  )

readr::write_tsv(out, args$out)
