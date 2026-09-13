#!/usr/bin/env Rscript
# emseq_annotate_methylkit.R — Annotate methylKit DB with CpG context, gene parts, TSS

suppressPackageStartupMessages({
  library(argparse)
  library(methylKit)
  library(GenomicRanges); library(GenomeInfoDb)
  library(dplyr); library(tidyr); library(readr)
  suppressWarnings(library(genomation))
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(Rsamtools)
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

md <- methylKit:::readMethylDiffDB(args$db)

chroms <- Rsamtools::seqnamesTabix(args$db)

# --- Local CpG island annotations (replaces annotatr::build_annotations) ---
local_cpg <- file.path(dirname(args$db), "hg38_cpgIslandExt.txt.gz")

if (!file.exists(local_cpg) || file.info(local_cpg)$size < 1000) {
    message("Downloading UCSC CpG Islands...")
    system(sprintf("wget -q -U 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)' -O %s https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz", local_cpg))
}

cpg_df <- suppressMessages(readr::read_tsv(local_cpg, col_names = c("chr", "start", "end"), col_types = "-cii-------"))
if (nrow(cpg_df) == 0) stop("CpG Island file is empty or corrupted! Please delete it and try again.")

cpg_isl <- GRanges(seqnames = cpg_df$chr, ranges = IRanges(start = cpg_df$start, end = cpg_df$end))
cpg_isl$type <- rep("hg38_cpg_islands", length(cpg_isl))

shores_raw <- cpg_isl + 2000
start(shores_raw) <- pmax(1, start(shores_raw))
cpg_shores <- setdiff(shores_raw, cpg_isl)
cpg_shores$type <- rep("hg38_cpg_shores", length(cpg_shores))

shelves_raw <- cpg_isl + 4000
start(shelves_raw) <- pmax(1, start(shelves_raw))
cpg_shelves <- setdiff(shelves_raw, shores_raw)
cpg_shelves$type <- rep("hg38_cpg_shelves", length(cpg_shelves))

ann_cpg <- c(cpg_isl, cpg_shores, cpg_shelves)

# --- Build gene-model features once ---
txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
ex    <- exons(txdb)
intr  <- unlist(intronsByTranscript(txdb), use.names=FALSE)
utr5  <- unlist(fiveUTRsByTranscript(txdb),  use.names=FALSE)
utr3  <- unlist(threeUTRsByTranscript(txdb), use.names=FALSE)
prom  <- promoters(txdb, upstream=args$prom_up, downstream=args$prom_dn)
tsspt <- promoters(txdb, upstream=0, downstream=1)

first_chunk <- TRUE

for (chrom in chroms) {
    message(sprintf("Annotating chromosome: %s", chrom))

    chrom_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = 1, end = .Machine$integer.max))

    diff_sub <- selectByOverlap(md, chrom_gr)

    if (is.null(diff_sub) || nrow(diff_sub) == 0) next

    gr <- as(diff_sub, "GRanges")
    mcols(gr)$qid <- seq_along(gr)

    # Harmonize seqlevels on first chromosome (idempotent)
    if (first_chunk) {
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
    }

    # --- CpG context via local overlap ---
    hits <- findOverlaps(gr, ann_cpg, ignore.strand=TRUE)

    overlap_tbl <- tibble(
        qid = mcols(gr)$qid[queryHits(hits)],
        type = ann_cpg$type[subjectHits(hits)]
    )

    keys <- tibble(qid = mcols(gr)$qid)

    cpg_flags <- overlap_tbl |>
      right_join(keys, by = "qid") |>
      group_by(qid) |>
      summarise(
        is_island  = any(type == "hg38_cpg_islands", na.rm = TRUE),
        is_shore   = any(type == "hg38_cpg_shores",  na.rm = TRUE),
        is_shelf   = any(type == "hg38_cpg_shelves", na.rm = TRUE),
        is_opensea = !any(type %in% c("hg38_cpg_islands", "hg38_cpg_shores", "hg38_cpg_shelves"), na.rm = TRUE),
        .groups = "drop"
      )

    # --- Gene-model annotation ---
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

    out <- as_tibble(getData(diff_sub)) |>
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

    write.table(out, file = args$out, append = !first_chunk,
                col.names = first_chunk, row.names = FALSE, sep = "\t", quote = FALSE)
    first_chunk <- FALSE
}
