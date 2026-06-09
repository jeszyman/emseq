#!/usr/bin/env Rscript
## emseq_checks.R — data-integrity helpers for emseq analysis scripts.
## Source with: source("scripts/emseq_checks.R")

## Log row count to stderr; return x invisibly for piping.
log_n <- function(x, label = "rows") {
  message(sprintf("[check] %s: %d", label, nrow(x)))
  invisible(x)
}

## Stop unless x has exactly n rows.
check_n <- function(x, n, label = "rows") {
  if (nrow(x) != n) stop(sprintf("[check] %s: expected %d, got %d", label, n, nrow(x)))
  invisible(x)
}

## Stop if x has zero rows.
check_nonempty <- function(x, label = "object") {
  if (nrow(x) == 0L) stop(sprintf("[check] %s is empty", label))
  invisible(x)
}

## Stop unless all `cols` are present in x.
assert_cols <- function(x, cols, label = "object") {
  miss <- setdiff(cols, colnames(x))
  if (length(miss) > 0) stop(sprintf("[check] %s missing columns: %s", label, paste(miss, collapse = ", ")))
  invisible(x)
}
