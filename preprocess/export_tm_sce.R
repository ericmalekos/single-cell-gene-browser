#!/usr/bin/env Rscript
# Export a SingleCellExperiment .rds to flat files with system R + Matrix only (no Bioconductor).
#
#   Rscript export_tm_sce.R --rds sce.rds --outdir DIR [--no-matrix] [--genes Gapdh,Zeb2]
#
# Writes to DIR:
#   genes.txt            rownames (gene symbols)
#   cells.tsv            cell_id + every colData column
#   umap.tsv             cell_id, umap_1, umap_2 (reducedDim UMAP_refined, else UMAP)
#   dims.txt             n_genes, n_cells, nnz of logcounts
#   x.f32 i.i32 p.i32    logcounts dgCMatrix slots (genes x cells, compressed by cell), little-endian
#   check_genes.tsv      dense logcounts of --genes, one column per gene
#
# Slots are reached through attributes(): `$` or slotNames() on the S4 objects would try to load
# SingleCellExperiment / GenomicRanges, which are not installed.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
rds <- get_arg("--rds")
outdir <- get_arg("--outdir")
with_matrix <- !("--no-matrix" %in% args)
check_genes <- get_arg("--genes", "")
check_genes <- if (nzchar(check_genes)) strsplit(check_genes, ",")[[1]] else character(0)
if (is.null(rds) || is.null(outdir)) {
  stop("usage: export_tm_sce.R --rds FILE --outdir DIR [--no-matrix] [--genes a,b]")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("reading", rds, "\n")
sce <- readRDS(rds)
a <- attributes(sce)
rm(sce)

assays <- attributes(attributes(a$assays)$data)$listData
cat("assays:", paste(names(assays), collapse = ", "), "\n")
lc <- assays$logcounts
if (is.null(lc)) stop("no logcounts assay")
# drop the counts assay before writing anything large
rm(assays)
a$assays <- NULL
invisible(gc())

dims <- attr(lc, "Dim")
dn <- attr(lc, "Dimnames")
rr_names <- attributes(attributes(a$rowRanges)$partitioning)$NAMES
genes <- if (is.null(dn[[1]])) rr_names else dn[[1]]
if (!is.null(rr_names) && !identical(genes, rr_names)) stop("logcounts rownames differ from rowRanges names")
stopifnot(length(genes) == dims[1])

cd <- attributes(a$colData)
cells <- if (is.null(cd$rownames)) dn[[2]] else cd$rownames
stopifnot(length(cells) == dims[2])
cols <- lapply(cd$listData, function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.character(v)) v <- gsub("[\t\r\n]", " ", v)
  v
})
coldata <- data.frame(cell_id = cells, cols, check.names = FALSE, stringsAsFactors = FALSE)
write.table(coldata, file.path(outdir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rdims <- attributes(attributes(a$int_colData)$listData$reducedDims)$listData
um_name <- if ("UMAP_refined" %in% names(rdims)) "UMAP_refined" else "UMAP"
um <- rdims[[um_name]]
if (is.null(um) || nrow(um) != dims[2]) stop("no usable UMAP among reducedDims: ", paste(names(rdims), collapse = ", "))
cat("reducedDim:", um_name, "\n")
write.table(data.frame(cell_id = cells, umap_1 = um[, 1], umap_2 = um[, 2]),
            file.path(outdir, "umap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(genes, file.path(outdir, "genes.txt"))
writeLines(as.character(c(dims, length(attr(lc, "x")))), file.path(outdir, "dims.txt"))

# writeBin writes at most 2^31 - 1 bytes per call, so long vectors go out in slices
write_slices <- function(v, path, size) {
  con <- file(path, "wb")
  on.exit(close(con))
  step <- 2^27
  s <- 1
  while (s <= length(v)) {
    e <- min(length(v), s + step - 1)
    writeBin(v[s:e], con, size = size, endian = "little")
    s <- e + 1
  }
}
if (with_matrix) {
  cat("writing logcounts slots: nnz =", length(attr(lc, "x")), "\n")
  write_slices(attr(lc, "x"), file.path(outdir, "x.f32"), 4L)   # double -> float32
  write_slices(attr(lc, "i"), file.path(outdir, "i.i32"), 4L)
  write_slices(attr(lc, "p"), file.path(outdir, "p.i32"), 4L)
}

if (length(check_genes)) {
  i <- attr(lc, "i")
  p <- attr(lc, "p")
  x <- attr(lc, "x")
  out <- data.frame(cell_id = cells)
  for (g in check_genes) {
    r <- match(g, genes)
    if (is.na(r)) next
    v <- numeric(dims[2])
    pos <- which(i == (r - 1L))          # 1-based positions of the gene's non-zeros
    v[findInterval(pos - 1L, p)] <- x[pos]
    out[[g]] <- v
  }
  write.table(out, file.path(outdir, "check_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("done\n")
