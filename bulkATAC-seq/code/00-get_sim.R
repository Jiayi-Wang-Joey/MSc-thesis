suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  source("~/chromVAR/R/getCounts.R")
})

m <- wcs$mtf
e <- wcs$eft

fd <- paste0("/mnt/plger/jwang/sim_data_es/",
    m, "_haploinsufficiency_", e, "_FALSE/seq_files")
fs <- list.files(fd, "\\.bam$|\\.bed$", full.names=TRUE)
lst <- as.list(fs)
names(lst) <- basename(fs)
x <- .importFragments(lst)
saveRDS(x, args$res)
