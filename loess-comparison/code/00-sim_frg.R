suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(SummarizedExperiment)
  source("~/chromVAR/R/getCounts.R")
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(rtracklayer)
  library(fields)
})

if (wcs$eft=="FLD") {
  dat <- "/mnt/plger/jwang/sim_data_fld/"
} else {
  dat <- "/mnt/plger/jwang/sim_data_gc/"
}
fd <- paste0(dat, wcs$mtf, "_haploinsufficiency_0_FALSE_", wcs$dsg,
  "/seq_files")

fs <- list.files(fd, "\\.bed$", full.names=TRUE)
lst <- as.list(fs)
names(lst) <- basename(fs)
x <- .importFragments(lst)
saveRDS(x, args$res)

#pf <- "/mnt/plger/esonder/R/tf_activity_benchmark/DTFAB/simulations/data/peaks/merged_summits.bed"
#peaks <- import.bed(con=pf)