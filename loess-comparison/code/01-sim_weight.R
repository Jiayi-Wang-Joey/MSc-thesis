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
atacFrag <- readRDS(args$frg)
fd <- paste0(dat, wcs$mtf, "_haploinsufficiency_0_FALSE_", wcs$dsg,"/peaks")

pf <- list.files(fd, "^merged*", full.names = TRUE)

peaks <- import.bed(con=pf)

se <- getCounts(atacFrag = atacFrag,
  ranges = peaks,
  rowType = "peak",
  mode = "weight",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "Homo sapiens",
  width = 300,
  nGCBins = 10,
  nWidthBins = 30,peakWeight = wcs$pkw)

saveRDS(se, args$res)