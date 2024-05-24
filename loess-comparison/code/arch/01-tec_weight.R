suppressPackageStartupMessages({
  library(SummarizedExperiment)
  source("code/utils.R")
})

# only test span

se <- readRDS("/mnt/plger/jwang/Cusanovich/weights_HS.rds")
low <- c("GSM5290881", "GSM5290882", "GSM5290883", 
  "GSM5290884", "GSM5290885", "GSM5290886")

high <- c("GSM5290887", "GSM5290888", "GSM5290889",
  "GSM5290890", "GSM5290891", "GSM5290892")

se <- se[,c(low,high)]
counts <- assay(se, "counts")

newCounts <- .weightPeaks(counts, 
  span=as.numeric(wcs$span))

assay(se,"counts") <- newCounts

saveRDS(se, args$res)
