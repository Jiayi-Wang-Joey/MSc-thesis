suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

mot <- wcs$mot
fragDt <- readRDS("/mnt/plger/jwang/scATAC-seq/fragDt.rds")
nSample <- length(unique(fragDt[fragDt$motif==mot,]$sample))
ctrl <- paste0("sgNT_", seq_len(nSample))
dt <- fragDt %>%
    filter(motif==mot | sample %in% ctrl)
frg <- split(dt, dt$barcode)

saveRDS(frg, args$res)