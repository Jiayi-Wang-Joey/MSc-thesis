suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg38)
    source("~/chromVAR/R/getCounts.R")
    library(data.table)
    library(SummarizedExperiment)
})

m <- wcs$mot
atacFrag <- readRDS(args$frg)
pf <- list.files(paste0("/mnt/plger/plger/DTFAB/fullFrags/",m,"/peaks"),
    full.names = TRUE)
peaks <- import.bed(con=pf)

se <- getCounts(atacFrag = atacFrag,
    ranges = peaks,
    mode = "total",
    genome = BSgenome.Hsapiens.UCSC.hg38,
    species = "Homo sapiens",
    width = 300, singleCell = TRUE)

saveRDS(se, args$res)
