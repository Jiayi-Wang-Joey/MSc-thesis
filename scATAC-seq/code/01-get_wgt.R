suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg38)
    source("~/chromVAR/R/getCounts.R")
    library(data.table)
    library(SummarizedExperiment)
    library(fields)
})

m <- wcs$mot
atacFrag <- readRDS(args$frg)
smt <- strsplit(wcs$smt, "_")[[1]]
mdr <- ifelse(wcs$mdr=="unmoderated", FALSE, TRUE)
pf <- list.files(paste0("/mnt/plger/plger/DTFAB/fullFrags/",m,"/peaks"),
    full.names = TRUE)
peaks <- import.bed(con=pf)
atacFrag <- lapply(atacFrag, data.table)
if (smt[1]=="none") aRange <- 0 else aRange <- as.numeric(smt[2])

se <- getCounts(atacFrag = atacFrag,
    ranges = peaks,
    mode = "weight",
    genome = BSgenome.Hsapiens.UCSC.hg38,
    species = "Homo sapiens",
    width = 300,
    smooth = smt[1],
    nGCBins = 10,
    nWidthBins = 30,
    aRange = aRange,
    moderating = mdr, 
    peakWeight = wcs$pkw,
    singleCell = TRUE)

saveRDS(se, args$res)