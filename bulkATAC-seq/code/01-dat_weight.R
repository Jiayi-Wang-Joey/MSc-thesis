suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(fields)
    library(SummarizedExperiment)
    source("~/chromVAR/R/getCounts.R")
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(rtracklayer)
    library(edgeR)
    library(limma)
    library(dplyr)
    library(rsample)
    library("profmem")
})

m <- wcs$mot
atacFrag <- readRDS(args$frg)
smt <- strsplit(wcs$smt, "_")[[1]]
pkw <- wcs$pkw
mdr <- ifelse(wcs$mdr=="unmoderated", FALSE, TRUE)
mice <- c("BANP", "NR1H3", "NR1H4")
pf <- list.files(paste0("/mnt/plger/plger/DTFAB/fullFrags/", m,"/peaks"),
    full.names = TRUE)
peaks <- import.bed(con=pf)
if (smt[1]=="none") aRange <- 0 else aRange <- as.numeric(smt[2])
if (m %in% mice) {
    species <- "Mus_musculus"
    genome <- BSgenome.Mmusculus.UCSC.mm10
} else {
    species <- "Homo sapiens"
    genome <- BSgenome.Hsapiens.UCSC.hg38
}

se <- getCounts(atacFrag = atacFrag,
        ranges = peaks,
        mode = "weight",
        genome = genome,
        species = species,
        width = 300,
        smooth = smt[1],
        nGCBins = 10,
        nWidthBins = 30,
        aRange = aRange,
        peakWeight = wcs$pkw,
        moderating = mdr)


saveRDS(se, args$res)
