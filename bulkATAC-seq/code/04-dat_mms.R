suppressPackageStartupMessages({
    library(chromVAR)
    source("~/chromVAR/R/motifMatch.R")
    library(reshape2)
    library(rtracklayer)
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("code/utils.R")
    library(MotifDb)
    library(SummarizedExperiment)
    library(data.table)
    library(BiocParallel)
    library(Matrix)
})
setDTthreads(8)
m <- wcs$mot
nrm <- ifelse(wcs$nrm=="LBS", TRUE, FALSE)
sym <- ifelse(wcs$sym=="sym", TRUE, FALSE)
se <- readRDS(paste0("/mnt/plger/jwang/data/dat/01-total/total-",m,".rds"))
atacFrag <- readRDS(paste0("/mnt/plger/jwang/data/dat/00-frg/",m,".rds"))
motifRanges <- readRDS(paste0("/mnt/plger/plger/DTFAB/fullFrags/",m,"/runATAC_results/others/pmoi.rds"))
mice <- c("BANP", "NR1H3", "NR1H4")

if (m %in% mice) {
    res <- computeMotifActivityScore(se=se,
        atacFrag=atacFrag,
        motifRanges=motifRanges,
        species="Mus_musculus",
        genome=BSgenome.Mmusculus.UCSC.mm10,
        symmetric=sym, 
        libNorm=nrm,
        niterations=2000)
    
} else {
    res <- computeMotifActivityScore(se=se,
        atacFrag=atacFrag,
        motifRanges=motifRanges,
        species="Homo sapiens",
        genome=BSgenome.Hsapiens.UCSC.hg38,
        symmetric=sym,
        libNorm=nrm,
        niterations=2000)
}

saveRDS(list(wcs,res), args$res)