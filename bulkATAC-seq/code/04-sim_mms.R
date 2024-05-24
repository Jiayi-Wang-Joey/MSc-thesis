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
    library(rtracklayer)
})
setDTthreads(8)
m <- wcs$mtf
e <- wcs$eft
nrm <- ifelse(wcs$nrm=="LBS", TRUE, FALSE)
sym <- ifelse(wcs$sym=="sym", TRUE, FALSE)
#uni <- ifelse(wcs$uni=="uniform", TRUE, FALSE)
se <- readRDS(paste0("/mnt/plger/jwang/data/sim/01-total/total-",m,"-",e,".rds"))
atacFrag <- readRDS(paste0("/mnt/plger/jwang/data/sim/00-frg/",m,"-",e,".rds"))
folder <- paste0(m, "_haploinsufficiency_", e, "_FALSE/peaks")
if (e !="0") {
    pf <- list.files(paste0("/mnt/plger/jwang/sim_data_es/",folder), 
        "^merged*",
        full.names = TRUE)
} else {
    pf <- "/mnt/plger/esonder/R/tf_activity_benchmark/DTFAB/simulations/data/peaks/merged_summits.bed"
}
peaks <- import.bed(pf)
motifRanges <- readRDS(paste0("/mnt/plger/jwang/data/sim/motifs/",m,"-",e,".rds"))

res <- computeMotifActivityScore(se=se,
    atacFrag=atacFrag,
    motifRanges=motifRanges,
    species="Homo sapiens",
    genome=BSgenome.Hsapiens.UCSC.hg38,
    symmetric=sym,
    libNorm=nrm,
    niterations=2000)

saveRDS(list(wcs=wcs,res=res), args$res)