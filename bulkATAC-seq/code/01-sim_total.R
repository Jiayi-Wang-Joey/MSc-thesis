suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(SummarizedExperiment)
    source("~/chromVAR/R/getCounts.R")
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(rtracklayer)
    source("code/utils.R")
})

m <- wcs$mtf
e <- wcs$eft
atacFrag <- readRDS(args$frg)
folder <- paste0(m, "_haploinsufficiency_", e, "_FALSE/peaks")
if (e !="0") {
  pf <- list.files(paste0("/mnt/plger/jwang/sim_data_es/",folder), 
    "^merged*",
    full.names = TRUE)
} else {
  pf <- "/mnt/plger/esonder/R/tf_activity_benchmark/DTFAB/simulations/data/peaks/merged_summits.bed"
}
peaks <- import.bed(con=pf)
motifRanges <- getpmoi(genome=BSgenome.Hsapiens.UCSC.hg38,
         peaks=peaks,
         spec="Hsapiens",
         seqStyle="UCSC",
         srcFolder="/mnt/plger/plger/DTFAB/Scripts")
saveRDS(motifRanges, paste0("/mnt/plger/jwang/data/sim/motifs/",m,"-",e,".rds"))
rm(motifRanges)
se <- getCounts(atacFrag = atacFrag,
  ranges = peaks,
  mode = "total",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "Homo sapiens",
  width = 300)

saveRDS(se, args$res)
