#args <- list(list.files("/mnt/plger/jwang/data/dat/00-frg", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(data.table)
  source("~/chromVAR/R/getCounts.R")
  library(ggplot2)
  library(ggh4x)
})
atacFrag <- readRDS(args$frg)
mice <- c("BANP", "NR1H3", "NR1H4")
m <- wcs$mot

if (m %in% mice) {
    species <- "Mus_musculus"
    genome <- BSgenome.Mmusculus.UCSC.mm10
  } else {
    species <- "Homo sapiens"
    genome <- BSgenome.Hsapiens.UCSC.hg38
}

atacFrag <- lapply(atacFrag, function(dt) {
  gr <- dtToGr(dt)
  gr <- .standardChromosomes(gr, species = species)
  as.data.table(gr)
})

nWidthBins <- 30
nGCBins <- 10
fragDt <- .getBins(atacFrag, genome = genome, 
    nWidthBins = nWidthBins, nGCBins = nGCBins)
fragDt[, bin:=paste0(widthBin, GCBin)]
fragDt[, bin:=as.integer(as.factor(bin))]
fragDt[,count_bin:=.N, by=c("sample", "bin")]
dt <- unique(fragDt, by=c("sample", "bin"))
dt$motif <- m
dt[,uncorrected_Freq:=(count_bin+1L)/(sum(count_bin)+1L),by=sample]
dt[,mean_freq:=mean(uncorrected_Freq),by=bin]
dt[,weight:=mean_freq/uncorrected_Freq]

## moderated
dt$uncorrected_moderatedFreq <- moderateBinFrequencies(dt$bin, 
    dt$sample, dt$count_bin)
dt[,moderated_mean_freq:=mean(uncorrected_moderatedFreq),by=bin]
dt[,moderated_weight:=moderated_mean_freq/uncorrected_moderatedFreq]


saveRDS(dt, args$res)
