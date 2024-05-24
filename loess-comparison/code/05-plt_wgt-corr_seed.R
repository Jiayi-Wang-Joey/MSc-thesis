#args <- list(list.files("/mnt/plger/jwang/loess-comparison/dat", full.names=TRUE), "plts/wgt-corr_seed.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(SummarizedExperiment)
    library(reshape2)
    library(tidyr)
    library(data.table)
})

fs <- args[[1]]
info <- lapply(fs, \(x) unlist(strsplit(basename(x), ","))) 
info <- data.frame(do.call(rbind, info))
colnames(info) <- c("motif", "span", "nSample", "nBin", "family", "seed")
info$seed <- gsub("\\.rds$", "", info$seed)
se <- lapply(fs, readRDS)
allCounts <- lapply(se, \(x) as.data.frame(assay(x,"counts")))

params <- expand.grid(motif=unique(info$motif),
  span=unique(info$span),
  nSample=unique(info$nSample),
  nBin=unique(info$nBin),
  family=unique(info$family))

corrs <- lapply(seq_len(nrow(params)), \(i) {
  idx <- which(info$motif==params$motif[i] &
    info$span==params$span[i] &
    info$nBin==params$nBin[i] &
    info$nSample==params$nSample[i] &
    info$family==params$family[i])
  ses <- allCounts[idx]
  samples <- colnames(ses[[1]])
  sample_corr <- lapply(samples, \(sample_id) {
    sample_count <- sapply(ses, \(x) x[,sample_id])
    corr <- cor(sample_count, use="pairwise")
    corr[lower.tri(corr, diag=TRUE)] <- NA
    na.omit(reshape2::melt(corr))$value
  })
  names(sample_corr) <- samples
  as.data.frame(sample_corr)
})

corr_motif <- lapply(unique(params$motif), \(m) {
  idx <- which(params$motif==m)
  df <- do.call(rbind, corrs[idx]) %>%
    pivot_longer(everything(),
    names_to='sample',
    values_to='value') 
})

df <- rbindlist(corr_motif, idcol="motifs")
gg <- ggplot(df,aes(x=sample,y=value,col=sample)) + 
    geom_violin() + 
    theme_minimal()


ggsave(args[[2]], gg, width=20, height=15, units="cm")



