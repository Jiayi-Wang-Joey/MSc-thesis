#args <- list(list.files("/mnt/plger/jwang/loess-comparison/sim/01-weight", full.names=TRUE), "plts/wgt-corr_seed.pdf")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(SummarizedExperiment)
  library(reshape2)
  library(tidyr)
})

res <- lapply(args[[1]], readRDS)
names(res) <- gsub("\\.rds$", "",basename(args[[1]]))
info <- do.call(rbind, strsplit(names(res), ","))
info <- data.frame(do.call(rbind, strsplit(names(res), ",")),
  mode="weight")
colnames(info) <- c("motif", "effect", "design", "peakWeight")

.rename <- \(y) {
  c(sapply(seq_len(ncol(y)/2), \(x) paste0("ctrl", x)),
    sapply(seq_len(ncol(y)/2), \(x) paste0("treat", x)))
}

.rmDiag <- \(x, method = "spearman") {
  colnames(x) <- .rename(x)
  corr <- cor(x, use="pairwise", method = method)
  corr[lower.tri(corr, diag=TRUE)] <- NA
  na.omit(melt(corr))
}

df1 <- lapply(names(res), \(x){
  se <- res[[x]]
  motif <- unlist(strsplit(x, ","))[[1]]
  effect <- unlist(strsplit(x, ","))[[2]]
  peakWeight <- unlist(strsplit(x, ","))[[4]]
  #MAZ,FLD,group.rds
  ttl <- readRDS(paste0("/mnt/plger/jwang/loess-comparison/sim/01-total/",
    motif,",",effect,",group.rds"))
  rbind(data.frame(.rmDiag(assay(se, "counts")), 
    cor = "Spearman", mode="weight", motif=motif, 
    effect=effect, peakWeight=peakWeight),
    data.frame(.rmDiag(assay(se, "counts"), method="pearson"), 
      cor = "Pearson", mode="weight", motif=motif, 
      effect=effect, peakWeight=peakWeight),
    data.frame(.rmDiag(assay(ttl, "counts")), 
      cor = "Spearman", mode="total", motif=motif, 
      effect=effect, peakWeight="origin"),
    data.frame(.rmDiag(assay(ttl, "counts"), method="pearson"), 
      cor = "Pearson", mode="total", motif=motif, 
      effect=effect, peakWeight="origin"))
})
df <- do.call(rbind, df1)

df <- df %>%
  mutate(type = case_when(
    grepl("ctrl", Var1) & grepl("ctrl", Var2) ~ "intra",
    grepl("treat", Var1) & grepl("treat", Var2) ~ "intra",
    grepl("ctrl", Var1) & grepl("treat", Var2) ~ "inter",
    grepl("treat", Var1) & grepl("ctrl", Var2) ~ "inter",
    TRUE ~ NA_character_
  ))

ps <- lapply(split(df, df$motif), \(fd) 
  ggplot(fd, aes(x=peakWeight, y=value, col=peakWeight)) +
    geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
    geom_boxplot(width = 0.05, position = position_dodge(width = 0.8), alpha = 0.1) +
    facet_grid(cor~type+effect, scales="free") + ggtitle(fd$motif[1]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme_bw() +
    scale_color_brewer(palette = "Paired")
)

pdf(args[[2]], onefile=TRUE, width=12, height=6)
for (p in ps) print(p); dev.off()

