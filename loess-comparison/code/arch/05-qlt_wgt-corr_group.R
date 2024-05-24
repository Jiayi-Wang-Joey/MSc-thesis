#args <- list(list.files("/mnt/plger/jwang/Cusanovich/loess_span", full.names=TRUE), "plts/wgt-corr_seed.pdf")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(SummarizedExperiment)
  library(reshape2)
  library(tidyr)
})

res <- lapply(args[[1]], readRDS)
names(res) <- gsub("\\.rds$", "",basename(args[[1]]))

.rename <- \(y) {
  c(sapply(seq_len(ncol(y)/2), \(x) paste0("low", x)),
    sapply(seq_len(ncol(y)/2), \(x) paste0("high", x)))
}

.rmDiag <- \(x, method = "spearman") {
  colnames(x) <- .rename(x)
  corr <- cor(x, use="pairwise", method = method)
  corr[lower.tri(corr, diag=TRUE)] <- NA
  na.omit(melt(corr))
}

df1 <- lapply(names(res), \(x){
  se <- res[[x]]
  span <- as.numeric(unlist(strsplit(x, "_"))[[2]])
  rbind(data.frame(.rmDiag(assay(se, "counts")), 
    cor = "Spearman", span=span),
    data.frame(.rmDiag(assay(se, "counts"), method="pearson"), 
      cor = "Pearson", span=span))
})
df <- do.call(rbind, df1)
df$span <- as.factor(df$span)

df <- df %>%
  mutate(type = case_when(
    grepl("low", Var1) & grepl("low", Var2) ~ "within_low",
    grepl("high", Var1) & grepl("high", Var2) ~ "within_high",
    grepl("low", Var1) & grepl("high", Var2) ~ "between",
    grepl("high", Var1) & grepl("low", Var2) ~ "between",
    TRUE ~ NA_character_
  ))

gg <- ggplot(df, aes(span, value, col=span)) + 
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme_bw() +
  facet_wrap(type~cor, scales="free") +
  scale_color_brewer(palette = "Paired")

ggsave(args[[2]], gg, width=25, height=12, units="cm")

