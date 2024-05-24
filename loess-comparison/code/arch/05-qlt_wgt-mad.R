#args <- list(list.files("/mnt/plger/jwang/Cusanovich/loess_span", full.names=TRUE), "plts/wgt-corr_seed.pdf")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(SummarizedExperiment)
  library(reshape2)
  library(tidyr)
  library(pheatmap)
})

res <- lapply(args[[1]], readRDS)
names(res) <- gsub("\\.rds$", "",basename(args[[1]]))

.pairwiseMAD <- function(x) {
  sapply(seq_len(ncol(x)), \(i) {
    sapply(seq_len(ncol(x)), \(j) {
      median(abs(x[,i]-x[,j]))
    })
  })
}

.rename <- \(y) {
  c(sapply(seq_len(ncol(y)/2), \(x) paste0("low", x)),
    sapply(seq_len(ncol(y)/2), \(x) paste0("high", x)))
}

mads <- lapply(names(res), \(i) {
  x <- res[[i]]
  span <- unlist(strsplit(i, "_"))[[2]]
  y <- assay(x, "counts")
  colnames(y) <- .rename(y)
  mad <- .pairwiseMAD(y)
  colnames(mad) <- rownames(mad) <- colnames(y)
  mad[lower.tri(mad, diag=TRUE)] <- NA
  data.frame(na.omit(melt(mad)), span=span)
}) 

df <- do.call(rbind, mads)
df <- df %>%
  mutate(type = case_when(
    grepl("low", Var1) & grepl("low", Var2) ~ "within_low",
    grepl("high", Var1) & grepl("high", Var2) ~ "within_high",
    grepl("low", Var1) & grepl("high", Var2) ~ "between",
    grepl("high", Var1) & grepl("low", Var2) ~ "between",
    TRUE ~ NA_character_
  ))
colnames(df) <- c("sample1", "sample2", "MAD", "span", "type")


gg <- ggplot(df, aes(span, MAD, col=span)) + 
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme_bw() +
  facet_wrap(~type, scales="free") +
  scale_color_brewer(palette = "Paired")

ggsave(args[[2]], gg, width=20, height=6, units="cm")
