#args <- list(list.files("/mnt/plger/jwang/data/dat/01-weight", "^weight-.*", full.names=TRUE), "plts/wgt-corr.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(SummarizedExperiment)
    library(reshape2)
    library(ggh4x)
})

# read data
info <- lapply(args[[1]], \(x) unlist(strsplit(basename(x), "-"))) 
info <- data.frame(do.call(rbind, info))
colnames(info) <- c("mode", "motif", "smooth", "peakWeight", "moderated")
info$moderated <- gsub("\\.rds$", "", info$moderated)
se <- lapply(args[[1]], readRDS)


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

df <- lapply(seq_len(length(se)), \(i) {
  x <- se[[i]]
  m <- info[i,"motif"]
  total <- readRDS(paste0("/mnt/plger/jwang/data/dat/01-total/total-", m, ".rds"))
  z <- assay(total, "counts")
  y <- assay(x, "counts")
  rbind(data.frame(.rmDiag(y, method="pearson"),
    cor = "Pearson", peakWeight = info[i, "peakWeight"], 
    motif = m, status=info[i,"moderated"]),
    data.frame(.rmDiag(y, method="spearman"), 
      cor = "Spearman", peakWeight = info[i, "peakWeight"], 
      motif = m, status=info[i,"moderated"]),
    data.frame(.rmDiag(z, method="pearson"),
      cor = "Pearson", peakWeight = "origin", 
      motif = m, status=info[i,"moderated"]),
    data.frame(.rmDiag(z, method="spearman"),
      cor = "Spearman", peakWeight = "origin", 
      motif = m, status=info[i,"moderated"]))
}) %>% bind_rows()

df <- df %>%
  mutate(type = case_when(
    grepl("ctrl", Var1) & grepl("ctrl", Var2) ~ "intra",
    grepl("treat", Var1) & grepl("treat", Var2) ~ "intra",
    grepl("ctrl", Var1) & grepl("treat", Var2) ~ "inter",
    grepl("treat", Var1) & grepl("ctrl", Var2) ~ "inter",
    TRUE ~ NA_character_
  ))

gg <- lapply(split(df, df$motif), \(fd) 
  ggplot(fd, aes(x=peakWeight, y=value, col=peakWeight)) +
      geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
      geom_boxplot(width = 0.05, position = position_dodge(width = 0.8), alpha = 0.1) +
      facet_grid2(cor~type+status, scales="free", independent = "y") + ggtitle(fd$motif[1]) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      scale_color_brewer(palette = "Paired") +
      theme_bw()
) |> wrap_plots(n=2) + plot_layout(guides = "collect")

ggsave(args[[2]], gg, width=60, height=40, units="cm")




