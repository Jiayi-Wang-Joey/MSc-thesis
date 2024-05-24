#args <- list(list.files("/mnt/plger/jwang/data/sim/01-weight", "^weight-.*", full.names=TRUE), "plts/wgt-corr.pdf")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(SummarizedExperiment)
  library(reshape2)
})

# read data

fs <- args[[1]]
info <- lapply(fs, \(x) unlist(strsplit(basename(x), "-"))) 
info <- data.frame(do.call(rbind, info))
colnames(info) <- c("mode", "motif", "effect", "smooth", "peakWeight")
info$peakWeight <- gsub("\\.rds$", "", info$peakWeight)
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
  me <- paste0(info[i,"motif"], "-", info[i,"effect"])
  total <- readRDS(paste0("/mnt/plger/jwang/data/sim/01-total/total-", 
    me, ".rds"))
  z <- assay(total, "counts")
  y <- assay(x, "counts")
  rbind(data.frame(.rmDiag(y, method="pearson"),
    cor = "Pearson", peakWeight = info[i, "peakWeight"], 
    motif = info[i,"motif"], effect=info[i,"effect"]),
    data.frame(.rmDiag(y, method="spearman"), 
      cor = "Spearman", peakWeight = info[i, "peakWeight"], 
      motif = info[i,"motif"], effect=info[i,"effect"]),
    data.frame(.rmDiag(z, method="pearson"),
      cor = "Pearson", peakWeight = "origin", 
      motif = info[i,"motif"], effect=info[i,"effect"]),
    data.frame(.rmDiag(z, method="spearman"),
      cor = "Spearman", peakWeight = "origin", 
      motif = info[i,"motif"], effect=info[i,"effect"]))
}) %>% bind_rows()

df <- df %>%
  mutate(type = case_when(
    grepl("ctrl", Var1) & grepl("ctrl", Var2) ~ "intra",
    grepl("treat", Var1) & grepl("treat", Var2) ~ "intra",
    grepl("ctrl", Var1) & grepl("treat", Var2) ~ "inter",
    grepl("treat", Var1) & grepl("ctrl", Var2) ~ "inter",
    TRUE ~ NA_character_
  ))

ps <- lapply(split(df, df$motif), \(fd) {
  lapply(split(fd, fd$effect), \(d) 
    ggplot(d, aes(x=peakWeight, y=value, col=peakWeight)) +
      geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
      geom_boxplot(width = 0.05, 
        position = position_dodge(width = 0.8), alpha = 0.1) +
      facet_grid(cor~type, scales="free") + ggtitle(d$effect[1]) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        scale_color_brewer(palette = "Paired") +
        theme_bw()
  ) |> wrap_plots(n=2) + plot_layout(guides = "collect") + plot_annotation(fd$motif[1]) 
})


pdf(args[[2]], onefile=TRUE, width=15, height=12)
for (p in ps) print(p); dev.off()