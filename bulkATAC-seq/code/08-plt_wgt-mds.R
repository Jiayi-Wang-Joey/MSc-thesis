#args <- list(list.files("/mnt/plger/jwang/data/dat/01-weight", "^weight-.*", full.names=TRUE), "plts/wgt-corr.pdf")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(SummarizedExperiment)
  library(reshape2)
  library(edgeR)
})

fs <- args[[1]]
res <- lapply(fs, readRDS)
names(res) <- basename(fs)
mds1 <- lapply(names(res), \(x) {
  se <- res[[x]]
  info <- unlist(strsplit(basename(x), "-"))
  xy <- plotMDS.SummarizedExperiment(se, plot=FALSE)
  group_id <- rep(LETTERS[1:2],each=ncol(se)/2)
  data.frame(x=xy$x, y=xy$y, mode=info[[1]], motif=info[[2]],
    smooth=info[[3]], group_id=group_id, sample_id=colnames(se),
    peakWeight=info[[4]], status=gsub("\\.rds$", "", info[[5]]))
})
mds1 <- do.call(rbind, mds1)
mds2 <- lapply(unique(mds1$motif), \(m) {
  se <- readRDS(paste0("/mnt/plger/jwang/data/dat/01-total/total-", 
    m, ".rds"))
  xy <- plotMDS.SummarizedExperiment(se, plot=FALSE)
  group_id <- rep(LETTERS[1:2],each=ncol(se)/2)
  data.frame(x=xy$x, y=xy$y, mode="total", motif=m, 
    smooth="none", group_id=group_id, sample_id=colnames(se),
    peakWeight="none", status="unmoderated")
})
mds2 <- do.call(rbind, mds2)
mds <- rbind(mds1, mds2)
mds$sample_id <- gsub("\\.bam$|\\.bed$", "", mds$sample_id)
ps <- lapply(split(mds, mds$motif), \(df) {
  ggplot(df, aes(x,y,col=group_id,label=sample_id, shape=group_id)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(size=2) + 
    facet_grid(status~peakWeight+mode, scales="free", 
      labeller = \(.) label_both(.)) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() + theme(aspect.ratio = 1,
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey")) +
    ggtitle(df$motif[1])
})

pdf(args[[2]], onefile=TRUE, width=12, height=7)
for (p in ps) print(p); dev.off()

