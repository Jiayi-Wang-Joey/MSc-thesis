#args <- list(list.files("outs/sim", "^dif-.*", full.names=TRUE), "plts/sim/dif-heatmap.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
    library(dplyr)
})


res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, rank, t, peakWeight, mode, 
    test, name, dif, smooth, effect, adj.P.Val) |>
    do.call(what=rbind) |>
    mutate(test=ifelse(test!="ZNF143", test, "ZN143"),
        significance=ifelse(adj.P.Val < 0.05, TRUE, FALSE),
        method=paste(mode,smooth,peakWeight,dif,sep=",")) |>
    filter(name==test)

dfw <- df[df$mode=="weight",]
ps <- lapply(split(dfw, dfw$test), \(fd) {
  ggplot(fd, aes(effect, abs(t), col=peakWeight)) +
    geom_line(stat = "identity", alpha=0.8, aes(group=peakWeight)) +
    geom_point(stat = "identity", alpha=0.8, size = 1) +
    geom_label_repel(aes(label = round(t, 2))) +
    facet_grid(smooth~dif, scales = "free", labeller = \(.) label_both(.)) + 
    theme_bw() +
    ggtitle(fd$test[1]) +
    scale_color_brewer(palette = "Set1") & theme(
      plot.margin=margin(),
      panel.border=element_rect(fill=NA),
      legend.key.size=unit(0.25, "lines")) 
  
})

pdf(args[[2]], onefile=TRUE, width=8, height=8)
for (p in ps) print(p); dev.off()