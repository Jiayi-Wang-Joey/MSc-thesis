# args <- args <- list(list.files("/mnt/plger/jwang/scATAC-seq/02-fil_wgt/", "-cd\\.rds$", full.names=TRUE), "plts/cd-UMAP.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggrastr)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, sample_id, group_id, 
    umap_1, umap_2, motif, mode, status, smooth, peakWeight) |>
    do.call(what=rbind) |>
    mutate(status = ifelse(status == "none", "unmoderated", status))

aes <- list(
    facet_grid(status~mode+smooth+peakWeight, scales="free", 
        labeller = \(.) label_both(.)),
    geom_point_rast(aes(umap_1, umap_2), alpha=0.5, size=1),
    guides(col=guide_legend(
        ncol=4, title.position="top",
        override.aes=list(alpha=1, size=1))),
    theme_bw(), theme(
        aspect.ratio=1,
        plot.margin=margin(),
        axis.text=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.25, "lines")))


ps <- lapply(split(df,df$motif), \(fd){
    ggplot(fd, aes(col=group_id, shape=sample_id)) + 
        aes +
        scale_color_manual(values=c("royalblue", "tomato"))
})

pdf(args[[2]], onefile=TRUE, width=15, height=10)
for (p in ps) print(p); dev.off()

