# args <- list(list.files("/mnt/plger/jwang/scATAC-seq/waag", "-cd\\.rds$", full.names=TRUE), "plts/cd-UMAP.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggpubr)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, cell_type, mode,
    umap_1, umap_2, status, smooth, peakWeight) |>
    do.call(what=rbind) |>
    mutate(status = ifelse(status == "none", "unmoderated", status))

aes <- list(
    facet_grid(status~mode+smooth+peakWeight, scales="free", 
        labeller = \(.) label_both(.)),
    geom_point_rast(aes(umap_1, umap_2), alpha=0.5, size=0.2),
    guides(col=guide_legend(
        ncol=4, title.position="top",
        override.aes=list(alpha=1, size=1))),
    theme_minimal(6), theme(
        aspect.ratio=1,
        plot.margin=margin(),
        axis.text=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.25, "lines")))



gg <- ggplot(df, aes(col=cell_type)) + 
        aes +
        scale_color_manual(values=get_palette("Paired",
            k=length(unique(df$cell_type))))


ggsave(args[[2]], gg, width=30, height=18, units="cm")